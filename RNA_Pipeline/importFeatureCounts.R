library("dplyr")
library("readr")
library("devtools")
library("SummarizedExperiment")
library("cqn")

load_all("seqUtils")

#Import sample names
design_matrix = read.table("analysis/data/treg_design_matrix.txt", sep ="\t",comment.char = "", stringsAsFactors = FALSE, header=TRUE)

#Remove bad samples
design_matrix = subset(design_matrix, cell_type != "Tcon") #Remove Tcons
design_matrix = subset(design_matrix, sample_id != "I0032") #Not stimulated
design_matrix = subset(design_matrix, sample_id != "I0034") #Not stimulated
design_matrix = subset(design_matrix, sample_id != "I0038") #Not stimulated
design_matrix = subset(design_matrix, sample_id != "I0041") #Not stimulated
design_matrix = subset(design_matrix, sample_id != "I0043") #Not stimulated
design_matrix = subset(design_matrix, sample_id != "I0046") #Not stimulated
design_matrix = subset(design_matrix, sample_id != "I0049") #Not stimulated
design_matrix = subset(design_matrix, sample_id != "I0356") #Low quality
design_matrix = subset(design_matrix, sample_id != "I0304") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0273") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0288") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0367") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0279") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0316") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0313") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0271") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0348") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0285") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0346") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0269") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0324") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0362") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0292") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0263") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0282") #Replicate
design_matrix = subset(design_matrix, sample_id != "I1408") #Replicate
design_matrix = subset(design_matrix, sample_id != "I0370") #Kinship

#Import featureCounts
count_matrix = loadCounts("processed/featureCounts/", design_matrix$sample_id, sub_dir = FALSE, counts_suffix = ".featureCounts.txt")

#Import transcript metadata
transcript_data = tbl_df(readRDS("processed/annotations/Homo_sapiens.GRCh38.87.transcript_data.rds")) %>%
  dplyr::rename(gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id, gene_name = external_gene_name, chr = chromosome_name)

#Filter transcript metadata
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9","MT","X","Y")
valid_gene_biotypes = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                        "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                        "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                        "antisense","sense_intronic","sense_overlapping")
filtered_tx_data = dplyr::filter(transcript_data, gene_biotype %in% valid_gene_biotypes, chr %in% valid_chromosomes)

#Extract length
length_df = dplyr::select(count_matrix, gene_id, length)

#Construct transcript metadata
gene_metadata = filtered_tx_data %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(start = min(transcript_start), end = max(transcript_end)) %>%
  dplyr::select(gene_id, gene_biotype, chr, start, end, gene_name, strand, percentage_gc_content) %>% 
  dplyr::left_join(length_df, by = "gene_id") %>%
  unique() %>%
  dplyr::ungroup() %>%
  as.data.frame()
rownames(gene_metadata) = gene_metadata$gene_id

#Process counts
filtered_data = dplyr::filter(count_matrix, gene_id %in% gene_metadata$gene_id)
counts = dplyr::select(filtered_data, -gene_id, -length)
rownames(counts) = filtered_data$gene_id
counts = counts[gene_metadata$gene_id,] #Reoder counts

#CQN normalize counts
cqn_matrix = calculateCQN(counts, gene_metadata)

#Prepare sample metadata
sample_metadata = readRDS("analysis/data/compiled_treg_metadata.rds")
sample_meta = dplyr::left_join(design_matrix, sample_metadata, by = c("donor","process")) %>% as.data.frame()
rownames(sample_meta) = sample_meta$sample_id

#Construct a SummarizedExperiment object
se = SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = counts, cqn = cqn_matrix), 
  colData = sample_meta, 
  rowData = gene_metadata)
saveRDS(se, "results/SummarizedExperiments/treg_featureCounts.rds")
