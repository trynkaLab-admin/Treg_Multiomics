library("dplyr")
library("devtools")
library("SummarizedExperiment")
library("purrr")
load_all("seqUtils")

data_info <- read.table("analysis/data/treg_design_matrix.txt", header=TRUE)

#Export featureCounts cqns for qtl mapping
se_gene = readRDS("results/SummarizedExperiments/treg_featureCounts.rds")
event_metadata = rowData(se_gene)

event_metadata_rm = event_metadata[event_metadata$chr == "6" & event_metadata$start > 20000000 & event_metadata$end < 40000000,]
event_metadata = event_metadata[!event_metadata$gene_id %in% event_metadata_rm$gene_id,]

event_dataset = se_gene[event_metadata$gene_id,]
event_dataset = se_gene[event_metadata[event_metadata$chr %in% c(1:22),]$gene_id,]

#Keep only expressed genes
counts_matrix = assays(event_dataset)$counts
mean_across = rowMeans(counts_matrix[,c(35:134)])
mean_across <- mean_across[mean_across > 25]

expressed_dataset = event_dataset[rownames(event_dataset) %in% names(mean_across),]

#Extract lists for each condition
on_list = idVectorToList(c("day0","rest","PMA"))
first <- expressed_dataset[,expressed_dataset$process == "day0"]
second <- expressed_dataset[,expressed_dataset$process == "rest"]
third <- expressed_dataset[,expressed_dataset$process == "PMA"]
event_conditions <- list(first,second,third)
names(event_conditions) <- c("day0","rest","PMA")

#Construct gene positions for QTL mapping
fastqtl_genepos = tbl_df2(rowData(expressed_dataset)) %>%
  dplyr::filter(gene_id %in% rownames(expressed_dataset)) %>%
  dplyr::mutate(transcript_id = gene_id) %>%
  constructQTLtoolsGenePos()
output_path = "processed/qtltools/input/featureCounts/"

#Extract cqn matrices
cqn_list = purrr::map(event_conditions, ~assays(.)$cqn)

colnames(cqn_list[[1]]) <- data_info$donor[match(colnames(cqn_list[[1]]),data_info$sample_id)]
colnames(cqn_list[[2]]) <- data_info$donor[match(colnames(cqn_list[[2]]),data_info$sample_id)]
colnames(cqn_list[[3]]) <- data_info$donor[match(colnames(cqn_list[[3]]),data_info$sample_id)]

qtltools_cqn_list = purrr::map(cqn_list, prepareQTLtoolsMatrix, fastqtl_genepos)
saveFastqtlMatrices(qtltools_cqn_list, output_path, file_suffix = "norm_prop")

#Calculate covariates
sample_meta = tbl_df2(colData(event_conditions$day0))
pca = prcomp(t(cqn_list[[1]]))
pca_mat = as.data.frame(pca$x[,1:13])
pca_mat$sample_id <- rownames(pca_mat)
covariates <- fastqtlMetadataToCovariates(pca_mat)

file_path = file.path(output_path, paste("day0","covariates_prop", "txt", sep = "."))
write.table(covariates, file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
