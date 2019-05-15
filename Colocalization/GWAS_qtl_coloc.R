library("dplyr")
library("tibble")
library("tidyr")
library("purrr")
library("coloc")
library("readr")
library("devtools")
library("optparse")
load_all("seqUtils")

#Parse command-line options
option_list <- list(
  make_option(c("-p", "--phenotype"), type="character", default=NULL,
              help="Type of QTLs used for coloc.", metavar = "type"),
  make_option(c("--qtl"), type="character", default=NULL,
              help="Path to the QTL directory.", metavar = "type"),
  make_option(c("-w", "--window"), type="character", default=NULL,
              help="Size of the cis window.", metavar = "type"),
  make_option(c("--gwas"), type="character", default=NULL,
              help="Name of the GWAS trait", metavar = "type"),
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Path to GWAS summary stats directory.", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Path to the output directory.", metavar = "type"),
  make_option(c("-s", "--samplesizes"), type="character", default=NULL,
              help="Path to the tab-separated text file with condition names and sample sizes.", metavar = "type"),
  make_option(c("--gwaslist"), type="character", default=NULL,
              help="Path to the list of GWAS studies.", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Extract parameters for CMD options
gwas_id = opt$gwas
cis_window = as.numeric(opt$w)
phenotype = opt$p
gwas_dir = opt$d
outdir = opt$o
sample_size_path = opt$s
gwas_list = opt$gwaslist
qtl_dir = opt$qtl

#Import variant information
GRCh38_variants = importVariantInformation("ref_files/AllDonors_Imputed_AllAutosomes_MAF0.1_GRCh38_b38.txt")
GRCh37_variants = importVariantInformation("ref_files/AllDonors_Imputed_AllAutosomes_MAF0.1_GRCh38_b37.txt")

#Import list of GWAS studies
gwas_stats_labeled = readr::read_tsv(gwas_list, col_names = c("trait","file_name","type","ratio_cases"), col_type = "cccn")

#Import sample sizes
sample_sizes = readr::read_tsv(sample_size_path, col_names = c("condition_name", "sample_size"), col_types = "cc")
sample_sizes_list = as.list(sample_sizes$sample_size)
names(sample_sizes_list) = sample_sizes$condition_name

#Construct a new QTL list 
phenotype_values = constructQtlListForColoc(phenotype, "processed/qtltools/output", sample_sizes_list)

#Spcecify the location of the GWAS summary stats file
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
gwas_prefix = file.path(gwas_dir, gwas_file_name)
ratio = gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$ratio_cases

#Prefilter coloc candidates
qtl_df_list = prefilterColocCandidates(phenotype_values$min_pvalues, gwas_prefix, 
                                       GRCh37_variants, fdr_thresh = 0.05, 
                                       overlap_dist = 1e5, gwas_thresh = 1e-5)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()

dim(qtl_pairs)

#Test for coloc
y <- data.frame()

for (i in 1:dim(qtl_pairs)[1]) {
   tryCatch({
    my_temp_df <- qtl_pairs[i,]
    qtl_ranges = constructVariantRanges(my_temp_df, GRCh38_variants, cis_dist = cis_window)
    gwas_ranges = constructVariantRanges(my_temp_df, GRCh37_variants, cis_dist = cis_window)
    
    qtl_summaries = qtltoolsTabixFetchPhenotypes(qtl_ranges, phenotype_values$qtl_summary_list$day0)[[1]] %>%
        dplyr::transmute(snp_id, chr = snp_chr, pos = snp_start, p_nominal, beta) 
    gwas_summaries = tabixFetchGWASSummary(gwas_ranges, paste0(gwas_prefix, ".sorted.txt.gz"))[[1]]
    
    qtl = summaryReplaceCoordinates(qtl_summaries, GRCh37_variants)
    gwas = summaryReplaceSnpId(gwas_summaries, GRCh37_variants)

    qtl_min = dplyr::arrange(qtl, p_nominal) %>% dplyr::filter(row_number() == 1)
    gwas_min = dplyr::arrange(gwas, p_nominal) %>% dplyr::filter(row_number() == 1)

    #Perform coloc analysis
    coloc_res = colocQtlGWAS(qtl, gwas, N_qtl = 123,ratio_case=ratio)
    coloc_summary = dplyr::tbl_df(t(data.frame(coloc_res$summary))) %>%
      dplyr::mutate(qtl_pval = qtl_min$p_nominal, gwas_pval = gwas_min$p_nominal,
                    qtl_lead = qtl_min$snp_id, gwas_lead = gwas_min$snp_id, phenotype_id = my_temp_df$phenotype_id) #Add minimal pvalues
    y <- rbind.data.frame(y, coloc_summary)
   }, error=function(e){})  
}

#Export results
coloc_output = file.path(outdir, paste(gwas_id, phenotype, opt$w, "txt", sep = "."))
write.table(y, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)
