library("dplyr")
library("tidyr")
library("purrr")
library("devtools")
library("SummarizedExperiment")

importAndFilterColocHits <- function(gwas_stats, coloc_suffix = ".eQTL.2e5.coloc.txt", coloc_prefix = "/coloc/",
                                     PP_power_thresh = 0.8, PP_coloc_thresh = .9, nsnps_thresh = 50, gwas_pval_thresh = 1e-6, mhc_phenotypes){
  #Import coloc hits
  
  #Name all files
  file_names = as.list(paste0(coloc_prefix, gwas_stats_labeled$trait, coloc_suffix))
  names(file_names) = gwas_stats_labeled$trait
  
  #Import enrichments
  coloc_df = purrr::map_df(file_names, ~readr::read_delim(., delim = "\t"), .id = "trait") %>%
    dplyr::mutate(PP_power = (PP.H4.abf + PP.H3.abf), PP_coloc = PP.H4.abf/PP_power)

  #Identify one overlap GWAS lead varaint
  coloc_hits = identifyColocHits(coloc_df, PP_power_thresh, PP_coloc_thresh, nsnps_thresh) %>%
    dplyr::filter(gwas_pval < 1e-6) %>%
    dplyr::filter(!(phenotype_id %in% mhc_phenotypes$phenotype_id))
  
  #Merge IBD overlaps together, keep only stronges association per gene
  coloc_hits_merged = dplyr::group_by(coloc_hits, summarised_trait, phenotype_id) %>% 
    dplyr::arrange(summarised_trait, phenotype_id, -PP.H4.abf) %>% 
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup()
  
  #Retreive posterior probabilitites in all conditions
  coloc_filtered = dplyr::semi_join(coloc_df, coloc_hits_merged, by = c("trait", "phenotype_id"))
  return(list(coloc_filtered = coloc_filtered, coloc_df = coloc_df))
}


#Import SummarizedExperiments
se_featureCounts = readRDS("results/SummarizedExperiments/treg_featureCounts.rds")

#Identify genes in the MHC region that should be excluded
mhc_featureCounts = dplyr::filter(tbl_df2(rowData(se_featureCounts)), chr == "6", start > 28510120, end < 33480577) %>%
  dplyr::rename(phenotype_id = gene_id)

#Gene names
featureCounts_name_map = dplyr::select(tbl_df2(rowData(se_featureCounts)), gene_id, gene_name) %>% 
  dplyr::rename(phenotype_id = gene_id)

#Import GWAS traits and only keep immune diseases
gwas_stats_labeled = readr::read_tsv("analysis/data/gwas/GWAS_summary_stat_list.labeled.txt",
                                     col_names = c("trait","file_name", "type")) %>%
  dplyr::filter(!(trait %in% c("MI","PLT","CRE","CRP","GLU","SSc","HGB","NAR","PCV","HOMA_b","HOMA_ir",
                                "TC","TG","UA","WBC","HDL","IL6","INS","MCV","RBC")))

featureCounts_200kb_hits = importAndFilterColocHits(gwas_stats_labeled, coloc_suffix = ".featureCounts.2e5.txt", 
                                          coloc_prefix = "processed/coloc/",
                                          PP_power_thresh = 0.8, PP_coloc_thresh = .8, nsnps_thresh = 50, 
                                          gwas_pval_thresh = 1e-5, mhc_phenotypes = mhc_tpm)$coloc_filtered %>%
  dplyr::left_join(tpm_name_map, by = "phenotype_id")

#Save as an RDS and a text file
saveRDS(featureCounts_200kb_hits, "results/coloc/treg_eQTL_coloc_hits.rds")
write.table(featureCounts_200kb_hits, "results/coloc/treg_eQTL_coloc_hits.txt", sep = "\t", quote = FALSE, row.names = FALSE)
