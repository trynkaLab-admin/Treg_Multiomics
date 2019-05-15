library("dplyr")
library("tidyr")
library("purrr")
library("devtools")
library("SummarizedExperiment")

#Import GWAS traits and only keep immune diseases
gwas_stats_labeled = readr::read_tsv("analysis/data/gwas/GWAS_summary_stat_list.labeled.txt",
                                     col_names = c("trait","file_name", "type")) %>%
  dplyr::filter(!(trait %in% c("MI","PLT","CRE","CRP","GLU","SSc","HGB","NAR","PCV","HOMA_b","HOMA_ir",
                                "TC","TG","UA","WBC","HDL","IL6","INS","MCV","RBC")))

#Import RNA SummarizedExperiments
se_featureCounts = readRDS("results/SummarizedExperiments/treg_featureCounts.rds")

#Identify genes in the MHC region that should be excluded
mhc_featureCounts = dplyr::filter(tbl_df2(rowData(se_featureCounts)), chr == "6", start > 28510120, end < 33480577) %>%
  dplyr::rename(phenotype_id = gene_id)

#Gene names
featureCounts_name_map = dplyr::select(tbl_df2(rowData(se_featureCounts)), gene_id, gene_name) %>% 
  dplyr::rename(phenotype_id = gene_id)

#Import RNA coloc output
file_names = as.list(paste0("processed/coloc_500/", gwas_stats_labeled$trait, ".RNA.2e5.txt"))
names(file_names) = gwas_stats_labeled$trait
  
#Import RNA coloc enrichments
coloc_df = purrr::map_df(file_names, ~readr::read_delim(., delim = "\t"), .id = "trait") %>%
  dplyr::mutate(PP_power = (PP.H4.abf + PP.H3.abf), PP_coloc = PP.H4.abf/PP_power)

coloc_hits = dplyr::filter(coloc_df, PP_power > 0.8) %>% 
  dplyr::filter(PP.H4.abf > 0.8) %>%
  dplyr::group_by(trait, gwas_lead) %>% 
  dplyr::arrange(trait, gwas_lead, -PP.H4.abf) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::filter(PP_coloc > 0.8) %>%
  dplyr::filter(nsnps > 50) %>%
  dplyr::filter(!(phenotype_id %in% mhc_phenotypes$phenotype_id)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(featureCounts_name_map, by = "phenotype_id")

write.table(featureCounts_200kb_hits, "results/coloc/treg_eQTL_coloc_hits.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Import H3K27ac coloc output
file_names = as.list(paste0("processed/coloc_500/", gwas_stats_labeled$trait, ".K27ac.2e5.txt"))
names(file_names) = gwas_stats_labeled$trait
  
#Import H3K27ac coloc enrichments
coloc_df = purrr::map_df(file_names, ~readr::read_delim(., delim = "\t"), .id = "trait") %>%
  dplyr::mutate(PP_power = (PP.H4.abf + PP.H3.abf), PP_coloc = PP.H4.abf/PP_power)

coloc_hits = dplyr::filter(coloc_df, PP_power > 0.8) %>% 
  dplyr::filter(PP.H4.abf > 0.8) %>%
  dplyr::group_by(trait, gwas_lead) %>% 
  dplyr::arrange(trait, gwas_lead, -PP.H4.abf) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::filter(PP_coloc > 0.8) %>%
  dplyr::filter(nsnps > 50) %>%
  dplyr::ungroup()

write.table(coloc_hits, "K27ac_treg_GWAS_coloc_hits.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Import H3K4me3 coloc output
file_names = as.list(paste0("processed/coloc_500/", gwas_stats_labeled$trait, ".K4me3.2e5.txt"))
names(file_names) = gwas_stats_labeled$trait
  
#Import H3K4me3 coloc enrichments
coloc_df = purrr::map_df(file_names, ~readr::read_delim(., delim = "\t"), .id = "trait") %>%
  dplyr::mutate(PP_power = (PP.H4.abf + PP.H3.abf), PP_coloc = PP.H4.abf/PP_power)

coloc_hits = dplyr::filter(coloc_df, PP_power > 0.8) %>% 
  dplyr::filter(PP.H4.abf > 0.8) %>%
  dplyr::group_by(trait, gwas_lead) %>% 
  dplyr::arrange(trait, gwas_lead, -PP.H4.abf) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::filter(PP_coloc > 0.8) %>%
  dplyr::filter(nsnps > 50) %>%
  dplyr::ungroup()

write.table(coloc_hits, "K4me3_treg_GWAS_coloc_hits.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Import ATAC coloc output
file_names = as.list(paste0("processed/coloc_500/", gwas_stats_labeled$trait, ".ATAC.2e5.txt"))
names(file_names) = gwas_stats_labeled$trait

#Import ATAC coloc enrichments
coloc_df = purrr::map_df(file_names, ~readr::read_delim(., delim = "\t"), .id = "trait") %>%
  dplyr::mutate(PP_power = (PP.H4.abf + PP.H3.abf), PP_coloc = PP.H4.abf/PP_power)

coloc_hits = dplyr::filter(coloc_df, PP_power > 0.8) %>%
  dplyr::filter(PP.H4.abf > 0.8) %>%
  dplyr::group_by(trait, gwas_lead) %>% 
  dplyr::arrange(trait, gwas_lead, -PP.H4.abf) %>% 
  dplyr::filter(row_number() == 1) %>% 
  dplyr::filter(PP_coloc > 0.8) %>%
  dplyr::filter(nsnps > 50) %>%
  dplyr::ungroup()

write.table(coloc_hits, "ATAC_treg_GWAS_coloc_hits.txt", sep = "\t", quote = FALSE, row.names = FALSE)
