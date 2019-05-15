#This script calculates short versus long insert ratio (Number of inserts <=150bp)/ (Number of inserts > 150bp) in ATAC-seq samples. 
#It's part of a picard based pipeline
#We have experimentaly established a 1.5 threshold for this ratio.
##Lara Bossini-Castillo 8th Nov 2017 lbc@sanger.ac.uk

temp = list.files(pattern="*_insert_size_txt_2")
list2env(
  lapply(setNames(temp, make.names(gsub("\\..*", "_Ratio", temp))), 
         function(i){read.table(i, header = TRUE)}), envir = .GlobalEnv) 
all_samples <- ls(pattern="_Ratio")
sink('analysis-output.txt')
for (i in 1:length(all_samples)){
  sample = all_samples[i]
  reads <- get(sample)
  short_reads <- reads[reads$insert_size <= 150, ]
  long_reads <- reads[reads$insert_size > 150, ]
  ratio <- sum(short_reads$All_Reads.fr_count)/sum(long_reads$All_Reads.fr_count)
  cat(paste("The short read-long read ratio for ",sample," = ",ratio,"\n",sep=""))
}
sink()
