df <- read.table("RNA_nominal_results.txt")
df <- na.omit(df)
df$V15 <- p.adjust(df$V12, method="fdr")
write.table(df, file="RNA_nominal_results_FDRcor.txt", sep="\t", row.names = FALSE, quote = FALSE)
