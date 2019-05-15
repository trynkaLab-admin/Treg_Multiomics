library(UpSetR)
library(Hmisc)

df <- read.table("Summary_table_Test_LeadSlope_NoAllNAs_New.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
#Remove duplicated values in same cell
for (i in 1:length(df)){
  df[,i] <- as.character(df[,i])
  df[,i] <- sapply(strsplit(df[,i], ",", fixed = TRUE), function(x) 
    paste(unique(x), collapse = ","))
}

#Change NAs to zeros
df[df=="NA"] <- NA
df[] <- lapply(df, function(x) {
  is.na(levels(x)) <- levels(x) == "NA"
  x
})
df[is.na(df)] <- 0

###Write to file
write.table(df, file="Summary_table_Test_LeadSlope_NoAllNAs_Processed_New.txt", sep="\t", quote = FALSE, row.names = FALSE)

#This step is no longer necessary
#Remove SNPs with different position in tabix and in VCF (real minority, but enough to cause troubles)
#pos =c()
#toremove = c()

#for (i in 1:length(df$eQTL_Gene_featureCounts)){
#  x <- as.vector(df$eQTL_Gene_featureCounts[i])
#  x = unlist(strsplit(x, split=","))
#  z = as.vector(df$eQTL_Gene_featureCounts[i])
#  z = unlist(strsplit(z, split="_"))
#  z = unlist(strsplit(z, split=","))
#  rm(pos)
#  pos = c()
#  j = 1
#  if(z[j] != 0){
#    pos[j] <- match(z[j],x, nomatch = 9999)
#    y <- as.vector(df$gene_slope_featureCounts[i])
#    y = unlist(strsplit(y, split=","))
#    if (pos[j] == 9999){ 
#      toremove <- c(toremove,i)
#    }
#  }
#}

#df <- df[-c(toremove),]
#row.names(df) <- 1:nrow(df)
##Loading preprocessed table

df <- read.table(file="Summary_table_Test_LeadSlope_NoAllNAs_Processed_New.txt", sep="\t", header = TRUE)
df_test <- df[which(df$QTL_ATAC_Peak != 0 & df$eQTL_Gene_featureCounts !=0),]
atac_inQTLpeak <- df_test[which(is.element(df_test$ATAC_Peak_Overlap,df_test$QTL_ATAC_Peak)),]
atac_inQTLpeak <- atac_inQTLpeak$ID
df_test <- df[which(df$QTL_K27ac_Peak != 0 & df$eQTL_Gene_featureCounts !=0),]
k27ac_inQTLpeak <- df_test[which(is.element(df_test$K27ac_Peak_Overlap,df_test$QTL_K27ac_Peak)),]
k27ac_inQTLpeak <- k27ac_inQTLpeak$ID
df_test <- df[which(df$QTL_K4me3_Peak != 0 & df$eQTL_Gene_featureCounts !=0),]
k4me3_inQTLpeak <- df_test[which(is.element(df_test$K4me3_Peak_Overlap,df_test$QTL_K4me3_Peak)),]
k4me3_inQTLpeak <- k4me3_inQTLpeak$ID
atac_inQTLpeak <-as.data.frame(atac_inQTLpeak)
names(atac_inQTLpeak) <- c("V1")
k27ac_inQTLpeak <- as.data.frame(k27ac_inQTLpeak)
names(k27ac_inQTLpeak) <- c("V1")
k4me3_inQTLpeak <- as.data.frame(k4me3_inQTLpeak)
names(k4me3_inQTLpeak) <- c("V1")
length(atac_inQTLpeak$V1)
length(k27ac_inQTLpeak$V1)
length(k4me3_inQTLpeak$V1)
eQTL_epiQTL_inSamePeak <- rbind(atac_inQTLpeak,k27ac_inQTLpeak,k4me3_inQTLpeak)
eQTL_epiQTL_inSamePeak <- unique(eQTL_epiQTL_inSamePeak$V1)

###Finding out if slopes of QTLs have same direction

df$QTL_SlopeScore <- 0
df$QTL_SlopeScore2 <- 0
###4 QTLs same direction
categories <- c(7, 8, 9, 14, 20, 26, 32)

for (i in categories){
  df[,i] <- as.character(df[,i])
  df[,i]<- ifelse(df[,i] == 0 ,0 ,1)
}

for (i in 1:length(df$eQTL_Gene_featureCounts)){
  if(as.numeric(as.character(df$gene_nominal_pvalue_featureCounts[i])) + as.numeric(as.character(df$ATAC_peak_nominal_pvalue[i])) + as.numeric(as.character(df$K27ac_peak_nominal_pvalue[i])) + as.numeric(as.character(df$K4me3_peak_nominal_pvalue[i])) == 4){
    x <- as.vector(df$gene_slope_featureCounts[i])
    x = unlist(strsplit(x, split=","))
    a <- as.vector(df$ATAC_peak_slope[i])
    a = unlist(strsplit(a, split=","))  
    b <- as.vector(df$K27ac_peak_slope[i])
    b = unlist(strsplit(b, split=","))  
    c <- as.vector(df$K4me3_peak_slope[i])
    c = unlist(strsplit(c, split=","))
    QTL1 = c()
    QTL2 = c()
    QTL3 = c()
    QTL4 = c()
    for (xx in 1:length(x)){
      if (x[xx] > 0){QTL1 <- c(QTL1,x[xx])}}
    for (aa in 1: length(a)){
      if(a[aa] > 0){QTL2 <- c(QTL2,a[aa])}} 
    for (bb in 1: length(b)){
      if(b[bb] > 0){QTL3 <- c(QTL3,b[bb])}}
    for (cc in 1: length(c)){
      if(c[cc] > 0){QTL4 <- c(QTL4,c[cc])}}
    if (length(QTL1) > 0){
      QTL1 <- NULL
      QTL1 = 1
    }else{
      QTL1 <- NULL
      QTL1 = 0
    }
    if (length(QTL2) > 0){
      QTL2 <- NULL
      QTL2 = 1
    }else{
      QTL2 <- NULL
      QTL2 = 0      
    }
    if (length(QTL3) > 0){
      QTL3 <- NULL
      QTL3 = 1
    }else{
      QTL3 <- NULL
      QTL3 = 0      
    }
    if (length(QTL4) > 0){
      QTL4 <- NULL
      QTL4 = 1
    }else{
      QTL4 <- NULL
      QTL4 = 0
    }
    QTL1_2 = c()
    QTL2_2 = c()
    QTL3_2 = c()
    QTL4_2 = c()
    for (xx in 1:length(x)){
      if (x[xx] < 0){QTL1_2 <- c(QTL1_2,x[xx])}}
    for (aa in 1: length(a)){
      if(a[aa] < 0){QTL2_2 <- c(QTL2_2,a[aa])}} 
    for (bb in 1: length(b)){
      if(b[bb] < 0){QTL3_2 <- c(QTL3_2,b[bb])}}
    for (cc in 1: length(c)){
      if(c[cc] < 0){QTL4_2 <- c(QTL4_2,c[cc])}}
    if (length(QTL1_2) > 0){
      QTL1_2 <- NULL
      QTL1_2 = 1
    }else{
      QTL1_2 <- NULL
      QTL1_2 = 0
    }
    if (length(QTL2_2) > 0){
      QTL2_2 <- NULL
      QTL2_2 = 1
    }else{
      QTL2_2 <- NULL
      QTL2_2 = 0      
    }
    if (length(QTL3_2) > 0){
      QTL3_2 <- NULL
      QTL3_2 = 1
    }else{
      QTL3_2 <- NULL
      QTL3_2 = 0      
    }
    if (length(QTL4_2) > 0){
      QTL4_2 <- NULL
      QTL4_2 = 1
    }else{
      QTL4_2 <- NULL
      QTL4_2 = 0
    }
    if(QTL1 + QTL2 + QTL3 + QTL4 ==4){
      df$QTL_SlopeScore[i] <- paste(df$QTL_SlopeScore[i],4,sep="_")}
    rm(QTL1,QTL2,QTL3,QTL4)
    if(QTL1_2 + QTL2_2 + QTL3_2 + QTL4_2 ==4){
      df$QTL_SlopeScore[i] <- paste(df$QTL_SlopeScore[i],4,sep="_")}
    rm(QTL1_2,QTL2_2,QTL3_2,QTL4_2)
  }
}

for (i in 1:length(df$eQTL_Gene_featureCounts)){
  x <- as.vector(df$QTL_SlopeScore[i])
  x = unlist(strsplit(x, split="_"))
  z <- match(4,x,nomatch = 0)
  if (z != 0){
    df$QTL_SlopeScore2[i] <- 4
  }
}


###3 QTLs same direction

for (i in 1:length(df$eQTL_Gene_featureCounts)){
  if(as.numeric(as.character(df$gene_nominal_pvalue_featureCounts[i])) + as.numeric(as.character(df$ATAC_peak_nominal_pvalue[i])) + as.numeric(as.character(df$K27ac_peak_nominal_pvalue[i])) + as.numeric(as.character(df$K4me3_peak_nominal_pvalue[i])) == 3){
    x <- as.vector(df$gene_slope_featureCounts[i])
    x = unlist(strsplit(x, split=","))
    a <- as.vector(df$ATAC_peak_slope[i])
    a = unlist(strsplit(a, split=","))  
    b <- as.vector(df$K27ac_peak_slope[i])
    b = unlist(strsplit(b, split=","))  
    c <- as.vector(df$K4me3_peak_slope[i])
    c = unlist(strsplit(c, split=","))
    QTL1 = c()
    QTL2 = c()
    QTL3 = c()
    QTL4 = c()
    for (xx in 1:length(x)){
      if (x[xx] > 0){QTL1 <- c(QTL1,x[xx])}}
    for (aa in 1: length(a)){
      if(a[aa] > 0){QTL2 <- c(QTL2,a[aa])}} 
    for (bb in 1: length(b)){
      if(b[bb] > 0){QTL3 <- c(QTL3,b[bb])}}
    for (cc in 1: length(c)){
      if(c[cc] > 0){QTL4 <- c(QTL4,c[cc])}}
    if (length(QTL1) > 0){
      QTL1 <- NULL
      QTL1 = 1
    }else{
      QTL1 <- NULL
      QTL1 = 0
    }
    if (length(QTL2) > 0){
      QTL2 <- NULL
      QTL2 = 1
    }else{
      QTL2 <- NULL
      QTL2 = 0      
    }
    if (length(QTL3) > 0){
      QTL3 <- NULL
      QTL3 = 1
    }else{
      QTL3 <- NULL
      QTL3 = 0      
    }
    if (length(QTL4) > 0){
      QTL4 <- NULL
      QTL4 = 1
    }else{
      QTL4 <- NULL
      QTL4 = 0
    }
    QTL1_2 = c()
    QTL2_2 = c()
    QTL3_2 = c()
    QTL4_2 = c()
    for (xx in 1:length(x)){
      if (x[xx] < 0){QTL1_2 <- c(QTL1_2,x[xx])}}
    for (aa in 1: length(a)){
      if(a[aa] < 0){QTL2_2 <- c(QTL2_2,a[aa])}} 
    for (bb in 1: length(b)){
      if(b[bb] < 0){QTL3_2 <- c(QTL3_2,b[bb])}}
    for (cc in 1: length(c)){
      if(c[cc] < 0){QTL4_2 <- c(QTL4_2,c[cc])}}
    if (length(QTL1_2) > 0){
      QTL1_2 <- NULL
      QTL1_2 = 1
    }else{
      QTL1_2 <- NULL
      QTL1_2 = 0
    }
    if (length(QTL2_2) > 0){
      QTL2_2 <- NULL
      QTL2_2 = 1
    }else{
      QTL2_2 <- NULL
      QTL2_2 = 0      
    }
    if (length(QTL3_2) > 0){
      QTL3_2 <- NULL
      QTL3_2 = 1
    }else{
      QTL3_2 <- NULL
      QTL3_2 = 0      
    }
    if (length(QTL4_2) > 0){
      QTL4_2 <- NULL
      QTL4_2 = 1
    }else{
      QTL4_2 <- NULL
      QTL4_2 = 0
    }
    if(QTL1 + QTL2 + QTL3 + QTL4 ==3){
      df$QTL_SlopeScore[i] <- paste(df$QTL_SlopeScore[i],3,sep="_")}
    rm(QTL1,QTL2,QTL3,QTL4)
    if(QTL1_2 + QTL2_2 + QTL3_2 + QTL4_2 ==3){
      df$QTL_SlopeScore[i] <- paste(df$QTL_SlopeScore[i],3,sep="_")}
    rm(QTL1_2,QTL2_2,QTL3_2,QTL4_2)
  }
}

for (i in 1:length(df$eQTL_Gene_featureCounts)){
  x <- as.vector(df$QTL_SlopeScore[i])
  x = unlist(strsplit(x, split="_"))
  z <- match(3,x,nomatch = 0)
  if (z != 0){
    df$QTL_SlopeScore2[i] <- 3
  }
}


###2 QTLs same direction

for (i in 1:length(df$eQTL_Gene_featureCounts)){
  if(as.numeric(as.character(df$gene_nominal_pvalue_featureCounts[i])) + as.numeric(as.character(df$ATAC_peak_nominal_pvalue[i])) + as.numeric(as.character(df$K27ac_peak_nominal_pvalue[i])) + as.numeric(as.character(df$K4me3_peak_nominal_pvalue[i])) == 2){
    x <- as.vector(df$gene_slope_featureCounts[i])
    x = unlist(strsplit(x, split=","))
    a <- as.vector(df$ATAC_peak_slope[i])
    a = unlist(strsplit(a, split=","))  
    b <- as.vector(df$K27ac_peak_slope[i])
    b = unlist(strsplit(b, split=","))  
    c <- as.vector(df$K4me3_peak_slope[i])
    c = unlist(strsplit(c, split=","))
    QTL1 = c()
    QTL2 = c()
    QTL3 = c()
    QTL4 = c()
    for (xx in 1:length(x)){
      if (x[xx] > 0){QTL1 <- c(QTL1,x[xx])}}
    for (aa in 1: length(a)){
      if(a[aa] > 0){QTL2 <- c(QTL2,a[aa])}} 
    for (bb in 1: length(b)){
      if(b[bb] > 0){QTL3 <- c(QTL3,b[bb])}}
    for (cc in 1: length(c)){
      if(c[cc] > 0){QTL4 <- c(QTL4,c[cc])}}
    if (length(QTL1) > 0){
      QTL1 <- NULL
      QTL1 = 1
    }else{
      QTL1 <- NULL
      QTL1 = 0
    }
    if (length(QTL2) > 0){
      QTL2 <- NULL
      QTL2 = 1
    }else{
      QTL2 <- NULL
      QTL2 = 0      
    }
    if (length(QTL3) > 0){
      QTL3 <- NULL
      QTL3 = 1
    }else{
      QTL3 <- NULL
      QTL3 = 0      
    }
    if (length(QTL4) > 0){
      QTL4 <- NULL
      QTL4 = 1
    }else{
      QTL4 <- NULL
      QTL4 = 0
    }
    QTL1_2 = c()
    QTL2_2 = c()
    QTL3_2 = c()
    QTL4_2 = c()
    for (xx in 1:length(x)){
      if (x[xx] < 0){QTL1_2 <- c(QTL1_2,x[xx])}}
    for (aa in 1: length(a)){
      if(a[aa] < 0){QTL2_2 <- c(QTL2_2,a[aa])}} 
    for (bb in 1: length(b)){
      if(b[bb] < 0){QTL3_2 <- c(QTL3_2,b[bb])}}
    for (cc in 1: length(c)){
      if(c[cc] < 0){QTL4_2 <- c(QTL4_2,c[cc])}}
    if (length(QTL1_2) > 0){
      QTL1_2 <- NULL
      QTL1_2 = 1
    }else{
      QTL1_2 <- NULL
      QTL1_2 = 0
    }
    if (length(QTL2_2) > 0){
      QTL2_2 <- NULL
      QTL2_2 = 1
    }else{
      QTL2_2 <- NULL
      QTL2_2 = 0      
    }
    if (length(QTL3_2) > 0){
      QTL3_2 <- NULL
      QTL3_2 = 1
    }else{
      QTL3_2 <- NULL
      QTL3_2 = 0      
    }
    if (length(QTL4_2) > 0){
      QTL4_2 <- NULL
      QTL4_2 = 1
    }else{
      QTL4_2 <- NULL
      QTL4_2 = 0
    }
    if(QTL1 + QTL2 + QTL3 + QTL4 ==2){
      df$QTL_SlopeScore[i] <- paste(df$QTL_SlopeScore[i],2,sep="_")}
    rm(QTL1,QTL2,QTL3,QTL4)
    if(QTL1_2 + QTL2_2 + QTL3_2 + QTL4_2 ==2){
      df$QTL_SlopeScore[i] <- paste(df$QTL_SlopeScore[i],2,sep="_")}
    rm(QTL1_2,QTL2_2,QTL3_2,QTL4_2)
  }
}

for (i in 1:length(df$eQTL_Gene_featureCounts)){
  x <- as.vector(df$QTL_SlopeScore[i])
  x = unlist(strsplit(x, split="_"))
  z <- match(2,x,nomatch = 0)
  if (z != 0){
    df$QTL_SlopeScore2[i] <- 2
  }
}

df$SlopesSameDirection <- ifelse(df$QTL_SlopeScore2 > 0, 1, 0)

#####Save modified table

write.table(df, file = "Summary_table_Test_LeadSlope_NoAllNA_WithSlopeScore_New.txt", sep="\t", row.names = FALSE, quote = FALSE)

######Loading once written

df <- read.table(file = "Summary_table_Test_LeadSlope_NoAllNA_WithSlopeScore_New.txt", sep = "\t", header = TRUE)

coloc_hits <- read.table("~/Desktop/Tregs//COLOC/treg_GWAS_coloc_hits.txt", sep="\t", header=TRUE)
coloc_hits <- coloc_hits[, (names(coloc_hits) %in% c("assay","trait","gwas_lead", "qtl_lead","phenotype_id","classification"))]
coloc_hits$COLOC <- paste(coloc_hits$assay, coloc_hits$trait, coloc_hits$gwas_lead, coloc_hits$phenotype_id, coloc_hits$classification, sep=":")
coloc_hits <- coloc_hits[,c("qtl_lead", "COLOC")]
coloc_hits <- aggregate(COLOC ~ qtl_lead, data = coloc_hits, paste, collapse = ",")
df <-  merge(df,coloc_hits,by.x = "eQTL_LEAD_featureCounts", by.y = "qtl_lead", all.x = TRUE)
df <- df[c(2:11,1,12:37)]
df[is.na(df)] <- 0
write.table(df, file = "Summary_table_Test_LeadSlope_NoAllNA_WithCOLOCInfo_New.txt", sep="\t", row.names = FALSE, quote = FALSE)

####Plotting using UpSetR library
categories <- c(1, 7, 8, 9, 14, 20, 26, 32, 36)
df2 <- df[,categories]
df3 <- df2[,-1]
rownames(df3) <- df2[,1]
df7 <- df3
df7$Category <- ifelse(df7$gene_nominal_pvalue_featureCounts == 0, "no_eQTL", "eQTL")
df7$Category <- ifelse(df7$K27ac_Peak_Overlap == 0, paste(df7$Category,"_noK27acPeak", sep=""), paste(df7$Category,"_K27acPeak", sep=""))
df7$Category <- ifelse(df7$K4me3_Peak_Overlap == 0, paste(df7$Category,"_noK4me3Peak", sep=""), paste(df7$Category,"_K4me3Peak", sep=""))
df7$Category <- ifelse(df7$ATAC_Peak_Overlap == 0, paste(df7$Category,"_noATACPeak", sep=""), paste(df7$Category,"_ATACPeak", sep=""))
df7$Category <- ifelse(df7$K27ac_peak_nominal_pvalue == 0, paste(df7$Category,"_noK27acQTL", sep=""), paste(df7$Category,"_K27acQTL", sep=""))
df7$Category <- ifelse(df7$K4me3_peak_nominal_pvalue == 0, paste(df7$Category,"_noK4me3QTL", sep=""), paste(df7$Category,"_K4me3QTL", sep=""))
df7$Category <- ifelse(df7$ATAC_peak_nominal_pvalue == 0, paste(df7$Category,"_noATACQTL", sep=""), paste(df7$Category,"_ATACQTL", sep=""))
df7$Category <- ifelse(df7$SlopesSameDirection == 0, paste(df7$Category,"_DifferentSlopeDirection", sep=""), paste(df7$Category,"_SameSlopeDirection", sep=""))
df8 <- as.data.frame(table(df7$Category))
write.table(df8, file = "Categories_New.txt", sep="\t", quote = FALSE, row.names = FALSE)

###Plotting SNP categories, no slopeScore considered
df9 <- df3
df9$Category <- ifelse(df9$gene_nominal_pvalue_featureCounts != 0, "eQTL", "noeQTL")
df9$Category <- ifelse(df9$K27ac_Peak_Overlap != 0 | df9$K4me3_Peak_Overlap != 0 | df9$ATAC_Peak_Overlap !=0 , paste(df9$Category,"_OverlapsatLeast1Peak", sep=""), paste(df9$Category,"_noPeakOverlap", sep=""))
df9$Category <- ifelse(df9$K27ac_peak_nominal_pvalue != 0 | df9$K4me3_peak_nominal_pvalue != 0 | df9$ATAC_peak_nominal_pvalue !=0 , paste(df9$Category,"_epiQTL", sep=""), paste(df9$Category,"_noepiQTL", sep=""))
df10 <- as.data.frame(table(df9$Category))
df10$Selection <- c("No_Focus","No_Focus","Focus","Focus", 
                    "No_Focus", "Focus", "No_Focus")
levels(df10$Var1)[c("eQTL_noPeakOverlap_epiQTL")] <- c("eQTL & regQTL")
levels(df10$Var1)[c("eQTL_noPeakOverlap_noepiQTL")] <- c("eQTL")
levels(df10$Var1)[c("eQTL_OverlapsatLeast1Peak_epiQTL")] <- c("eQTL & regQTL (in a peak)")
levels(df10$Var1)[c("eQTL_OverlapsatLeast1Peak_noepiQTL")] <- c("eQTL (in a peak)")
levels(df10$Var1)[c("noeQTL_noPeakOverlap_epiQTL")] <- c("regQTL")
levels(df10$Var1)[c("noeQTL_OverlapsatLeast1Peak_epiQTL")] <- c("regQTL (in a peak)")
levels(df10$Var1)[c("noeQTL_OverlapsatLeast1Peak_noepiQTL")] <- c("In a peak")
df10$Var1[df10$Var1 == "eQTL_noPeakOverlap_epiQTL"] <- c("eQTL & regQTL")
df10$Var1[df10$Var1 == "eQTL_noPeakOverlap_noepiQTL"] <- c("eQTL")
df10$Var1[df10$Var1 == "eQTL_OverlapsatLeast1Peak_epiQTL"] <- c("eQTL & regQTL (in a peak)")
df10$Var1[df10$Var1 == "eQTL_OverlapsatLeast1Peak_noepiQTL"] <- c("eQTL (in a peak)")
df10$Var1[df10$Var1 == "noeQTL_noPeakOverlap_epiQTL"] <- c("regQTL")
df10$Var1[df10$Var1 == "noeQTL_OverlapsatLeast1Peak_epiQTL"] <- c("regQTL (in a peak)")
df10$Var1[df10$Var1 == "noeQTL_OverlapsatLeast1Peak_noepiQTL"] <- c("In a peak")

df10$Var1 <- factor(df10$Var1, levels = df10$Var1[order(-df10$Freq)])
df11 <- df10[!(df10$Var1 == "noeQTL_noPeakOverlap_noepiQTL"),]
df11$Var1 <- factor(df11$Var1, levels = df11$Var1[order(-df11$Freq)])

library(ggplot2)
ggplot(df11, aes(x = reorder(Var1, -Freq), y = Freq, color=Selection, fill = Selection)) +
  geom_bar(stat = "identity")+
  geom_text(aes(y= Freq+Freq/10, label=Freq))+
#  facet_wrap( ~ Selection, ncol=2, scales = "free")+  
  theme_bw()+
  theme(axis.text.x = element_text(), axis.title.x=element_blank(),
        strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour="#EAE8E8"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        text = element_text(size=16),
        legend.position = "top")+
  labs(y="Number of SNPs")+
  scale_fill_manual(values = c("#848383","#848383"))+
  scale_color_manual(values = c("#6f28ff","#848383"))
ggsave("SNP_categories_New.pdf", plot = last_plot(), width = 12, height = 10, units = "in",
       dpi = 300,device = "pdf", useDingbats=FALSE)


df3$gene_nominal_pvalue_featureCounts <- as.numeric(as.character((df3$gene_nominal_pvalue_featureCounts)))
png("All_FunctionalCategories_New.png", height = 20, width = 25, units = c("cm"), res = 300)
upset(df3, nsets = 8, nintersects = 200, sets = rev(c("gene_nominal_pvalue_featureCounts",
                        "K27ac_peak_nominal_pvalue",
                        "K4me3_peak_nominal_pvalue",
                        "ATAC_peak_nominal_pvalue","K27ac_Peak_Overlap",
                        "K4me3_Peak_Overlap","ATAC_Peak_Overlap",
                        "SlopesSameDirection")),
      #queries = list(list(query = intersects, params = list("gene_nominal_pvalue",
      #                                                      "ATAC_peak_nominal_pvalue",
      #                                                      "H3K27ac_peak_SNP_nominal_pvalue",
      #                                                      "H3K4me3_peak_SNP_nominal_pvalue",
      #                                                      "ATAC_Peak_Overlap",
      #                                                      "K27ac_Peak_Overlap",
      #                                                      "K4me3_Peak_Overlap",
      #                                                      "SlopesSameDirection"), 
      #                    color = "red", active = T),
      #               list(query = intersects, params = list("ATAC_peak_nominal_pvalue",
      #                                                      "H3K27ac_peak_SNP_nominal_pvalue",
      #                                                      "H3K4me3_peak_SNP_nominal_pvalue",
      #                                                      "ATAC_Peak_Overlap",
      #                                                      "K27ac_Peak_Overlap",
      #                                                      "K4me3_Peak_Overlap",
      #                                                      "SlopesSameDirection"), 
      #                    color = "orange", active = T)
      #),
      number.angles = 30, point.size = 1, line.size = 1, 
      mainbar.y.label = "Category Intersections", sets.x.label = "SNPs Per Category", 
      text.scale = c(1, 1, 1, 1, 1, 0.75), order.by = "freq",
      keep.order = TRUE)
dev.off()

df4 <- df3[ which( df3$gene_nominal_pvalue_featureCounts == 1) , ]
df5 <- df3[ which( df3$gene_nominal_pvalue_featureCounts == 0) , ]
df6 <- df3[ which( df3$K4me3_Peak_Overlap == 1) , ]

png("In_K4me3Peak_New.png", height = 20, width = 25, units = c("cm"), res = 300)
upset(df6, nsets = 8, nintersects = 200, sets = rev(c("gene_nominal_pvalue_featureCounts",
                                                      "K27ac_peak_nominal_pvalue","K4me3_peak_nominal_pvalue",
                                                      "ATAC_peak_nominal_pvalue","K27ac_Peak_Overlap",
                                                      "K4me3_Peak_Overlap","ATAC_Peak_Overlap",
                                                      "SlopesSameDirection")),
      queries = list(list(query = intersects, params = list("gene_nominal_pvalue_featureCounts",
                                                            "ATAC_peak_nominal_pvalue",
                                                            "K27ac_peak_nominal_pvalue",
                                                            "K4me3_peak_nominal_pvalue",
                                                            "ATAC_Peak_Overlap",
                                                            "K27ac_Peak_Overlap",
                                                            "K4me3_Peak_Overlap",
                                                            "SlopesSameDirection"), 
                          color = "red", active = T),
                     list(query = intersects, params = list("gene_nominal_pvalue_featureCounts",
                                                            "ATAC_peak_nominal_pvalue",
                                                            "K27ac_peak_nominal_pvalue",
                                                            "K4me3_peak_nominal_pvalue",
                                                            "ATAC_Peak_Overlap",
                                                            "K27ac_Peak_Overlap",
                                                            "K4me3_Peak_Overlap",
                                                            "SlopesSameDirection"), 
                          color = "green", active = T)),
      number.angles = 30, point.size = 1, line.size = 0.5, 
      mainbar.y.label = "Category Intersections", sets.x.label = "SNPs Per Category", 
      text.scale = c(1, 1, 1, 1, 1, 0.75), order.by = "freq",
      keep.order = TRUE)
dev.off()

png("Lead_LDBLOCK_FunctionalCategories_New.png", height = 20, width = 25, units = c("cm"), res = 300)
upset(df4, nsets = 8, nintersects = 200, sets = rev(c("gene_nominal_pvalue_featureCounts",
                    "K27ac_peak_nominal_pvalue","K4me3_peak_nominal_pvalue",
                    "ATAC_peak_nominal_pvalue","K27ac_Peak_Overlap",
                    "K4me3_Peak_Overlap","ATAC_Peak_Overlap",
                    "SlopesSameDirection")),
      queries = list(list(query = intersects, params = list("gene_nominal_pvalue_featureCounts",
                                                            "ATAC_peak_nominal_pvalue",
                                                            "K27ac_peak_nominal_pvalue",
                                                            "K4me3_peak_nominal_pvalue",
                                                            "ATAC_Peak_Overlap",
                                                            "K27ac_Peak_Overlap",
                                                            "K4me3_Peak_Overlap",
                                                            "SlopesSameDirection"), 
                          color = "red", active = T)),
      number.angles = 30, point.size = 1, line.size = 0.5, 
      mainbar.y.label = "Category Intersections", sets.x.label = "SNPs Per Category", 
      text.scale = c(1, 1, 1, 1, 1, 0.75), order.by = "freq",
      keep.order = TRUE)
dev.off()

png("SNP_NO_Lead_LDBLOCK_FunctionalCategories_New.png", height = 20, width = 25, units = c("cm"), res = 300)
upset(df5, nsets = 7, nintersects = 200, sets = rev(c("H3K27ac_peak_SNP_nominal_pvalue",
                        "H3K4me3_peak_SNP_nominal_pvalue",
                        "ATAC_peak_nominal_pvalue","K27ac_Peak_Overlap",
                        "K4me3_Peak_Overlap","ATAC_Peak_Overlap",
                        "SlopesSameDirection")),
      queries = list(list(query = intersects, params = list("ATAC_peak_nominal_pvalue",
                                                            "H3K27ac_peak_SNP_nominal_pvalue",
                                                            "H3K4me3_peak_SNP_nominal_pvalue",
                                                            "ATAC_Peak_Overlap",
                                                            "K27ac_Peak_Overlap",
                                                            "K4me3_Peak_Overlap",
                                                            "SlopesSameDirection"), 
                          color = "orange", active = T)),
      number.angles = 30, point.size = 1, line.size = 1, 
      mainbar.y.label = "Category Intersections", sets.x.label = "SNPs Per Category", 
      text.scale = c(1, 1, 1, 1, 1, 0.75), order.by = "freq",
      keep.order = TRUE)
dev.off()

###Select SNPs by rsID that have and don't eQTls
eQTL_SNPs <- df[ which( df$eQTL_LEAD_BLOCK_featureCounts == 1), ]
length(eQTL_SNPs$SNPposition)
noeQTL_SNPs <- df[ which( df$eQTL_LEAD_BLOCK_featureCounts == 0), ]
length(noeQTL_SNPs$SNPposition)

###Get the LD scores for those SNPs: x and y
eQTL_SNPs_Scores <- LDScores[ which(LDScores$SNP %in% eQTL_SNPs$ID),]
x <- eQTL_SNPs_Scores[c("SNP","LD.SCORE")]
length(x$SNP)
noeQTL_SNPs_Scores <- LDScores[ which(LDScores$SNP %in% noeQTL_SNPs$ID),]
y <- noeQTL_SNPs_Scores[c("SNP","LD.SCORE")]
length(y$SNP)
df6 <- matchCases(xcase = x$LD.SCORE, ycase = rep(0,length(x$LD.SCORE)), idcase = x$SNP,
                  xcontrol = y$LD.SCORE, ycontrol = rep(0,length(y$LD.SCORE)), idcontrol = y$SNP, tol=30,
                  maxmatch=1, which=c('closest','random'))
df7 <- df6[ which( df6$type == "control"),]
noeQTL_SNPs_Selection <- noeQTL_SNPs[ which(noeQTL_SNPs$ID %in% df7$id),]
length(noeQTL_SNPs_Selection$SNPposition)
df8 <- df6[ which( df6$type == "case"),]
eQTL_SNPs_Selection <- eQTL_SNPs[ which(eQTL_SNPs$ID %in% df8$id),]
length(eQTL_SNPs_Selection$SNPposition)

png("Lead_LDBLOCK_FunctionalCategories_Selection_New.png", height = 20, width = 25, units = c("cm"), res = 300)
upset(eQTL_SNPs_Selection, nsets = 8, nintersects = 200, sets = rev(c("gene_nominal_pvalue",
                                                      "H3K27ac_peak_SNP_nominal_pvalue","H3K4me3_peak_SNP_nominal_pvalue",
                                                      "ATAC_peak_nominal_pvalue","K27ac_Peak_Overlap",
                                                      "K4me3_Peak_Overlap","ATAC_Peak_Overlap",
                                                      "SlopesSameDirection")),
      queries = list(list(query = intersects, params = list("gene_nominal_pvalue",
                                                            "ATAC_peak_nominal_pvalue",
                                                            "H3K27ac_peak_SNP_nominal_pvalue",
                                                            "H3K4me3_peak_SNP_nominal_pvalue",
                                                            "ATAC_Peak_Overlap",
                                                            "K27ac_Peak_Overlap",
                                                            "K4me3_Peak_Overlap",
                                                            "SlopesSameDirection"), 
                          color = "red", active = T)),
      number.angles = 30, point.size = 0.5, line.size = 0.5, 
      mainbar.y.label = "Category Intersections", sets.x.label = "SNPs Per Category", 
      text.scale = c(1, 1, 1, 1, 1, 0.75), order.by = "freq",
      keep.order = TRUE)
dev.off()

png("SNP_NO_Lead_LDBLOCK_FunctionalCategories_Selection_New.png", height = 20, width = 25, units = c("cm"), res = 300)
upset(noeQTL_SNPs_Selection, nsets = 7, nintersects = 200, sets = rev(c("H3K27ac_peak_SNP_nominal_pvalue",
                                                      "H3K4me3_peak_SNP_nominal_pvalue",
                                                      "ATAC_peak_nominal_pvalue","K27ac_Peak_Overlap",
                                                      "K4me3_Peak_Overlap","ATAC_Peak_Overlap",
                                                      "SlopesSameDirection")),
      #queries = list(list(query = intersects, params = list("ATAC_peak_nominal_pvalue",
      #                                                      "H3K27ac_peak_SNP_nominal_pvalue",
      #                                                      "H3K4me3_peak_SNP_nominal_pvalue",
      #                                                      "ATAC_Peak_Overlap",
      #                                                      "K27ac_Peak_Overlap",
      #                                                      "K4me3_Peak_Overlap",
      #                                                      "SlopesSameDirection"), 
      #                    color = "orange", active = T)),
      number.angles = 30, point.size = 1, line.size = 1, 
      mainbar.y.label = "Category Intersections", sets.x.label = "SNPs Per Category", 
      text.scale = c(1, 1, 1, 1, 1, 0.75), order.by = "freq",
      keep.order = TRUE)
dev.off()

