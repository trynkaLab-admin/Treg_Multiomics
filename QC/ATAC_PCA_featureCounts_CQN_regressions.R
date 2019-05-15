library("ggplot2")
library("gplots")
library("cqn")
library("magrittr")

input=c("ATAC_Counts_featureCounts")
file=paste(input,".txt",sep="")
data <- read.table(file=file, row.names=1, header=TRUE, sep="\t")
counts_matrix <- data[rowSums(data) > 0, ]
gene_metadata <- read.table("Peak_MetaData.txt", header=TRUE, sep="\t")

ggplot(gene_metadata, aes(x=length))+
  geom_histogram(binwidth=100,color="#bb814f", alpha = 0.9, fill="#bb814f")+ 
  geom_vline(xintercept = median(gene_metadata$length),color="#bb814f")+
  theme(strip.background=element_blank(),panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text=element_text(size=12))+
  annotate("text", x = median(gene_metadata$length)+2000, y = 3000, label = paste("Median peak length = ", round(median(gene_metadata$length),digits=2),sep=""))

ggsave(paste(input,"_PeakLength_Paper.pdf",sep=""), height=5, width=5, dpi=300, units = c("cm"))



###Number of reads per sample and peak
breaks=c(seq(0,5000,100),max(rowMeans(data)))
ggplot(data,aes(x=rowMeans(data)))+
  geom_histogram(breaks=breaks, color="#bb814f", fill="#bb814f", alpha = 0.2)+
  coord_cartesian(xlim = c(0, 6000)) +
  geom_vline(xintercept = median(rowMeans(data)), color="#bb814f")+
  theme(text = element_text(size=12),strip.background=element_blank(),panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  annotate("text", x = median(rowMeans(data))+1000, y = 7500, label = paste("Median reads per peak = ", round(median(rowMeans(data)),digits=2),sep=""))
ggsave(paste("MedianReads_",input,".pdf",sep=""), height=5, width=5, dpi=300, units = c("cm"))


###CQN method
calculateCQN <- function(counts_matrix, gene_metadata, return_type = "normalised"){
  #Normalize read counts using the CQN method.
  #sizeFactors: An optional vector of sizeFactors, 
  #ie. the sequencing effort of the various samples. 
  #If NULL this is calculated as the column sums of counts.
  expression_cqn = cqn::cqn(counts = counts_matrix[gene_metadata$gene_id,], 
                            x = gene_metadata$percentage_gc_content, 
                            lengths = gene_metadata$length, verbose = TRUE)
  #Choose return type
  if(return_type == "normalised"){
    expression_norm = expression_cqn$y + expression_cqn$offset
    return(expression_norm)
  }
  else{
    return(expression_cqn)
  }
}

counts_matrix_cqn <- calculateCQN(counts_matrix = counts_matrix,gene_metadata = gene_metadata)

###Save normalized counts
write.table(counts_matrix_cqn, file=paste(input,"_NoBatchCorrection_CQN.txt",
                                          sep=""), sep = "\t", row.names=T)
#####Covariates
covariates <- read.table("Covariates.txt",header=TRUE, row.names=1)
covariates_filtered <- covariates[which(row.names(covariates) %in% colnames(counts_matrix_cqn)),]
counts_matrix_cqn2 <- read.table(paste(input,"_NoBatchCorrection_CQN.txt",sep=""), sep="\t", header = TRUE)

##############################################
#PCA without batch correction

data.pca.nocor <- prcomp(t(counts_matrix_cqn2), center = TRUE, scale = TRUE)
write.table(data.frame(data.pca.nocor$x[,1:30]), file=paste(input,"_noBatchCorrection_CQN.PCA",sep=""), sep = "\t", row.names=T)
influence.nocor <- abs(data.pca.nocor$rotation)
SNPinfluence.nocor <- as.matrix(sweep(influence.nocor, 2, colSums(influence.nocor), "/"))
write.table(SNPinfluence.nocor, file=paste(input,"_NoBatchCorrection.BinContribution",sep=""), sep = "\t", row.names=T)
png(file=paste("PCs_Variance_NoBatchCorrection_",input,".png",sep=""), height=20, width=40, unit="cm", res=300)
plot(data.pca.nocor, type = "l")
dev.off()
summary(data.pca.nocor)

eigs <- data.pca.nocor$sdev^2

varpercent <- vector(mode="numeric", length=30)
for (i in 1:30){
  varpercent[i]<- (eigs[i] / sum(eigs))*100
}

varpercent2 <- data.frame(matrix(ncol = 2, nrow = 30))
x <- c("PCs", "PerCent")
colnames(varpercent2) <- x
varpercent2$PCs <- seq(1,30,1)
varpercent2$PerCent <- varpercent


ggplot(varpercent2, aes(x = PCs, y = PerCent,group=1)) +
  geom_point(shape=21, fill="blue", size=2)+geom_line()+ylab("% Variance Explained")+
  geom_hline(yintercept = 1)+
  theme_bw()

ggsave(paste(input,"_PercentVarianceExplained_PCs_onlyCQN.png", sep = ""), height = 10, width = 10, units = c("cm"), dpi = 300)

scores.nocor <- data.frame(data.pca.nocor$x[,1:30])
Forplotting.nocor <- merge(scores.nocor, covariates, by="row.names")
percentVar.nocor <- round(100*data.pca.nocor$sdev^2/sum(data.pca.nocor$sdev^2),1)

ggplot(Forplotting.nocor,aes(scores, x=PC1, y=PC2, colour=as.factor(Forplotting.nocor$Sex)))+ 
  geom_point() + labs(color='Gender') +
  geom_text(aes(label= Forplotting.nocor$Row.names),hjust=0, vjust=0, size = 4)+
  theme(text = element_text(size=24))
ggsave(paste("PC1_PC2_",input,"_ColoredBySex_NoBatchCor.png"), height=10, width=10, dpi=300)

ggplot(Forplotting.nocor,aes(scores, x=PC1, y=PC2, colour=as.factor(Forplotting.nocor$Month_batch)))+ 
  geom_point() + labs(color='Batch') +geom_text(aes(label= Forplotting.nocor$Row.names),hjust=0, vjust=0, size = 4)
ggsave(paste("PC1_PC2_",input,"_ColoredByBatch_NoBatchCor.png"), height=10, width=10, dpi=300)

qplot(data=Forplotting.nocor, PC1, PC2, color= Month_batch,main="PC1 vs PC2", size =I(4))+
  labs(x=paste0("PC1,VarExp:",round(percentVar.nocor[1],4)), y=paste0("PC2,VarExp:",round(percentVar.nocor[2],4)))+ 
  theme(strip.background=element_blank(),panel.grid.major.x = element_line(colour="gray50"), panel.grid.major.y = element_line(colour="gray50"), panel.grid.minor = element_blank(),panel.background = element_blank(),text=element_text(size=20))+scale_color_manual(values=c("blue","purple","darkorange","red","black","green","orange"),guide=guide_legend(title = "Batch"))
ggsave(paste(input,"_PerCentVariance_Batch_NoBatchCor.png",sep=""), height=10, width=10, dpi=300)

qplot(data=Forplotting.nocor, PC1, PC2, color= Sex,main="PC1 vs PC2", size =I(4))+ 
  labs(x=paste0("PC1,VarExp:",round(percentVar.nocor[1],4)), y=paste0("PC2,VarExp:",round(percentVar.nocor[2],4)))+ 
  theme(strip.background=element_blank(),panel.grid.major.x = element_line(colour="gray50"), panel.grid.major.y = element_line(colour="gray50"), panel.grid.minor = element_blank(),panel.background = element_blank(),text=element_text(size=20))+scale_color_manual(values=c("blue","purple","darkorange","red","black","green","orange"),guide=guide_legend(title = "Gender"))
ggsave(paste(input,"_PerCentVariance_Sex_NoBatchCor.png",sep=""), height=10, width=10, dpi=300)

aov1 = aov(Forplotting.nocor$PC1 ~ Forplotting.nocor$Month_batch)
summary(aov1)
aov2 = aov(Forplotting.nocor$PC1 ~ Forplotting.nocor$Sex)
summary(aov2)

aov3 = aov(Forplotting.nocor$PC2 ~ Forplotting.nocor$Month_batch)
summary(aov3)
aov4 = aov(Forplotting.nocor$PC2 ~ Forplotting.nocor$Sex)
summary(aov4)

###Check effect of regression out PCs

regressPrinicpalComponents <- function(data_matrix, n_pcs){
  #Regress out first n principal components from the data matrix
  if(n_pcs > 0){
    pca = prcomp(t(data_matrix))
    pca_explained = (pca$x[,1:n_pcs] %*% t(pca$rotation[,1:n_pcs])) %>%
      scale(center = -1 * pca$center, scale = FALSE) %>% t()
    result = data_matrix - pca_explained
  } else{
    result = data_matrix
  }
  return(result)
}

regressedPCs <- regressPrinicpalComponents(data_matrix = counts_matrix_cqn2, n_pcs = 30)
###Save CQN counts with PCs regressed

write.table(regressedPCs, file=paste(input,"_CQN_30PCsRegressed.txt",
                                     sep=""), sep = "\t", row.names=T)


data.pcaregressed <- prcomp(t(regressedPCs), center = TRUE, scale = TRUE)
influence <- abs(data.pcaregressed$rotation)

percentVar.regressed <- round(100*data.pcaregressed$sdev^2/sum(data.pcaregressed$sdev^2),1)
#PLOT

scores_2 <- data.frame(data.pcaregressed$x[,1:30])
Forplotting2 <- merge(scores_2, covariates, by="row.names")


qplot(data=Forplotting2, PC1, PC2, color= Month_batch,main="PC1 vs PC2", size =I(4))+
  labs(x=paste0("PC1,VarExp:",round(percentVar.nocor[1],4)), y=paste0("PC2,VarExp:",round(percentVar.nocor[2],4)))+ 
  theme(strip.background=element_blank(),panel.grid.major.x = element_line(colour="gray50"), panel.grid.major.y = element_line(colour="gray50"), panel.grid.minor = element_blank(),panel.background = element_blank(),text=element_text(size=20))+scale_color_manual(values=c("blue","purple","darkorange","red","black","green","orange"),guide=guide_legend(title = "Batch"))

qplot(data=Forplotting2, PC1, PC2, color= Sex,main="PC1 vs PC2", size =I(4))+ 
  labs(x=paste0("PC1,VarExp:",round(percentVar.nocor[1],4)), y=paste0("PC2,VarExp:",round(percentVar.nocor[2],4)))+ 
  theme(strip.background=element_blank(),panel.grid.major.x = element_line(colour="gray50"), panel.grid.major.y = element_line(colour="gray50"), panel.grid.minor = element_blank(),panel.background = element_blank(),text=element_text(size=20))+scale_color_manual(values=c("blue","purple","darkorange","red","black","green","orange"),guide=guide_legend(title = "Gender"))

aov13 = aov(Forplotting2$PC1 ~ Forplotting2$Month_batch)
summary(aov13)
aov14 = aov(Forplotting2$PC2 ~ Forplotting2$Sex)
summary(aov14)
aov15 = aov(Forplotting2$PC1 ~ Forplotting2$Sex)
summary(aov15)
aov16 = aov(Forplotting2$PC2 ~ Forplotting2$Month_batch)
summary(aov16)
aov17 = aov(Forplotting2$PC3 ~ Forplotting2$Month_batch)
summary(aov17)
aov18 = aov(Forplotting2$PC3 ~ Forplotting2$Sex)
summary(aov18)
aov19 = aov(Forplotting2$PC4 ~ Forplotting2$Month_batch)
summary(aov19)
aov20 = aov(Forplotting2$PC4 ~ Forplotting2$Sex)
summary(aov20)
aov21 = aov(Forplotting2$PC5 ~ Forplotting2$Month_batch)
summary(aov21)
aov22 = aov(Forplotting2$PC5 ~ Forplotting2$Sex)
summary(aov22)
aov23 = aov(Forplotting2$PC6 ~ Forplotting2$Month_batch)
summary(aov23)
aov24 = aov(Forplotting2$PC6 ~ Forplotting2$Sex)
summary(aov24)
aov25 = aov(Forplotting2$PC7 ~ Forplotting2$Month_batch)
summary(aov25)
aov26 = aov(Forplotting2$PC7 ~ Forplotting2$Sex)
summary(aov26)
aov27 = aov(Forplotting2$PC8 ~ Forplotting2$Month_batch)
summary(aov27)
aov28 = aov(Forplotting2$PC8 ~ Forplotting2$Sex)
summary(aov28)
aov29 = aov(Forplotting2$PC9 ~ Forplotting2$Month_batch)
summary(aov29)
aov30 = aov(Forplotting2$PC9 ~ Forplotting2$Sex)
summary(aov30)
aov31 = aov(Forplotting2$PC10 ~ Forplotting2$Month_batch)
summary(aov31)
aov32 = aov(Forplotting2$PC10 ~ Forplotting2$Sex)
summary(aov32)

#Paper figures
ggplot(Forplotting.nocor, aes(x=PC1, y=PC2)) + 
  geom_point(color="#bb814f", alpha = 0.5, size = 3)+ 
  labs(x=paste0("PC1,VarExp:",round(percentVar.nocor[1],4)), y=paste0("PC2,VarExp:",round(percentVar.nocor[2],4)))+ 
  theme(strip.background=element_blank(),panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text=element_text(size=12))
ggsave(paste(input,"_PerCentVariance_NoBatchCor_Paper.pdf",sep=""), height=5, width=5, dpi=300, units = c("cm"))

ggplot(Forplotting2, aes(x=PC1, y=PC2)) + 
  geom_point(color="#bb814f", alpha = 0.9, size = 3)+ 
  labs(x=paste0("PC1,VarExp:",round(percentVar.regressed[1],4)), y=paste0("PC2,VarExp:",round(percentVar.regressed[2],4)))+ 
  theme(strip.background=element_blank(),panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text=element_text(size=12))
ggsave(paste(input,"_PerCentVariance_RegressedPCs_Paper.pdf",sep=""), height=5, width=5, dpi=300, units = c("cm"))
