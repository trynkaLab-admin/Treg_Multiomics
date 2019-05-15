library(ggplot2)
library(reshape)

D1=read.table("density.ATAC.around.All.txt", head=FALSE, stringsAsFactors=FALSE)
names(D1) <- c("Bin1","Bin2","ATAC")
D2=read.table("density.K27ac.around.All.txt", head=FALSE, stringsAsFactors=FALSE)
names(D2) <- c("Bin1","Bin2","K27ac")
D3=read.table("density.K4me3.around.All.txt", head=FALSE, stringsAsFactors=FALSE)
names(D3) <- c("Bin1","Bin2","K4me3")
D4 = merge(D1, D2, by=c("Bin1","Bin2"))
D4 = merge(D4, D3, by=c("Bin1","Bin2"))
D5 <- melt(D4, id=c("Bin1", "Bin2"))

annotation <- ggplot(data=D5, aes(x=(D5$Bin1+D5$Bin2)/2, y=D5$value, colour=D5$variable)) +
  geom_line()+
  xlab("Distance to QTLs")+ ylab("#annotations/kb")+
  scale_colour_manual(values=c("#bb814f","#7bbcb5","#ee7465"), name = "Annotation")+
  theme_bw()+
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
ggsave("Annotation_Distribution_All.pdf", plot = annotation, width = 15, height = 10, units = "cm",
       dpi = 300,device = "pdf", useDingbats=FALSE)