library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("ChIPpeakAnno")
library("VennDiagram")
library("ggplot2")
library("RColorBrewer")
library("scales")
data("TSS.human.GRCh38")

peak1 <- read.table("Merged_Day0_treg_ChM_K27ac_sorted-broad_peaks.broadPeak", sep="\t", header=TRUE) 

gr <- GRanges(seqnames = peak1$chr, strand = peak1$strand,
              ranges = IRanges(start = peak1$start, width = peak1$end-peak1$start))
values(gr) <- DataFrame(score = peak1$score, name = peak1$name)
gr

annotatedPeak <- annotatePeakInBatch(gr,
                                     AnnotationData = TSS.human.GRCh38)

peak2 <- read.table("Merged_Day0_treg_ChM_K27ac_sorted-broad_peaks_qval_0.001_FE_2.broadPeak", sep="\t", header=TRUE) 
gr2 <- GRanges(seqnames = peak2$chr, strand = peak2$strand,
               ranges = IRanges(start = peak2$start, width = peak2$end-peak2$start))
values(gr2) <- DataFrame(score = peak2$score, name = peak2$name)
gr2
annotatedPeak2 <- annotatePeakInBatch(gr2, AnnotationData = TSS.human.GRCh38)

peak3 <- read.table("Merged_Day0_treg_ChM_K4me3_sorted-sharp_peaks.narrowPeak", sep="\t", header=TRUE) 

gr3 <- GRanges(seqnames = peak3$chr, strand = peak3$strand,
              ranges = IRanges(start = peak3$start, width = peak3$end-peak3$start))
values(gr3) <- DataFrame(score = peak3$score, name = peak3$name)
gr3

annotatedPeak3 <- annotatePeakInBatch(gr3,
                                     AnnotationData = TSS.human.GRCh38)

peak4 <- read.table("Merged_Day0_treg_ChM_K4me3_sorted-sharp_peaks_qval_0.001_FE_2.narrowPeak", sep="\t", header=TRUE) 
gr4 <- GRanges(seqnames = peak4$chr, strand = peak4$strand,
               ranges = IRanges(start = peak4$start, width = peak4$end-peak4$start))
values(gr4) <- DataFrame(score = peak4$score, name = peak4$name)
gr4
annotatedPeak4 <- annotatePeakInBatch(gr4, AnnotationData = TSS.human.GRCh38)

#Peak length
peak1$length <-peak1$end - peak1$start
peak2$length <-peak2$end - peak2$start
peak3$length <-peak3$end - peak3$start
peak4$length <-peak4$end - peak4$start
plot(hist(as.numeric(peak1$length), breaks=seq(0,max(peak1$length)+500,50)), xlim=c(0,5000))
abline(v=median(peak1$length))
text(median(peak1$length)+500,3000,label=paste("Median H3K27ac All =",median(peak1$length),sep=""))
plot(hist(as.numeric(peak2$length), breaks=seq(0,max(peak2$length)+500,50)), xlim=c(0,10000))
abline(v=median(peak2$length))
text(median(peak2$length)+500,200,label=paste("Median H3K27ac Selected =",median(peak2$length),sep=""))
plot(hist(as.numeric(peak3$length), breaks=seq(0,max(peak3$length)+500,50)), xlim=c(0,5000))
abline(v=median(peak3$length))
text(median(peak3$length)+500,3000,label=paste("Median H3K4me3 All =",median(peak3$length),sep=""))
plot(hist(as.numeric(peak4$length), breaks=seq(0,max(peak4$length)+500,50)), xlim=c(0,5000))
abline(v=median(peak4$length))
text(median(peak4$length)+500,3000,label=paste("Median H3K4me3 Selected =",median(peak4$length),sep=""))


#Annotation
write.table(annotatedPeak, file="AnnotatedPeak_H3K27ac_all.txt", sep="\t")
write.table(annotatedPeak2, file="AnnotatedPeak_H3K27ac_selected.txt", sep="\t")
write.table(annotatedPeak3, file="AnnotatedPeak_H3K4me3_all.txt", sep="\t")
write.table(annotatedPeak4, file="AnnotatedPeak_H3K4me3_selected.txt", sep="\t")

pie1(table(as.data.frame(annotatedPeak)$insideFeature))
pie1(table(as.data.frame(annotatedPeak2)$insideFeature))
pie1(table(as.data.frame(annotatedPeak3)$insideFeature))
pie1(table(as.data.frame(annotatedPeak4)$insideFeature))

dist1 <- hist(annotatedPeak$distancetoFeature, breaks=c(as.numeric(min(annotatedPeak$distancetoFeature)), -50000,-10000,-3000,-500,500,3000,10000,50000,as.numeric(max(annotatedPeak$distancetoFeature))))
bin <- as.data.frame(dist1$counts)
bin$labels<-c("<-50000", "-50000:-10000", "-10000:-3000","-3000:-500","-500:500","500:3000","3000:10000","10000:50000",">50000")
bin$labels <- factor(bin$labels, levels = c("<-50000", "-50000:-10000", "-10000:-3000","-3000:-500","-500:500","500:3000","3000:10000","10000:50000",">50000"))
names(bin) <- c("counts", "labels")

ggplot(data=bin,aes(x=bin$labels, y= (as.numeric(bin$counts)/sum(as.numeric(bin$counts)))*100)) +
  geom_bar(stat="identity", aes(fill=bin$labels))+
  scale_x_discrete(limits = bin$labels) +
  coord_flip()+
  xlab("Peak location")+
  ylab("%")+
  labs(fill="Location")+
  theme_bw()

dist2 <- hist(annotatedPeak2$distancetoFeature, breaks=c(as.numeric(min(annotatedPeak2$distancetoFeature)), -50000,-10000,-3000,-500,500,3000,10000,50000,as.numeric(max(annotatedPeak2$distancetoFeature))))
bin2 <- as.data.frame(dist2$counts)
bin2$labels<-c("<-50000", "-50000:-10000", "-10000:-3000","-3000:-500","-500:500","500:3000","3000:10000","10000:50000",">50000")
bin2$labels <- factor(bin2$labels, levels = c("<-50000", "-50000:-10000", "-10000:-3000","-3000:-500","-500:500","500:3000","3000:10000","10000:50000",">50000"))
names(bin2) <- c("counts", "labels")

ggplot(data=bin2,aes(x=bin2$labels, y= (as.numeric(bin2$counts)/sum(as.numeric(bin2$counts)))*100)) +
  geom_bar(stat="identity", aes(fill=bin2$labels))+
  scale_x_discrete(limits = bin2$labels) +
  coord_flip()+
  xlab("Peak location")+
  ylab("%")+
  labs(fill="Location")+
  theme_bw()

dist3 <- hist(annotatedPeak3$distancetoFeature, breaks=c(as.numeric(min(annotatedPeak3$distancetoFeature)),
                                                         -50000,-10000,-3000,-500,500,3000,10000,50000,
                                                         as.numeric(max(annotatedPeak3$distancetoFeature))))
bin3 <- as.data.frame(dist3$counts)
bin3$labels<-c("<-50000", "-50000:-10000", "-10000:-3000","-3000:-500","-500:500","500:3000","3000:10000","10000:50000",">50000")
bin3$labels <- factor(bin3$labels, levels = c("<-50000", "-50000:-10000", "-10000:-3000","-3000:-500","-500:500","500:3000","3000:10000","10000:50000",">50000"))
names(bin3) <- c("counts", "labels")

ggplot(data=bin3,aes(x=bin3$labels, y= (as.numeric(bin3$counts)/sum(as.numeric(bin3$counts)))*100)) +
  geom_bar(stat="identity", aes(fill=bin3$labels))+
  scale_x_discrete(limits = bin3$labels) +
  coord_flip()+
  xlab("Peak location")+
  ylab("%")+
  labs(fill="Location")+
  theme_bw()

dist4 <- hist(annotatedPeak4$distancetoFeature, breaks=c(as.numeric(min(annotatedPeak4$distancetoFeature)),
                                                         -50000,-10000,-3000,-500,500,3000,10000,50000,
                                                         as.numeric(max(annotatedPeak4$distancetoFeature))))
bin4 <- as.data.frame(dist4$counts)
bin4$labels<-c("<-50000", "-50000:-10000", "-10000:-3000","-3000:-500","-500:500","500:3000","3000:10000","10000:50000",">50000")
bin4$labels <- factor(bin4$labels, levels = c("<-50000", "-50000:-10000", "-10000:-3000","-3000:-500","-500:500","500:3000","3000:10000","10000:50000",">50000"))
names(bin4) <- c("counts", "labels")

ggplot(data=bin4,aes(x=bin4$labels, y= (as.numeric(bin4$counts)/sum(as.numeric(bin4$counts)))*100)) +
  geom_bar(stat="identity", aes(fill=bin4$labels))+
  scale_x_discrete(limits = bin4$labels) +
  coord_flip()+
  xlab("Peak location")+
  ylab("%")+
  labs(fill="Location")+
  theme_bw()



genes1 = as.data.frame(annotatedPeak)$feature
genes2 = as.data.frame(annotatedPeak2)$feature
genes3 = as.data.frame(annotatedPeak3)$feature
genes4 = as.data.frame(annotatedPeak4)$feature


overlap1 <- calculate.overlap(x=list("K27ac_all"=unique(genes1),
                                    "K27ac_Selected"=unique(genes2)))
v1 <- draw.pairwise.venn(length(overlap1$a1),length(overlap1$a2),length(overlap1$a3),
                        c("K27ac_all","K27ac_Selected"),
                        fill = c("lightblue","darkblue"))

overlap2 <- calculate.overlap(x=list("K4me3_all"=unique(genes3),
                                     "K4me3_Selected"=unique(genes4)))
v2 <- draw.pairwise.venn(length(overlap2$a1),length(overlap2$a2),length(overlap2$a3),
                         c("K4me3_all","K4me3_Selected"),
                         fill = c("lightgreen","darkgreen"))

overlap3 <- calculate.overlap(x=list("K27ac_Selected"=unique(genes2),
                                     "K4me3_Selected"=unique(genes4)))
v3 <- draw.pairwise.venn(length(overlap3$a1),length(overlap3$a2),length(overlap3$a3),
                         c("K27ac_Selected","K4me3_Selected"),
                         fill = c("darkblue","darkgreen"))