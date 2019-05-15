#!/usr/bin/env Rscript
##Example of running command:
###while read a b c d e f g h i j; do bsub -q normal -G teamtrynka -M100 -R'select[mem>100] rusage[mem=100]' -o testing.%J Rscript Testing_Arguments.R "$a" "$b" "$c" "$d" "$e" "$f" "$g" "$h" "$i" "$j"; done < COLOC_loci_assigned_4PlottingScript.txt

args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(tidyr)
library(ggplot2)
library(sva)
library(limma) 
library(reshape2)
#devtools::install_github("kauralasoo/wiggleplotr")
library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)
library(biomaRt)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)
library(Gviz)
library(gridExtra)

#Helper function to make wiggle plots

readCoverageFromBigWig <- function(bigwig_path, gene_range){
  #Read coverage over a region from a bigWig file
  sel = rtracklayer::BigWigSelection(gene_range)
  coverage_ranges = rtracklayer::import.bw(bigwig_path, selection = sel)
  GenomeInfoDb::seqlevels(coverage_ranges) = S4Vectors::as.vector.Rle(GenomicRanges::seqnames(gene_range), mode = "character")
  coverage_rle = GenomicRanges::coverage(coverage_ranges, weight = GenomicRanges::score(coverage_ranges))[[1]]
  coverage_rle = coverage_rle[(GenomicRanges::start(gene_range)):(GenomicRanges::end(gene_range))] #Keep the region of interest
}

joinExons <- function(exons) {
  #Join a list of exons into one GRanges object
  
  #Test that all transcripts are on the same chromosome
  chrs = purrr::map_chr(as.list(exons), ~GenomicRanges::seqnames(.)[1] %>% 
                          S4Vectors::as.vector.Rle(mode = "character"))
  if (!all(chrs == chrs[1])){
    stop("Some transcripts are on different chromosomes.")
  }
  
  #Join all exons together
  transcript_ids = names(exons)
  joint_exons = c()
  for(tx_id in transcript_ids){
    tx = exons[[tx_id]]
    if(length(joint_exons) == 0){
      joint_exons = tx
    }
    else{
      joint_exons = c(joint_exons, tx)
    }
  }
  joint_exons = GenomicRanges::reduce(joint_exons)
  return(joint_exons)
}

extractStrandsFromGrangesList <- function(granges_list){
  strands = purrr::map(as.list(granges_list), ~(GenomicRanges::strand(.) %>%
                                                  S4Vectors::as.vector.Rle(.,"character"))[1])
  return(unlist(strands))
}

prepareTranscriptAnnotations <- function(transcript_annotations){
  assertthat::assert_that(assertthat::has_name(transcript_annotations, "transcript_id"))
  assertthat::assert_that(assertthat::has_name(transcript_annotations, "strand"))
  
  
  #Make sure that the strand information is represented correctly
  transcript_annotations = dplyr::mutate(transcript_annotations,
                                         strand = ifelse(strand %in% c("+","*") | strand == 1, 1, -1))
  
  #Add transcript label
  if(assertthat::has_name(transcript_annotations, "gene_name")){
    transcript_annotations = dplyr::select_(transcript_annotations, "transcript_id", "gene_name", "strand") %>% 
      dplyr::mutate(transcript_label = ifelse(strand == 1, 
                                              paste(paste(gene_name, transcript_id, sep = ":")," >",sep =""), 
                                              paste("< ",paste(gene_name, transcript_id, sep = ":"),sep ="")))
  } else{
    transcript_annotations = dplyr::mutate(transcript_annotations, transcript_label = ifelse(strand == 1, 
                                                                                             paste(paste(transcript_id, sep = ":")," >",sep =""), 
                                                                                             paste("< ",paste(transcript_id, sep = ":"),sep =""))) 
  }
  return(transcript_annotations)
}

prepareTranscriptStructureForPlotting <- function(exon_ranges, cds_ranges, transcript_annotations){
  #Combine exon_ranges and cds_ranges into a single data.frame that also contains transcript rank
  
  #Convert exon ranges into data.frame and add transcript rank
  exons_df = purrr::map_df(exon_ranges, data.frame, .id = "transcript_id")
  exons_df = dplyr::mutate(exons_df, transcript_rank = as.numeric(factor(exons_df$transcript_id)), type = "")
  transcript_rank = unique(exons_df[,c("transcript_id", "transcript_rank", "type")])
  
  #Convert CDS ranges into a data.frame
  cds_df = purrr::map_df(cds_ranges, data.frame, .id = "transcript_id")
  cds_df = dplyr::left_join(cds_df, transcript_rank, by = "transcript_id") #Add matching transcript rank
  
  #Join exons and cdss together
  exons_df = dplyr::mutate(exons_df, feature_type = "exon")
  cds_df = dplyr::mutate(cds_df, feature_type = "cds")
  transcript_struct = rbind(exons_df, cds_df)
  
  #Add transcript label to transcript structure
  transcript_struct = dplyr::left_join(transcript_struct, transcript_annotations, by = "transcript_id")
  return(transcript_struct)
}

intronsFromJointExonRanges <- function(joint_exon_ranges, flanking_length){
  #Construct intron ranges from joint exon ranges
  introns = IRanges::gaps(joint_exon_ranges, 
                          start = min(IRanges::start(joint_exon_ranges)) - flanking_length[1], 
                          end = max(IRanges::end(joint_exon_ranges)) + flanking_length[2])
  return(introns)
}

# Find the start and end coordinates of the whole gene form joint exons. 
constructGeneRange <- function(joint_exon_ranges, flanking_length){
  gene_range = GenomicRanges::reduce(c(joint_exon_ranges, GenomicRanges::gaps(joint_exon_ranges, start = NA, end = NA)))
  GenomeInfoDb::seqlevels(gene_range) = S4Vectors::as.vector.Rle(GenomicRanges::seqnames(gene_range), mode = "character")[1]
  GenomicRanges::start(gene_range) = GenomicRanges::start(gene_range) - flanking_length[1]
  GenomicRanges::end(gene_range) = GenomicRanges::end(gene_range) + flanking_length[2]
  return(gene_range)
}

#' Paste two factors together and preserved their joint order.
#'
#' @param factor1 First factor
#' @param factor2 Second factor
#' 
#' @return Factors factor1 and factor2 pasted together.
pasteFactors <- function(factor1, factor2){
  #Extract levels
  levels1 = levels(factor1)
  levels2 = levels(factor2)
  
  #Construct joint levels
  new1 = rep(levels1, length(levels2))
  new2 = rep(levels2, each = length(levels1))
  new_levels = paste(new1, new2, sep = "_")
  
  new_factor = factor(paste(as.character(factor1), as.character(factor2), sep = "_"), levels = new_levels)
  return(new_factor)
}

# Calculate mean coverage within each track_id and colour_group
meanCoverage <- function(coverage_df){
  coverage_df = dplyr::group_by_(coverage_df, "track_id", "colour_group", "bins") %>% 
    dplyr::summarise_(.dots = stats::setNames(list(~mean(coverage)), c("coverage"))) %>%
    dplyr::ungroup() %>% # It's important to do ungroup before mutate, or you get unexpected factor results
    dplyr::mutate_(.dots = stats::setNames(list(~pasteFactors(as.factor(track_id), as.factor(colour_group))),c("sample_id")) ) #Construct a new sample id for mean vector
  return(coverage_df)
}

# Choose a subsample of points to make plotting faster
# Makes sure that intron-exon boundaries are well samples.
subsamplePoints <- function(tx_annotations, plot_fraction){
  #Define the start and end coorinates of the region
  region_start = min(IRanges::start(tx_annotations$new_introns))
  region_end = max(IRanges::end(tx_annotations$new_introns))
  region_length = region_end - region_start
  
  #Take a subsample of points that's easier to plot
  points = sample(region_length, floor(region_length*plot_fraction))
  #Subtract the start coordinate of the region
  exon_starts = unique(unlist(lapply(tx_annotations$exon_ranges, IRanges::start))) - (region_start -1)
  exon_ends = unique(unlist(lapply(tx_annotations$exon_ranges, IRanges::end))) - (region_start - 1)
  points = unique(sort(c(points, exon_starts, exon_ends, 
                         exon_starts -3, exon_starts +3, 
                         exon_ends + 3, exon_ends -3)))
  points = points[points >= 0]
  return(points)
}


#' Returns a three-colour palette suitable for visualising read coverage stratified by genotype
#'
#' @return Vector of three colours.
#' @export
#'
#' @examples
#' getGenotypePalette()
getGenotypePalette <- function(old = FALSE){
  if(old){
    c("#d7191c","#fdae61","#1a9641")
  } else{
    #Borrowed from Kumasaka, et al 2015
    c("#E9181D","#51BEEE","#18354B") 
  }
}

#Common theme for all data track plots
dataTrackTheme <- function(){
  theme = theme(axis.text.x = element_blank(), 
                axis.title.x = element_blank(), 
                axis.ticks.x = element_blank(),
                plot.margin=unit(c(0.1,1,0.1,1),"line"),
                legend.position="none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.text.y = element_text(colour = "grey10"),
                strip.background = element_rect(fill = "grey85"))
  return(theme)
}

shrinkIntronsCoverage <- function(coverage, old_introns, new_introns){
  
  #Covert coverage vector from Rle to normal vector
  coverage = S4Vectors::as.vector.Rle(coverage, mode = "double")
  
  #Calculate full annotations
  old_annot = sort(c(old_introns, IRanges::gaps(old_introns)))
  new_annot = sort(c(new_introns, IRanges::gaps(new_introns)))
  
  #If new and old annotations are identical then return coverage as data frame
  if(all(IRanges::width(old_annot) == IRanges::width(new_annot))){
    bins = seq(min(IRanges::start(new_annot)), max(IRanges::end(new_annot)))
    #Make sure that coverage vector and bins vector have equal length
    assertthat::assert_that(assertthat::are_equal(length(bins), length(coverage)))
    new_coverage = dplyr::data_frame(bins = bins, coverage = coverage)
    return(new_coverage)
    
  } else{ #Otherwise shrink intron converage
    
    #Calculate the width of each annotation bin
    bin_width = ceiling(IRanges::width(old_annot)/IRanges::width(new_annot))
    #Build summarisation groups
    s_coord = IRanges::start(new_annot)
    e_coord = IRanges::end(new_annot)
    w_old = IRanges::width(old_annot)
    
    bins = c()
    
    for (i in seq_along(new_annot)){
      bin_id = rep(c(s_coord[i]:e_coord[i]),each = bin_width[i])[1:w_old[i]]
      bins = c(bins, bin_id)
    }
    
    #Calculate mean coverage in bins
    df = data.frame(coverage, bins)
    new_coverage = dplyr::summarize(dplyr::group_by(df, bins), coverage = mean(coverage))
    return(new_coverage)
  }
}

translateExonCoordinates <- function(exons, old_introns, new_introns){
  #Tranlate exon coordinates by shortening introns
  old_exon_starts = IRanges::start(exons)
  old_intron_ends = IRanges::end(old_introns)
  new_intron_ends = IRanges::end(new_introns)
  
  #Translate old exon coordinates to new exon coordinates
  new_exon_starts = rep(0,length(old_exon_starts))
  for (i in seq_along(old_exon_starts)){
    #Find the nearest upstream intron for the current gene
    nearest_intron_number = max(which(old_exon_starts[i] > old_intron_ends))
    new_exon_starts[i] = old_exon_starts[i] - old_intron_ends[nearest_intron_number] + new_intron_ends[nearest_intron_number]
  }
  
  #Create new exon coordinates
  new_exons = IRanges::IRanges(start = new_exon_starts, width = IRanges::width(exons))
  return(new_exons)
}

rescaleIntrons <- function(exons, cdss, joint_exons, new_intron_length, flanking_length){
  
  #Convert exons and cds objects to ranges
  exon_ranges = lapply(exons, GenomicRanges::ranges)
  cds_ranges = lapply(cdss, GenomicRanges::ranges)
  
  #Shorten introns and translate exons into the new exons
  old_introns = intronsFromJointExonRanges(GenomicRanges::ranges(joint_exons), flanking_length = flanking_length)
  new_introns = shortenIntrons(old_introns,new_intron_length)
  new_exon_ranges = lapply(exon_ranges, translateExonCoordinates, old_introns, new_introns)
  new_cds_ranges = lapply(cds_ranges, translateExonCoordinates, old_introns, new_introns)
  
  return(list(exon_ranges = new_exon_ranges, cds_ranges = new_cds_ranges, 
              old_introns = old_introns, new_introns = new_introns))
}

plotTranscriptStructure <- function(exons_df, limits = NA, connect_exons = TRUE,  
                                    xlabel = "Distance from gene start (bp)", transcript_label = TRUE){
  
  #Extract the position for plotting transcript name
  transcript_annot = dplyr::group_by_(exons_df, ~transcript_id) %>% 
    dplyr::filter_(~feature_type == "exon") %>%
    dplyr::arrange_('transcript_id', 'start') %>%
    dplyr::filter(row_number() == 1)
  
  #Create a plot of transcript structure
  plot = ggplot(exons_df) + geom_blank()
  if(connect_exons){ #Print line connecting exons
    plot = plot + geom_line(aes_(x = ~start, y = ~transcript_rank, group = ~transcript_rank, color = ~feature_type))
  }
  plot = plot + 
    geom_rect(aes_(xmin = ~start, 
                   xmax = ~end, 
                   ymax = ~transcript_rank + 0.25, 
                   ymin = ~transcript_rank - 0.25, 
                   fill = ~feature_type)) + 
    theme_light() +
    theme(plot.margin=unit(c(0,1,1,1),"line"), 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(colour = "grey10"),
          strip.background = element_rect(fill = "grey85")) +
    xlab(xlabel) +
    facet_grid(type~.) +
    scale_y_continuous(expand = c(0.2,0.15)) +
    scale_fill_manual(values = c("#2c7bb6","#abd9e9")) + 
    scale_colour_manual(values = c("#2c7bb6","#abd9e9"))
  if(all(!is.na(limits))){
    plot = plot + scale_x_continuous(limits = limits, expand = c(0,0))
  }
  if(transcript_label){
    plot = plot + geom_text(aes_(x = ~start, 
                                 y = ~transcript_rank + 0.30, 
                                 label = ~transcript_label), 
                            data = transcript_annot, hjust = 0, vjust = 0, size = 4)
    
  }
  return(plot)
}

makeCoveragePlot2 <- function(coverage_df, limits, alpha, fill_palette, coverage_type){
  #Plot coverage over a region
  coverage_plot = ggplot(coverage_df, aes_(~bins, ~coverage, group = ~sample_id, alpha = ~alpha)) + 
    geom_blank() +
    theme_light()
  #Choose between plotting a line and plotting area
  if(coverage_type == "line"){
    coverage_plot = coverage_plot + 
      geom_line(aes_(colour = ~colour_group), alpha = alpha, position = "identity") 
  } else if (coverage_type == "area"){
    coverage_plot = coverage_plot + 
      geom_area(aes_(fill = ~colour_group), alpha = alpha, position = "identity")
  } else if (coverage_type == "both"){
    coverage_plot = coverage_plot + 
      geom_area(aes_(fill = ~colour_group), alpha = alpha, position = "identity") +
      geom_line(aes_(colour = ~colour_group), alpha = alpha, position = "identity") 
  } else{
    stop("Coverage type not supported.")
  }
  coverage_plot = coverage_plot +
    facet_grid(track_id~.) +
    dataTrackTheme() + 
    scale_x_continuous(limits = limits, expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_color_manual(values = fill_palette) +
    scale_fill_manual(values = fill_palette) +
    ylab("SPMR")
  return(coverage_plot)
}


plotCoverage2 <- function(exons, cdss = NULL, transcript_annotations = NULL, track_data, rescale_introns = TRUE,
                          new_intron_length = 50, flanking_length = c(50,50),
                          plot_fraction = 0.1, heights = c(0.75, 0.25), alpha = 0.7,
                          fill_palette = c("#a1dab4","#41b6c4","#225ea8"), mean_only = TRUE, 
                          connect_exons = TRUE, transcript_label = TRUE, return_subplots_list = FALSE,
                          region_coords = NULL, coverage_type = "area"){
  
  #IF cdss is not specified then use exons instead on cdss
  if(is.null(cdss)){
    cdss = exons
  }
  
  #Make some assertions about the input data
  #Check track_data
  assertthat::assert_that(assertthat::has_name(track_data, "sample_id"))
  assertthat::assert_that(assertthat::has_name(track_data, "track_id"))
  assertthat::assert_that(assertthat::has_name(track_data, "bigWig"))
  assertthat::assert_that(assertthat::has_name(track_data, "scaling_factor"))
  assertthat::assert_that(assertthat::has_name(track_data, "colour_group"))
  
  #Make sure that bigWig column is not a factor
  if(is.factor(track_data$bigWig)){
    warning("bigWig column in track_data data.frame is a factor, coverting to a character vector.")
    track_data = dplyr::mutate_(track_data, .dots = stats::setNames(list(~as.character(bigWig)), c("bigWig")))
  }
  
  #Check transcript annotation
  #If transcript annotations are not supplied then construct them manually from the GRanges list
  if(is.null(transcript_annotations)){
    plotting_annotations = dplyr::data_frame(transcript_id = names(exons),
                                             strand = extractStrandsFromGrangesList(exons)) %>%
      prepareTranscriptAnnotations()
  } else{
    assertthat::assert_that(assertthat::has_name(transcript_annotations, "transcript_id"))
    assertthat::assert_that(assertthat::has_name(transcript_annotations, "gene_name"))
    assertthat::assert_that(assertthat::has_name(transcript_annotations, "strand"))
    plotting_annotations = prepareTranscriptAnnotations(transcript_annotations)
  }
  
  #Check exons and cdss
  assertthat::assert_that(is.list(exons) || class(exons) == "GRangesList") #Check that exons and cdss objects are lists
  assertthat::assert_that(is.list(cdss) || class(exons) == "GRangesList")
  #TODO: Check that the names of the exons and cdss list match that of the transcript_annotations data.frame
  
  #Find the start and end cooridinates of the whole region spanning the gene
  joint_exons = joinExons(exons)
  
  #If region_coords is specificed, then ignore the flanking_length attrbute and compute
  # flanking_length form region_coords
  if(!is.null(region_coords)){
    gene_range = constructGeneRange(joint_exons, c(0,0))
    min_start = min(GenomicRanges::start(gene_range))
    max_end = max(GenomicRanges::end(gene_range))
    flanking_length = c(min_start - region_coords[1], region_coords[2] - max_end)
    
    gene_range = constructGeneRange(joint_exons, flanking_length)
  } else{
    gene_range = constructGeneRange(joint_exons, flanking_length)
  }
  assertthat::assert_that(length(flanking_length) == 2) #flanking_length is a vector of two elements
  
  #Extract chromosome name
  chromosome_name = as.vector(GenomicRanges::seqnames(gene_range)[1])
  
  #Read coverage tracks from BigWig file
  sample_list = as.list(track_data$bigWig)
  names(sample_list) = track_data$sample_id
  coverage_list = lapply(sample_list, readCoverageFromBigWig, gene_range)
  
  #Shorten introns and translate exons into the new introns
  if(rescale_introns){
    #Recale transcript annotations
    tx_annotations = rescaleIntrons(exons, cdss, joint_exons, 
                                    new_intron_length = new_intron_length, flanking_length = flanking_length)
    #Make a label for gene structure plot
    xlabel = "Distance from region start (bp)"
  }
  else{ #Do not rescale transcript annotationn
    #Need to calculate joint intron coordinates for transcript annotations
    old_introns = intronsFromJointExonRanges(GenomicRanges::ranges(joint_exons), flanking_length = flanking_length)
    tx_annotations = list(exon_ranges = lapply(exons, GenomicRanges::ranges), cds_ranges = lapply(cdss, GenomicRanges::ranges),
                          old_introns = old_introns, new_introns = old_introns)
    
    #Make a label for gene structure plot
    xlabel = paste("Chromosome", chromosome_name, "position (bp)")
  }
  #Shrink intron coverage and convert coverage vectors into data frames
  coverage_list = lapply(coverage_list, shrinkIntronsCoverage, tx_annotations$old_introns, tx_annotations$new_introns)
  
  #Take a subsample of points that is easier to plot
  points = subsamplePoints(tx_annotations, plot_fraction)
  coverage_list = lapply(coverage_list, function(x) {x[points,]} )
  
  #Convert to data frame and plot
  coverage_df = purrr::map_df(coverage_list, identity, .id = "sample_id") %>% 
    as.data.frame() %>%
    dplyr::mutate_(.dots = stats::setNames(list(~as.character(sample_id)), c("sample_id")) ) #Convert factor to character
  coverage_df = dplyr::left_join(coverage_df, track_data, by = "sample_id") %>%
    dplyr::mutate_(.dots = stats::setNames(list(~coverage/scaling_factor), c("coverage")) ) #Normalize by library size
  
  #Calculate mean coverage within each track and colour group
  if(mean_only){  coverage_df = meanCoverage(coverage_df) }
  
  #Make plots
  #Construct transcript structure data.frame from ranges lists
  limits = c( min(IRanges::start(tx_annotations$new_introns)), max(IRanges::end(tx_annotations$new_introns)))
  transcript_struct = prepareTranscriptStructureForPlotting(tx_annotations$exon_ranges, 
                                                            tx_annotations$cds_ranges, plotting_annotations)
  tx_structure = plotTranscriptStructure(transcript_struct, limits, connect_exons, xlabel, transcript_label)
  
  coverage_plot = makeCoveragePlot2(coverage_df, limits, alpha, fill_palette, coverage_type)
  
  #Choose between returning plot list or a joint plot using plot_grid
  if(return_subplots_list){
    plot_list = list(coverage_plot = coverage_plot, tx_structure = tx_structure)
    return(c(plot_list,coverage_df))
  } else {
    plot = cowplot::plot_grid(coverage_plot, tx_structure, align = "v", rel_heights = heights, ncol = 1)
    return(plot)
  }
}

###Accepting arguments

args = commandArgs(trailingOnly=TRUE)
window = 250000
chr = args[1]
snp = args[2]
snp_bp  = as.numeric(args[3])
GWAS  = args[4]
slope  = args[5]
atacpeak = args[6]
k27acpeak = args[7]
k4me3peak = args[8]
gene = args[9]
gene_name = args[10]

###Extracting regions

if(k27acpeak != "chrNA"){
  k27acpeak_region <- as.vector(unlist(strsplit(k27acpeak, "_")))
}
if(k4me3peak != "chrNA"){
  k4me3peak_region <- as.vector(unlist(strsplit(k4me3peak, "_")))
}
if(atacpeak != "chrNA"){
  atacpeak_region <- as.vector(unlist(strsplit(atacpeak, "_")))
}
dir.create(gene_name)

## Load SNP file
SNP_file_name <- paste("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/TitanTable/AllDonors_Imputed_AllAutosomes_MAF0.1_GRCh38_Setted_",chr,".vcf", sep="")
snp_table = read.table(SNP_file_name,sep="\t", header=T,row.names=3, stringsAsFactors=FALSE)

## Load gene expression file
expression_file_name <- c("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/Variance_Decomposition/day0.norm_prop_CQN_16PCsRegressed.txt")
expr = read.table(expression_file_name,sep="\t", header=T,row.names=1)

##Violin plots if slope is possitive
if(slope > 0){
  if(gene != "Peak"){
## plot your eQTL Violing plot
    e = subset(expr,rownames(expr)== gene);
    s = subset(snp_table,rownames(snp_table)==snp)
    variant <- as.vector(unlist(strsplit(rownames(s), "_")))
    s[s == "1|1"] <- paste(variant[4],variant[4],sep = "-")
    s[s == "0|1"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "1|0"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "0|0"] <- paste(variant[3],variant[3],sep = "-")


    mysamples <- colnames(e) 
    s <- s[mysamples]
    es=as.data.frame(t(rbind(e,s)))
    colnames(es) <- c("gene","snp")
    es$gene<-as.numeric(as.character(es$gene))
    es$snp<-as.factor(es$snp)
    es$snp <- factor(es$snp, levels = c(paste(variant[3],variant[3],sep = "-"),paste(variant[3],variant[4],sep = "-"),paste(variant[4],variant[4],sep = "-")))


    pals <- rev(c("#deaff2","#d387f4","#d46ff9"))
    eQTL_plot <- ggplot(es,aes(x=snp,y=gene, color= snp)) +
      geom_violin(draw_quantiles = c(0.5), position = "dodge") +
      geom_boxplot(width=.1, outlier.colour=NA, position = "dodge") +
      geom_jitter(position=position_jitter(width=0.1), size=1.5) +
      xlab(paste(chr,snp_bp,sep=":")) + ylab(paste("norm CQN ", gene_name )) +
      scale_color_manual(values = pals) +
      theme_minimal(base_size=14) +
      theme(legend.title=element_blank(),panel.grid.minor = element_blank(),legend.position="none",panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line( size=.05, color="#D9D9D9"),legend.text = element_text(size = 12),
            strip.background = element_blank(),
            strip.text.x = element_blank())
    ggsave(paste(gene_name, "/treg_eQTL_",gene_name, "_", gene, "_", snp, "_Violin.pdf", sep =""), plot = eQTL_plot, width = 7, height = 6, units = "cm",
         dpi = 300,device = "pdf", useDingbats=FALSE)


    Major_eQTL <- es[es$snp == paste(variant[3],variant[3],sep = "-"),]
    Major_eQTL <- head(Major_eQTL[order(Major_eQTL$gene),],5)
    Minor_eQTL <- es[es$snp == paste(variant[4],variant[4],sep = "-"),]
    Minor_eQTL <- head(Minor_eQTL[order(-Minor_eQTL$gene),],5)
    heteroz_eQTL <- es[es$snp == paste(variant[3],variant[4],sep = "-"),]
    heteroz_eQTL <- heteroz_eQTL[order(-heteroz_eQTL$gene),]
    heteroz_mean_eQTL <- heteroz_eQTL[(ceiling(length(heteroz_eQTL$gene)/2)-2):(ceiling(length(heteroz_eQTL$gene)/2)+2),]

    print(Major_eQTL)
    print(heteroz_mean_eQTL)
    print(Minor_eQTL)
}
## Load K27ac file
  if(k27acpeak != "chrNA"){
    k27ac_file_name <- c("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/Variance_Decomposition/ChM_K27ac_Peak_featureCounts_NoBatchCor_CQN_22PCsRegressed.txt")
    k27ac = read.table(k27ac_file_name,sep="\t", header=T,row.names=4)
    k27ac <- k27ac[,-(1:3)] 
    k27ac <- k27ac[,-(1:2)] 

## plot your K27ac peak and SNP of interest Violing plot
    k27ac_e = subset(k27ac,rownames(k27ac)== k27acpeak);
    s = subset(snp_table,rownames(snp_table)==snp)
    variant <- as.vector(unlist(strsplit(rownames(s), "_")))
    s[s == "1|1"] <- paste(variant[4],variant[4],sep = "-")
    s[s == "0|1"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "1|0"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "0|0"] <- paste(variant[3],variant[3],sep = "-")


    mysamples <- colnames(k27ac_e) 
    s <- s[mysamples]
    rm(es)
    es=as.data.frame(t(rbind(k27ac_e,s)))
    colnames(es) <- c("gene","snp")
    es$gene<-as.numeric(as.character(es$gene))
    es$snp<-as.factor(es$snp)
    es$snp <- factor(es$snp, levels = c(paste(variant[3],variant[3],sep = "-"),paste(variant[3],variant[4],sep = "-"),paste(variant[4],variant[4],sep = "-")))


    pals <- rev(c("#9df4e9","#85d3ca","#78bbb5"))
    k27acQTL_plot <- ggplot(es,aes(x=snp,y=gene, color= snp)) +
      geom_violin(draw_quantiles = c(0.5), position = "dodge") +
      geom_boxplot(width=.1, outlier.colour=NA, position = "dodge") +
      geom_jitter(position=position_jitter(width=0.1), size=1.5) +
      xlab(paste(chr,snp_bp,sep=":")) + ylab(paste(k27acpeak," norm CQN",sep="")) +
      scale_color_manual(values = pals) +
      theme_minimal(base_size=14) +
      theme(legend.title=element_blank(),panel.grid.minor = element_blank(),legend.position="none",panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line( size=.05, color="#D9D9D9"),legend.text = element_text(size = 14),
            strip.background = element_blank(),
            strip.text.x = element_blank()) 
    ggsave(paste(gene_name,"/treg_actQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = k27acQTL_plot, width = 7, height = 6, units = "cm",
          dpi = 300,device = "pdf", useDingbats=FALSE)


    Major_actQTL <- es[es$snp == paste(variant[3],variant[3],sep = "-"),]
    Major_actQTL <- head(Major_actQTL[order(Major_actQTL$gene),],5)
    Minor_actQTL <- es[es$snp == paste(variant[4],variant[4],sep = "-"),]
    Minor_actQTL <- head(Minor_actQTL[order(-Minor_actQTL$gene),],5)
    heteroz_actQTL <- es[es$snp == paste(variant[3],variant[4],sep = "-"),]
    heteroz_actQTL <- heteroz_actQTL[order(-heteroz_actQTL$gene),]
    heteroz_actQTL_mean <- heteroz_actQTL[(ceiling(length(heteroz_actQTL$gene)/2)-2):(ceiling(length(heteroz_actQTL$gene)/2)+2),]

    print(Major_actQTL)
    print(heteroz_actQTL_mean)
    print(Minor_actQTL)
  }
## Load K4me3 file
  if(k4me3peak != "chrNA"){
    k4me3_file_name <- c("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/Variance_Decomposition/ChM_K4me3_Peak_featureCounts_NoBatchCor_CQN_33PCsRegressed.txt")
    k4me3 = read.table(k4me3_file_name,sep="\t", header=T,row.names=4)
    k4me3 <- k4me3[,-(1:3)] 
    k4me3 <- k4me3[,-(1:2)] 

## plot your K4me3 peak and SNP of interest Violing plot
    k4me3_e = subset(k4me3,rownames(k4me3)== k4me3peak);
    s = subset(snp_table,rownames(snp_table)==snp)
    variant <- as.vector(unlist(strsplit(rownames(s), "_")))
    s[s == "1|1"] <- paste(variant[4],variant[4],sep = "-")
    s[s == "0|1"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "1|0"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "0|0"] <- paste(variant[3],variant[3],sep = "-")


    mysamples <- colnames(k4me3_e) 
    s <- s[mysamples]
    es=as.data.frame(t(rbind(k4me3_e,s)))
    colnames(es) <- c("gene","snp")
    es$gene<-as.numeric(as.character(es$gene))
    es$snp<-as.factor(es$snp)
    es$snp <- factor(es$snp, levels = c(paste(variant[3],variant[3],sep = "-"),paste(variant[3],variant[4],sep = "-"),paste(variant[4],variant[4],sep = "-")))


    pals <- rev(c("#f4adab", "#ef877f", "#f47362"))
    k4me3QTL_plot <- ggplot(es,aes(x=snp,y=gene, color= snp)) +
      geom_violin(draw_quantiles = c(0.5), position = "dodge") +
      geom_boxplot(width=.1, outlier.colour=NA, position = "dodge") +
      geom_jitter(position=position_jitter(width=0.1), size=1.5) +
      xlab(paste(chr,snp_bp,sep=":")) + ylab(paste(k4me3peak," norm CQN",sep="")) +
      scale_color_manual(values = pals) +
      theme_minimal(base_size=14) +
      theme(legend.title=element_blank(),panel.grid.minor = element_blank(),legend.position="none",panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line( size=.05, color="#D9D9D9"),legend.text = element_text(size = 14),
            strip.background = element_blank(),
            strip.text.x = element_blank()) +
      ylim(-1.3,1.3)
    ggsave(paste(gene_name, "/treg_K4me3QTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = k4me3QTL_plot, width = 7, height = 6, units = "cm",
          dpi = 300,device = "pdf", useDingbats=FALSE)


    Major_k4me3QTL <- es[es$snp == paste(variant[3],variant[3],sep = "-"),]
    Major_k4me3QTL <- head(Major_k4me3QTL[order(Major_k4me3QTL$gene),],5)
    Minor_k4me3QTL <- es[es$snp == paste(variant[4],variant[4],sep = "-"),]
    Minor_k4me3QTL <- head(Minor_k4me3QTL[order(-Minor_k4me3QTL$gene),],5)
    heteroz_k4me3QTL <- es[es$snp == paste(variant[3],variant[4],sep = "-"),]
    heteroz_k4me3QTL <- heteroz_k4me3QTL[order(-heteroz_k4me3QTL$gene),]
    heteroz_k4me3QTL_mean <- heteroz_k4me3QTL[(ceiling(length(heteroz_k4me3QTL$gene)/2)-2):(ceiling(length(heteroz_k4me3QTL$gene)/2)+2),]

    print(Major_k4me3QTL)
    print(heteroz_k4me3QTL_mean)
    print(Minor_k4me3QTL)
  }
## Load ATAC file
  if(atacpeak != "chrNA"){
    atac_file_name <- c("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/Variance_Decomposition/ATAC_Counts_featureCounts_CQN_30PCsRegressed.txt")
    atac = read.table(atac_file_name,sep="\t", header=T,row.names=4)
    atac <- atac[,-(1:3)] 
    atac <- atac[,-(1:2)] 

## plot your ATAC peak and SNP of interest Violing plot
    atac_e = subset(atac,rownames(atac)== atacpeak);
    s = subset(snp_table,rownames(snp_table)==snp)
    variant <- as.vector(unlist(strsplit(rownames(s), "_")))
    s[s == "1|1"] <- paste(variant[4],variant[4],sep = "-")
    s[s == "0|1"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "1|0"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "0|0"] <- paste(variant[3],variant[3],sep = "-")


    mysamples <- colnames(atac_e) 
    s <- s[mysamples]
    es=as.data.frame(t(rbind(atac_e,s)))
    colnames(es) <- c("gene","snp")
    es$gene<-as.numeric(as.character(es$gene))
    es$snp<-as.factor(es$snp)
    es$snp <- factor(es$snp, levels = c(paste(variant[3],variant[3],sep = "-"),paste(variant[3],variant[4],sep = "-"),paste(variant[4],variant[4],sep = "-")))


    pals <- rev(c("#efd6c5", "#b78b69", "#bb804f"))
    atacQTL_plot <- ggplot(es,aes(x=snp,y=gene, color= snp)) +
      geom_violin(draw_quantiles = c(0.5), position = "dodge") +
      geom_boxplot(width=.1, outlier.colour=NA, position = "dodge") +
      geom_jitter(position=position_jitter(width=0.1), size=1.5) +
      xlab(paste(chr,snp_bp,sep=":")) + ylab(paste(atacpeak," norm CQN",sep="")) +
      scale_color_manual(values = pals) +
      theme_minimal(base_size=14) +
      theme(legend.title=element_blank(),panel.grid.minor = element_blank(),legend.position="none",panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line( size=.05, color="#D9D9D9"),legend.text = element_text(size = 14),
            strip.background = element_blank(),
            strip.text.x = element_blank())
    ggsave(paste(gene_name, "/treg_ATACQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = atacQTL_plot, width = 7, height = 6, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)

    Major_atacQTL <- es[es$snp == paste(variant[3],variant[3],sep = "-"),]
    Major_atacQTL <- head(Major_atacQTL[order(Major_atacQTL$gene),],5)
    Minor_atacQTL <- es[es$snp == paste(variant[4],variant[4],sep = "-"),]
    Minor_atacQTL <- head(Minor_atacQTL[order(-Minor_atacQTL$gene),],5)
    heteroz_atacQTL <- es[es$snp == paste(variant[3],variant[4],sep = "-"),]
    heteroz_atacQTL <- heteroz_atacQTL[order(-heteroz_atacQTL$gene),]
    heteroz_atacQTL_mean <- heteroz_atacQTL[(ceiling(length(heteroz_atacQTL$gene)/2)-2):(ceiling(length(heteroz_atacQTL$gene)/2)+2),]

    print(Major_atacQTL)
    print(heteroz_atacQTL_mean)
    print(Minor_atacQTL)
  }
}

##Violin plots if slope is negative
if(slope < 0){
  if(gene != "Peak"){
  ## plot your eQTL Violing plot
    e = subset(expr,rownames(expr)== gene);
    s = subset(snp_table,rownames(snp_table)==snp)
    variant <- as.vector(unlist(strsplit(rownames(s), "_")))
    s[s == "1|1"] <- paste(variant[4],variant[4],sep = "-")
    s[s == "0|1"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "1|0"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "0|0"] <- paste(variant[3],variant[3],sep = "-") 
  
    mysamples <- colnames(e) 
    s <- s[mysamples]
    es=as.data.frame(t(rbind(e,s)))
    colnames(es) <- c("gene","snp")
    es$gene<-as.numeric(as.character(es$gene))
    es$snp<-as.factor(es$snp)
    es$snp <- factor(es$snp, levels = c(paste(variant[3],variant[3],sep = "-"),paste(variant[3],variant[4],sep = "-"),paste(variant[4],variant[4],sep = "-")))
  
  
    pals <- rev(c("#deaff2","#d387f4","#d46ff9"))
    eQTL_plot <- ggplot(es,aes(x=snp,y=gene, color= snp)) +
      geom_violin(draw_quantiles = c(0.5), position = "dodge") +
      geom_boxplot(width=.1, outlier.colour=NA, position = "dodge") +
      geom_jitter(position=position_jitter(width=0.1), size=1.5) +
      xlab(paste(chr,snp_bp,sep=":")) + ylab(paste("norm CQN ", gene_name )) +
      scale_color_manual(values = pals) +
      theme_minimal(base_size=14) +
      theme(legend.title=element_blank(),panel.grid.minor = element_blank(),legend.position="none",panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line( size=.05, color="#D9D9D9"),legend.text = element_text(size = 12),
            strip.background = element_blank(),
            strip.text.x = element_blank())
    ggsave(paste(gene_name, "/treg_eQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = eQTL_plot, width = 7, height = 6, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
  
  
  Major_eQTL <- es[es$snp == paste(variant[3],variant[3],sep = "-"),]
  Major_eQTL <- head(Major_eQTL[order(-Major_eQTL$gene),],5)
  Minor_eQTL <- es[es$snp == paste(variant[4],variant[4],sep = "-"),]
  Minor_eQTL <- head(Minor_eQTL[order(Minor_eQTL$gene),],5)
  heteroz_eQTL <- es[es$snp == paste(variant[3],variant[4],sep = "-"),]
  heteroz_eQTL <- heteroz_eQTL[order(-heteroz_eQTL$gene),]
  heteroz_mean_eQTL <- heteroz_eQTL[(ceiling(length(heteroz_eQTL$gene)/2)-2):(ceiling(length(heteroz_eQTL$gene)/2)+2),]
  
  print(Major_eQTL)
  print(heteroz_mean_eQTL)
  print(Minor_eQTL)
}
  
  ## Load K27ac file
  if(k27acpeak != "chrNA"){
    k27ac_file_name <- c("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/Variance_Decomposition/ChM_K27ac_Peak_featureCounts_NoBatchCor_CQN_22PCsRegressed.txt")
    k27ac = read.table(k27ac_file_name,sep="\t", header=T,row.names=4)
    k27ac <- k27ac[,-(1:3)] 
    k27ac <- k27ac[,-(1:2)] 
  
    ## plot your K27ac peak and SNP of interest Violing plot
    k27ac_e = subset(k27ac,rownames(k27ac)== k27acpeak);
    s = subset(snp_table,rownames(snp_table)==snp)
    variant <- as.vector(unlist(strsplit(rownames(s), "_")))
    s[s == "1|1"] <- paste(variant[4],variant[4],sep = "-")
    s[s == "0|1"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "1|0"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "0|0"] <- paste(variant[3],variant[3],sep = "-")
 
    mysamples <- colnames(k27ac_e) 
    s <- s[mysamples]
    rm(es)
    es=as.data.frame(t(rbind(k27ac_e,s)))
    colnames(es) <- c("gene","snp")
    es$gene<-as.numeric(as.character(es$gene))
    es$snp<-as.factor(es$snp)
    es$snp <- factor(es$snp, levels = c(paste(variant[3],variant[3],sep = "-"),paste(variant[3],variant[4],sep = "-"),paste(variant[4],variant[4],sep = "-")))
  
  
    pals <- rev(c("#9df4e9","#85d3ca","#78bbb5"))
    k27acQTL_plot <- ggplot(es,aes(x=snp,y=gene, color= snp)) +
      geom_violin(draw_quantiles = c(0.5), position = "dodge") +
      geom_boxplot(width=.1, outlier.colour=NA, position = "dodge") +
      geom_jitter(position=position_jitter(width=0.1), size=1.5) +
      xlab(paste(chr,snp_bp,sep=":")) + ylab(paste(k27acpeak," norm CQN",sep="")) +
      scale_color_manual(values = pals) +
      theme_minimal(base_size=14) +
      theme(legend.title=element_blank(),panel.grid.minor = element_blank(),legend.position="none",panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line( size=.05, color="#D9D9D9"),legend.text = element_text(size = 14),
            strip.background = element_blank(),
            strip.text.x = element_blank()) 
    ggsave(paste(gene_name, "/treg_actQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = k27acQTL_plot, width = 7, height = 6, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
  
  
    Major_actQTL <- es[es$snp == paste(variant[3],variant[3],sep = "-"),]
    Major_actQTL <- head(Major_actQTL[order(-Major_actQTL$gene),],5)
    Minor_actQTL <- es[es$snp == paste(variant[4],variant[4],sep = "-"),]
    Minor_actQTL <- head(Minor_actQTL[order(Minor_actQTL$gene),],5)
    heteroz_actQTL <- es[es$snp == paste(variant[3],variant[4],sep = "-"),]
    heteroz_actQTL <- heteroz_actQTL[order(-heteroz_actQTL$gene),]
    heteroz_actQTL_mean <- heteroz_actQTL[(ceiling(length(heteroz_actQTL$gene)/2)-2):(ceiling(length(heteroz_actQTL$gene)/2)+2),]
  
    print(Major_actQTL)
    print(heteroz_actQTL_mean)
    print(Minor_actQTL)
  }   
  ## Load K4me3 file
  if(k4me3peak != "chrNA"){
    k4me3_file_name <- c("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/Variance_Decomposition/ChM_K4me3_Peak_featureCounts_NoBatchCor_CQN_33PCsRegressed.txt")
    k4me3 = read.table(k4me3_file_name,sep="\t", header=T,row.names=4)
    k4me3 <- k4me3[,-(1:3)] 
    k4me3 <- k4me3[,-(1:2)] 
  
  ## plot your K4me3 peak and SNP of interest Violing plot
    k4me3_e = subset(k4me3,rownames(k4me3)== k4me3peak);
    s = subset(snp_table,rownames(snp_table)==snp)
    variant <- as.vector(unlist(strsplit(rownames(s), "_")))
    s[s == "1|1"] <- paste(variant[4],variant[4],sep = "-")
    s[s == "0|1"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "1|0"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "0|0"] <- paste(variant[3],variant[3],sep = "-") 
  
    mysamples <- colnames(k4me3_e) 
    s <- s[mysamples]
    es=as.data.frame(t(rbind(k4me3_e,s)))
    colnames(es) <- c("gene","snp")
    es$gene<-as.numeric(as.character(es$gene))
    es$snp<-as.factor(es$snp)
    es$snp <- factor(es$snp, levels = c(paste(variant[3],variant[3],sep = "-"),paste(variant[3],variant[4],sep = "-"),paste(variant[4],variant[4],sep = "-")))
  
  
    pals <- rev(c("#f4adab", "#ef877f", "#f47362"))
    k4me3QTL_plot <- ggplot(es,aes(x=snp,y=gene, color= snp)) +
      geom_violin(draw_quantiles = c(0.5), position = "dodge") +
      geom_boxplot(width=.1, outlier.colour=NA, position = "dodge") +
      geom_jitter(position=position_jitter(width=0.1), size=1.5) +
      xlab(paste(chr,snp_bp,sep=":")) + ylab(paste(k4me3peak," norm CQN",sep="")) +
      scale_color_manual(values = pals) +
      theme_minimal(base_size=14) +
      theme(legend.title=element_blank(),panel.grid.minor = element_blank(),legend.position="none",panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line( size=.05, color="#D9D9D9"),legend.text = element_text(size = 14),
            strip.background = element_blank(),
            strip.text.x = element_blank()) +
      ylim(-1.3,1.3)
    ggsave(paste(gene_name, "/treg_K4me3QTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = k4me3QTL_plot, width = 7, height = 6, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
  
  
    Major_k4me3QTL <- es[es$snp == paste(variant[3],variant[3],sep = "-"),]
    Major_k4me3QTL <- head(Major_k4me3QTL[order(-Major_k4me3QTL$gene),],5)
    Minor_k4me3QTL <- es[es$snp == paste(variant[4],variant[4],sep = "-"),]
    Minor_k4me3QTL <- head(Minor_k4me3QTL[order(Minor_k4me3QTL$gene),],5)
    heteroz_k4me3QTL <- es[es$snp == paste(variant[3],variant[4],sep = "-"),]
    heteroz_k4me3QTL <- heteroz_k4me3QTL[order(-heteroz_k4me3QTL$gene),]
    heteroz_k4me3QTL_mean <- heteroz_k4me3QTL[(ceiling(length(heteroz_k4me3QTL$gene)/2)-2):(ceiling(length(heteroz_k4me3QTL$gene)/2)+2),]
  
    print(Major_k4me3QTL)
    print(heteroz_k4me3QTL_mean)
    print(Minor_k4me3QTL)
  }
  
  ## Load ATAC file
  if(atacpeak != "chrNA"){
    atac_file_name <- c("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Integration/Variance_Decomposition/ATAC_Counts_featureCounts_CQN_30PCsRegressed.txt")
    atac = read.table(atac_file_name,sep="\t", header=T,row.names=4)
    atac <- atac[,-(1:3)] 
    atac <- atac[,-(1:2)] 
  
  ## plot your ATAC peak and SNP of interest Violing plot
    atac_e = subset(atac,rownames(atac)== atacpeak);
    s = subset(snp_table,rownames(snp_table)==snp)
    variant <- as.vector(unlist(strsplit(rownames(s), "_")))
    s[s == "1|1"] <- paste(variant[4],variant[4],sep = "-")
    s[s == "0|1"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "1|0"] <- paste(variant[3],variant[4],sep = "-")
    s[s == "0|0"] <- paste(variant[3],variant[3],sep = "-")  
  
    mysamples <- colnames(atac_e) 
    s <- s[mysamples]
    es=as.data.frame(t(rbind(atac_e,s)))
    colnames(es) <- c("gene","snp")
    es$gene<-as.numeric(as.character(es$gene))
    es$snp<-as.factor(es$snp)
    es$snp <- factor(es$snp, levels = c(paste(variant[3],variant[3],sep = "-"),paste(variant[3],variant[4],sep = "-"),paste(variant[4],variant[4],sep = "-")))
  
  
    pals <- rev(c("#efd6c5", "#b78b69", "#bb804f"))
    atacQTL_plot <- ggplot(es,aes(x=snp,y=gene, color= snp)) +
      geom_violin(draw_quantiles = c(0.5), position = "dodge") +
      geom_boxplot(width=.1, outlier.colour=NA, position = "dodge") +
      geom_jitter(position=position_jitter(width=0.1), size=1.5) +
      xlab(paste(chr,snp_bp,sep=":")) + ylab(paste(atacpeak," norm CQN",sep="")) +
      scale_color_manual(values = pals) +
      theme_minimal(base_size=14) +
      theme(legend.title=element_blank(),panel.grid.minor = element_blank(),legend.position="none",panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line( size=.05, color="#D9D9D9"),legend.text = element_text(size = 14),
            strip.background = element_blank(),
            strip.text.x = element_blank())
    ggsave(paste(gene_name, "/treg_ATACQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = atacQTL_plot, width = 7, height = 6, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
  
    Major_atacQTL <- es[es$snp == paste(variant[3],variant[3],sep = "-"),]
    Major_atacQTL <- head(Major_atacQTL[order(-Major_atacQTL$gene),],5)
    Minor_atacQTL <- es[es$snp == paste(variant[4],variant[4],sep = "-"),]
    Minor_atacQTL <- head(Minor_atacQTL[order(Minor_atacQTL$gene),],5)
    heteroz_atacQTL <- es[es$snp == paste(variant[3],variant[4],sep = "-"),]
    heteroz_atacQTL <- heteroz_atacQTL[order(-heteroz_atacQTL$gene),]
    heteroz_atacQTL_mean <- heteroz_atacQTL[(ceiling(length(heteroz_atacQTL$gene)/2)-2):(ceiling(length(heteroz_atacQTL$gene)/2)+2),]
  
    print(Major_atacQTL)
    print(heteroz_atacQTL_mean)
    print(Minor_atacQTL)
  }
}

##combined violin plots
if(gene != "Peak" & atacpeak != "chrNA" & k4me3peak != "chrNA" & k27acpeak != "chrNA"){
allviolin <- arrangeGrob(arrangeGrob(eQTL_plot), arrangeGrob(atacQTL_plot), arrangeGrob(k27acQTL_plot),arrangeGrob(k4me3QTL_plot),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = allviolin, width = 7, height = 6*as.numeric(length(ls(pattern= "QTL_plot", all.names = TRUE))), units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}
if(gene != "Peak" & atacpeak == "chrNA" & k4me3peak != "chrNA" & k27acpeak != "chrNA"){
allviolin <- arrangeGrob(arrangeGrob(eQTL_plot), arrangeGrob(k27acQTL_plot),arrangeGrob(k4me3QTL_plot),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = allviolin, width = 7, height = 6*as.numeric(length(ls(pattern= "QTL_plot", all.names = TRUE))), units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}
if(gene != "Peak" & atacpeak == "chrNA" & k4me3peak == "chrNA" & k27acpeak != "chrNA"){
allviolin <- arrangeGrob(arrangeGrob(eQTL_plot), arrangeGrob(k27acQTL_plot),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = allviolin, width = 7, height = 6*as.numeric(length(ls(pattern= "QTL_plot", all.names = TRUE))), units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene != "Peak" & atacpeak == "chrNA" & k4me3peak != "chrNA" & k27acpeak == "chrNA"){
allviolin <- arrangeGrob(arrangeGrob(eQTL_plot), arrangeGrob(k4me3QTL_plot),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = allviolin, width = 7, height = 6*as.numeric(length(ls(pattern= "QTL_plot", all.names = TRUE))), units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene != "Peak" & atacpeak != "chrNA" & k4me3peak == "chrNA" & k27acpeak == "chrNA"){
allviolin <- arrangeGrob(arrangeGrob(eQTL_plot), arrangeGrob(atacQTL_plot),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = allviolin, width = 7, height = 6*as.numeric(length(ls(pattern= "QTL_plot", all.names = TRUE))), units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene != "Peak" & atacpeak != "chrNA" & k4me3peak == "chrNA" & k27acpeak != "chrNA"){
allviolin <- arrangeGrob(arrangeGrob(eQTL_plot), arrangeGrob(atacQTL_plot), arrangeGrob(k27acQTL_plot),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = allviolin, width = 7, height = 6*as.numeric(length(ls(pattern= "QTL_plot", all.names = TRUE))), units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene == "Peak" & atacpeak != "chrNA" & k4me3peak != "chrNA" & k27acpeak != "chrNA"){
allviolin <- arrangeGrob(arrangeGrob(atacQTL_plot), arrangeGrob(k27acQTL_plot),arrangeGrob(k4me3QTL_plot),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = allviolin, width = 7, height = 6*3, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}
if(gene == "Peak" & atacpeak == "chrNA" & k4me3peak != "chrNA" & k27acpeak != "chrNA"){
allviolin <- arrangeGrob(arrangeGrob(k27acQTL_plot),arrangeGrob(k4me3QTL_plot),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = allviolin, width = 7, height = 6*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}
if(gene == "Peak" & atacpeak == "chrNA" & k4me3peak == "chrNA" & k27acpeak != "chrNA"){
allviolin <- arrangeGrob(arrangeGrob(k27acQTL_plot),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = allviolin, width = 7, height = 6, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene == "Peak" & atacpeak == "chrNA" & k4me3peak != "chrNA" & k27acpeak == "chrNA"){
allviolin <- arrangeGrob(arrangeGrob(k4me3QTL_plot),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = allviolin, width = 7, height = 6, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene == "Peak" & atacpeak != "chrNA" & k4me3peak == "chrNA" & k27acpeak == "chrNA"){
allviolin <- arrangeGrob(arrangeGrob(atacQTL_plot),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = allviolin, width = 7, height = 6, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene == "Peak" & atacpeak != "chrNA" & k4me3peak == "chrNA" & k27acpeak != "chrNA"){
allviolin <- arrangeGrob(arrangeGrob(atacQTL_plot), arrangeGrob(k27acQTL_plot),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_","_Violin.pdf", sep=""), plot = allviolin, width = 7, height = 6*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

#Plot K27ac-seq wide region coverage plots with genotype stratification
# make annotation for peaks

extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                         qValue = "numeric")
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
# define regions of interest 
chr = chr
from = as.numeric(snp_bp) - as.numeric(window) 
to = as.numeric(snp_bp) + as.numeric(window)

#Load K27ac peaks
if(k27acpeak != "chrNA"){
  k27acpeaks <- import("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/Merged_Day0_treg_ChM_K27ac_sorted-broad_peaks_qval_0.001_FE_2.broadPeak", format = "BED",extraCols = extraCols_broadPeak)
  k27ac_selected_peaks <- k27acpeaks[k27acpeaks@seqnames == chr & start(k27acpeaks@ranges) > from & end(k27acpeaks@ranges) < to, ]
  k27ac_peak_list <- split(k27ac_selected_peaks, as.factor(k27ac_selected_peaks))

#get SNP info
  selected_ids_k27ac <- c(row.names(Major_actQTL))
  s = subset(snp_table,rownames(snp_table)==snp)
  s = s[ , which(names(s) %in% selected_ids_k27ac)]
  s[s == "0|0"] <- c("1_Hom_Maj")
  s[s == "0|1"] <- c("2_Het")
  s[s == "1|0"] <- c("2_Het")
  s[s == "1|1"] <- c("3_Hom_Min")

#Create sample data table for wiggleplotr
  ids <- as.vector(s)
  colnames(ids) <- NULL
  row.names(ids) <- NULL 
  sample_data = dplyr::data_frame(sample_id = names(s), condition = as.character(ids),scaling_factor = 1)
  sample_data = sample_data %>%
    dplyr::mutate(bigWig = paste0("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/BigWig/",sample_id, "_treg_ChM_K27ac_Merged.trimmed.ReAligned_GRCh38.sorted2_uniq_nodup_Autosomes_SPMR.bw",sep=""))
  as.data.frame(sample_data)
  track_data = dplyr::mutate(sample_data, track_id = condition, colour_group = condition)
  track_data2 = dplyr::mutate(sample_data, track_id = "H3K27ac", colour_group = condition)

### and then to make the plot
#line
  q1 <- plotCoverage2(track_data = track_data2,
               exons = k27ac_peak_list, 
               cdss = k27ac_peak_list,
               fill_palette = c("#78bbb5"),
               mean_only = TRUE, region_coords = c(from,to),
               connect_exons = FALSE,  rescale_introns = FALSE,
               transcript_label = FALSE, 
               coverage_type="both", return_subplots_list = TRUE,
               plot_fraction = 0.005)
#  ggsave(paste(gene_name,"/treg_actQTL_",gene_name, "_", gene, "_", snp, "_",window, "_coverage.pdf", sep=""), plot = q1$coverage_plot, width = 15, height = 2, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)

###Smoothed coverage plot

coverage <- as.data.frame(q1$bins)
coverage$coverage <- q1$coverage
coverage$colour_group <- q1$colour_group
coverage$sample_id <- q1$sample_id
names(coverage) <- c("bins","coverage", "colour_group", "sample_id")
coverage2 <- as.data.frame(tapply(coverage$coverage, cut(coverage$bins, seq(min(coverage$bins), max(coverage$bins), by=5000)), max))
names(coverage2) <- c("coverage")
coverage3 <- as.data.frame(tapply(coverage$bins, cut(coverage$bins, seq(min(coverage$bins), max(coverage$bins), by=5000)), mean))
names(coverage3) <- c("new_bin")
coverage4 <- as.data.frame(tapply(coverage$colour_group, cut(coverage$bins, seq(min(coverage$bins), max(coverage$bins), by=5000)), max))
names(coverage4) <- c("colour_group")
coverage5 <- as.data.frame(tapply(coverage$sample_id, cut(coverage$bins, seq(min(coverage$bins), max(coverage$bins), by=5000)), unique))
names(coverage5) <- c("sample_id")
coverage6 <- merge(coverage3,coverage2, by="row.names")
coverage7 <- merge(coverage6,coverage4, by.x = "Row.names", by.y="row.names")
coverage8_k27ac <- merge(coverage7,coverage5, by.x = "Row.names", by.y="row.names")
names(coverage8_k27ac) <- c("bin","new_bin","coverage", "colour_group", "sample_id")


spline_int_k27ac <- as.data.frame(spline(coverage8_k27ac$new_bin, coverage8_k27ac$coverage))
spline_int_k27ac$y <- ifelse(spline_int_k27ac$y < 0, 0, spline_int_k27ac$y)

q1 = ggplot(coverage8_k27ac, aes(coverage8_k27ac$new_bin, coverage8_k27ac$coverage)) + 
  geom_blank() +
  theme_light()+
  geom_line(data = spline_int_k27ac, aes(x = x, y = y), colour=c("#78bbb5"))+ 
  geom_area(data = spline_int_k27ac, aes(x = x, y = y),fill=c("#78bbb5")) +
  ylab("SPMR")+xlab(paste(chr, "bp"))+
  scale_y_continuous(labels = function(x) round(as.numeric(x), digits=0))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "#d9d9d9"),
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())+
                     ylim(0, max(as.numeric(spline_int_k27ac$y))+0.2)
#ggsave(paste(gene_name,"/treg_actQTL_",gene_name, "_", gene, "_", snp, "_",window, "_coverage.pdf", sep=""), plot = q1, width = 15, height = 2, units = "cm",
#       dpi = 300,device = "pdf", useDingbats=FALSE)
}
#Load K4me3 peaks
if(k4me3peak != "chrNA"){
  k4me3peaks <- import("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/Merged_Day0_treg_ChM_K4me3_sorted-sharp_peaks_qval_0.001_FE_2.narrowPeak", format = "BED",extraCols = extraCols_narrowPeak)
  k4me3_selected_peaks <- k4me3peaks[k4me3peaks@seqnames == chr & start(k4me3peaks@ranges) > from & end(k4me3peaks@ranges) < to, ]
  k4me3_peak_list <- split(k4me3_selected_peaks, as.factor(k4me3_selected_peaks))
#get SNP info
  selected_ids_k4me3 <- c(row.names(Major_k4me3QTL))
  s = subset(snp_table,rownames(snp_table)==snp)
  s = s[ , which(names(s) %in% selected_ids_k4me3)]
  s[s == "0|0"] <- c("1_Hom_Maj")
  s[s == "0|1"] <- c("2_Het")
  s[s == "1|0"] <- c("2_Het")
  s[s == "1|1"] <- c("3_Hom_Min")

#Create data table for wiggleplotr
  ids <- as.vector(s)
  colnames(ids) <- NULL
  row.names(ids) <- NULL 
  sample_data = dplyr::data_frame(sample_id = names(s), condition = as.character(ids),scaling_factor = 1)
  sample_data = sample_data %>%
    dplyr::mutate(bigWig = paste0("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/BigWig/",sample_id, "_treg_ChM_K4me3_Merged.trimmed.ReAligned_GRCh38.sorted2_uniq_nodup_Autosomes-sharp_SPMR.bw",sep=""))
  as.data.frame(sample_data)
  track_data = dplyr::mutate(sample_data, track_id = condition, colour_group = condition)
  track_data2 = dplyr::mutate(sample_data, track_id = "H3K4me3", colour_group = condition)

### and then to make the plot
#line
  q2 <- plotCoverage2(track_data = track_data2,
                     exons = k4me3_peak_list, 
                     cdss = k4me3_peak_list,
                     fill_palette = c("#f9aba5"),
                     mean_only = TRUE, region_coords = c(from,to),
                     connect_exons = FALSE,  rescale_introns = FALSE,
                     transcript_label = FALSE, 
                     coverage_type="both", return_subplots_list = TRUE,
                     plot_fraction = 0.005)
#  ggsave(paste(gene_name, "/treg_K4me3QTL_",gene_name, "_", gene, "_", snp, "_",window, "_coverage.pdf", sep =""), plot = q2$coverage_plot, width = 15, height = 2, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)

###Smoothed coverage plot

coverage <- as.data.frame(q2$bins)
coverage$coverage <- q2$coverage
coverage$colour_group <- q2$colour_group
coverage$sample_id <- q2$sample_id
names(coverage) <- c("bins","coverage", "colour_group", "sample_id")
coverage2 <- as.data.frame(tapply(coverage$coverage, cut(coverage$bins, seq(min(coverage$bins), max(coverage$bins), by=5000)), max))
names(coverage2) <- c("coverage")
coverage3 <- as.data.frame(tapply(coverage$bins, cut(coverage$bins, seq(min(coverage$bins), max(coverage$bins), by=5000)), mean))
names(coverage3) <- c("new_bin")
coverage4 <- as.data.frame(tapply(coverage$colour_group, cut(coverage$bins, seq(min(coverage$bins), max(coverage$bins), by=5000)), max))
names(coverage4) <- c("colour_group")
coverage5 <- as.data.frame(tapply(coverage$sample_id, cut(coverage$bins, seq(min(coverage$bins), max(coverage$bins), by=5000)), unique))
names(coverage5) <- c("sample_id")
coverage6 <- merge(coverage3,coverage2, by="row.names")
coverage7 <- merge(coverage6,coverage4, by.x = "Row.names", by.y="row.names")
coverage8_k4me3 <- merge(coverage7,coverage5, by.x = "Row.names", by.y="row.names")
names(coverage8_k4me3) <- c("bin","new_bin","coverage", "colour_group", "sample_id")


spline_int_k4me3 <- as.data.frame(spline(coverage8_k4me3$new_bin, coverage8_k4me3$coverage))
spline_int_k4me3$y <- ifelse(spline_int_k4me3$y < 0, 0, spline_int_k4me3$y)

q2 = ggplot(coverage8_k4me3, aes(coverage8_k4me3$new_bin, coverage8_k4me3$coverage)) +  
  geom_blank() +   
  theme_light()+
  geom_line(data = spline_int_k4me3, aes(x = x, y = y), colour=c("#f9aba5"))+  
  geom_area(data = spline_int_k4me3, aes(x = x, y = y),fill=c("#f9aba5")) +
  ylab("SPMR")+xlab(paste(chr, "bp"))+
  scale_y_continuous(labels = function(x) round(as.numeric(x), digits=0))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),  
                     axis.line = element_line(colour = "#d9d9d9"),
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())+
                     ylim(0, max(as.numeric(spline_int_k4me3$y))+0.2)
#ggsave(paste(gene_name,"/treg_K4me3QTL_",gene_name, "_", gene, "_", snp, "_",window, "_coverage.pdf", sep=""), plot = q2, width = 15, height = 2, units = "cm",
#       dpi = 300,device = "pdf", useDingbats=FALSE)
}
#Load ATAC peaks
if(atacpeak != "chrNA"){
  atacpeaks <- import("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ATACseq/QTLtools/featureCounts/2MRandomReads_Autosomes_sorted_10Reads58Samples.narrowPeak", format = "BED",extraCols = extraCols_narrowPeak)
  atac_selected_peaks <- atacpeaks[atacpeaks@seqnames == chr & start(atacpeaks@ranges) > from & end(atacpeaks@ranges) < to, ]
  atac_peak_list <- split(atac_selected_peaks, as.factor(atac_selected_peaks))
#get SNP info
  selected_ids_atac <- c(row.names(Major_atacQTL))

  s = subset(snp_table,rownames(snp_table)==snp)
  s = s[ , which(names(s) %in% selected_ids_atac)]
  s[s == "0|0"] <- c("1_Hom_Maj")
  s[s == "1|0"] <- c("2_Het")
  s[s == "0|1"] <- c("2_Het")
  s[s == "1|1"] <- c("3_Hom_Min")

  ids <- as.vector(s)
  colnames(ids) <- NULL
  row.names(ids) <- NULL
  sample_data = dplyr::data_frame(sample_id = names(s), condition = as.character(ids),scaling_factor = 1)
  sample_data = sample_data %>%
    dplyr::mutate(bigWig = paste0("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ATACseq/BigWigs/",sample_id, "_ATAC.bw",sep=""))
  as.data.frame(sample_data)
  track_data = dplyr::mutate(sample_data, track_id = condition, colour_group = condition)
  track_data2 = dplyr::mutate(sample_data, track_id = "ATAC", colour_group = condition)

### and then to make the plot
#line
  q6 <- plotCoverage2(track_data = track_data2,
                     exons = atac_peak_list,
                     cdss = atac_peak_list,
                     fill_palette = c("#efd6c5"),
                     mean_only = TRUE, region_coords = c(from,to),
                     connect_exons = FALSE,  rescale_introns = FALSE,
                     transcript_label = FALSE,
                     coverage_type="both", return_subplots_list = TRUE,
                     plot_fraction = 0.005)
#  ggsave(paste(gene_name, "/treg_atacQTL_",gene_name, "_", gene, "_", snp, "_",window, "_coverage.pdf", sep =""), plot = q6$coverage_plot, width = 15, height = 2, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)

###Smoothed coverage plot

coverage <- as.data.frame(q6$bins)
coverage$coverage <- q6$coverage
coverage$colour_group <- q6$colour_group
coverage$sample_id <- q6$sample_id
names(coverage) <- c("bins","coverage", "colour_group", "sample_id")
coverage2 <- as.data.frame(tapply(coverage$coverage, cut(coverage$bins, seq(min(coverage$bins), max(coverage$bins), by=5000)), max))
names(coverage2) <- c("coverage")
coverage3 <- as.data.frame(tapply(coverage$bins, cut(coverage$bins, seq(min(coverage$bins), max(coverage$bins), by=5000)), mean))
names(coverage3) <- c("new_bin")
coverage4 <- as.data.frame(tapply(coverage$colour_group, cut(coverage$bins, seq(min(coverage$bins), max(coverage$bins), by=5000)), max))
names(coverage4) <- c("colour_group")
coverage5 <- as.data.frame(tapply(coverage$sample_id, cut(coverage$bins, seq(min(coverage$bins), max(coverage$bins), by=5000)), unique))
names(coverage5) <- c("sample_id")
coverage6 <- merge(coverage3,coverage2, by="row.names")
coverage7 <- merge(coverage6,coverage4, by.x = "Row.names", by.y="row.names")
coverage8_atac <- merge(coverage7,coverage5, by.x = "Row.names", by.y="row.names")
names(coverage8_atac) <- c("bin","new_bin","coverage", "colour_group", "sample_id")


spline_int_atac <- as.data.frame(spline(coverage8_atac$new_bin, coverage8_atac$coverage))
spline_int_atac$y <- ifelse(spline_int_atac$y < 0, 0, spline_int_atac$y)

q6 = ggplot(coverage8_atac, aes(coverage8_atac$new_bin, coverage8_atac$coverage)) +
  geom_blank() +
  theme_light()+
  geom_line(data = spline_int_atac, aes(x = x, y = y), colour=c("#efd6c5"))+
  geom_area(data = spline_int_atac, aes(x = x, y = y),fill=c("#efd6c5")) +
  ylab("SPMR")+xlab(paste(chr, "bp"))+
  scale_y_continuous(labels = function(x) round(as.numeric(x), digits=0))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "#d9d9d9"),
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())+
                     ylim(0, max(as.numeric(spline_int_atac$y))+0.2)
#ggsave(paste(gene_name,"/treg_atacQTL_",gene_name, "_", gene, "_", snp, "_",window, "_coverage.pdf", sep=""), plot = q6, width = 15, height = 2, units = "cm",
#       dpi = 300,device = "pdf", useDingbats=FALSE)

}

##combined coverage plots
if(atacpeak != "chrNA" & k4me3peak != "chrNA" & k27acpeak != "chrNA"){
allcor <- arrangeGrob(arrangeGrob(q6), arrangeGrob(q1), arrangeGrob(q2),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_coverage.pdf", sep =""), plot = allcor, width = 15, height = 2*3, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(atacpeak == "chrNA" & k4me3peak != "chrNA" & k27acpeak != "chrNA"){
allcor <- arrangeGrob(arrangeGrob(q1), arrangeGrob(q2),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_coverage.pdf", sep =""), plot = allcor, width = 15, height = 2*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(atacpeak != "chrNA" & k4me3peak == "chrNA" & k27acpeak != "chrNA"){
allcor <- arrangeGrob(arrangeGrob(q6), arrangeGrob(q1),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_coverage.pdf", sep =""), plot = allcor, width = 15, height = 2*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(atacpeak != "chrNA" & k4me3peak != "chrNA" & k27acpeak == "chrNA"){
allcor <- arrangeGrob(arrangeGrob(q6), arrangeGrob(q2),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_coverage.pdf", sep =""), plot = allcor, width = 15, height = 2*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(atacpeak != "chrNA" & k4me3peak == "chrNA" & k27acpeak == "chrNA"){
allcor <- arrangeGrob(arrangeGrob(q6),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_coverage.pdf", sep =""), plot = allcor, width = 15, height = 2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(atacpeak == "chrNA" & k4me3peak != "chrNA" & k27acpeak == "chrNA"){
allcor <- arrangeGrob(arrangeGrob(q2),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_coverage.pdf", sep =""), plot = allcor, width = 15, height = 2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(atacpeak == "chrNA" & k4me3peak == "chrNA" & k27acpeak != "chrNA"){
allcor <- arrangeGrob(arrangeGrob(q1),ncol=1)
ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_coverage.pdf", sep =""), plot = allcor, width = 15, height = 2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

#Get GWAS SNP LD block and plot plot them in the zoomed region coordinates
#SNPs
LD <- read.table(paste(snp, ".ld", sep=""), sep="\t")
names(LD) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
Add_autoLD <- c(chr, snp_bp, snp, chr, snp_bp, snp, "1", "NA")
LD <- rbind(LD, Add_autoLD)
LD2 <- LD[ which(LD$bp1==snp_bp),]
LD3 <- LD[ which(LD$bp2==snp_bp),]
LD3 <- LD3[,c(1,5,6,4,2,3,7,8)]
names(LD3) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
LD4 <- rbind(LD2,LD3)
LD4 <- LD4[ which(LD4$r2 >= 0.8),]
if(k27acpeak != "chrNA"){
LD4$color <- ifelse((LD4$bp2 >= as.numeric(k27acpeak_region[2])) & (LD4$bp2 <= as.numeric(k27acpeak_region[3])), "#f78484", "#a8a8a8")
}
if(k27acpeak == "chrNA"){
LD4$color <- c(rep("#a8a8a8",length(LD4$r2)))
}
LD4 <- LD4[order(LD4$bp2),]
from_zoomed = min(as.numeric(LD4$bp2))-55000
to_zoomed = max(as.numeric(LD4$bp2))+55000

if(k27acpeak != "chrNA"){
  SNP_k27ac <-  ggplot(LD4, aes(x=as.numeric(LD4$bp2), y=rep(1,length(LD4$chr_snp))))+
                      annotate("rect", xmin=as.numeric(k27acpeak_region[2]), xmax=as.numeric(k27acpeak_region[3]), ymin=0, ymax=Inf, alpha=0.5, fill="#f4d7c3")+
                      geom_bar(stat="identity", color = as.character(LD4$color))+
                      coord_cartesian(xlim = c(from_zoomed, to_zoomed), ylim =  c(0,1))+
                      xlab(NULL)+ ylab(NULL)+
                      theme_classic()+
                      theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line = element_line(colour = "#d9d9d9"),
                      legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
  ggsave(paste(gene_name, "/treg_",gene_name, "_", gene, "_", snp, "_","_Zoomed_SNPs_K27ac.pdf", sep =""), plot = SNP_k27ac , width = 15, height = 2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)  
}


if(k4me3peak != "chrNA"){
  LD4$color <- ifelse((LD4$bp2 >= as.numeric(k4me3peak_region[2])) & (LD4$bp2 <= as.numeric(k4me3peak_region[3])), "#f78484", "#a8a8a8")
  SNP_k4me3 <-  ggplot(LD4, aes(x=as.numeric(LD4$bp2), y=rep(1,length(LD4$chr_snp))))+
                      annotate("rect", xmin=as.numeric(k4me3peak_region[2]), xmax=as.numeric(k4me3peak_region[3]), ymin=0, ymax=Inf, alpha=0.5, fill="#f4d7c3")+
                      geom_bar(stat="identity", color = as.character(LD4$color))+
                      coord_cartesian(xlim = c(from_zoomed, to_zoomed), ylim =  c(0,1))+
                      xlab(NULL)+ ylab(NULL)+
                      theme_classic()+
                      theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line = element_line(colour = "#d9d9d9"),
                      legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())

  ggsave(paste(gene_name, "/treg_",gene_name, "_", gene, "_", snp, "_","_Zoomed_SNPs_K4me3.pdf", sep =""), plot = SNP_k4me3 , width = 15, height = 2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)  

}

if(atacpeak != "chrNA"){
  LD4$color <- ifelse((LD4$bp2 >= as.numeric(atacpeak_region[2])) & (LD4$bp2 <= as.numeric(atacpeak_region[3])), "#f78484", "#a8a8a8")
  SNP_ATAC <-  ggplot(LD4, aes(x=as.numeric(LD4$bp2), y=rep(1,length(LD4$chr_snp))))+
                      annotate("rect", xmin=as.numeric(atacpeak_region[2]), xmax=as.numeric(atacpeak_region[3]), ymin=0, ymax=Inf, alpha=0.5, fill="#f4d7c3")+
                      geom_bar(stat="identity", color = as.character(LD4$color))+
                      coord_cartesian(xlim = c(from_zoomed, to_zoomed), ylim =  c(0,1))+
                      xlab(NULL)+ ylab(NULL)+
                      theme_classic()+
                      theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line = element_line(colour = "#d9d9d9"),
                      legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
  ggsave(paste(gene_name, "/treg_",gene_name, "_", gene, "_", snp, "_","_Zoomed_SNPs_ATAC.pdf", sep =""), plot = SNP_ATAC , width = 15, height = 2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)  

}

# Plot genes in region
# set gene track
# With get names
ucscGenes <- UcscTrack(genome="hg38", table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=chr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)
z <- ranges(ucscGenes)

mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z
df2 <- data.frame(chr=seqnames(ucscGenes2),
                  starts=start(ucscGenes2),
                  ends=end(ucscGenes2),
                  names=symbol(ucscGenes2),
                  strands=strand(ucscGenes2),
                  symbol = symbol(ucscGenes2), transcript = transcript(ucscGenes2),
                  feature=feature(ucscGenes2))
md2<- df2 %>% group_by(transcript) %>% filter(starts==min(starts))
md3<- df2 %>% group_by(transcript) %>% filter(ends==max(ends))
md4 <- md2[ , -which(names(md2) %in% c("ends"))]
md4$ends <- md3$ends
md4$length <- md4$ends - md4$starts
md5 <- md4 %>% 
  group_by(names) %>%
  filter(length == max(length)) %>%
  arrange(chr,starts,names,strands,symbol,transcript,feature,ends)
df3 <- df2[df2$transcript %in% md5$transcript ,]
df3$feature <- ifelse(df3$symbol == gene_name, 'red', 'black')
names(df3) <- c("chr", "start","end","name", "strand", "symbol","transcript","feature")
ucscGenes3 <- AnnotationTrack(chromosome = chr, start = df3$start,end = df3$end, strand = df3$strand,
                              feature = as.vector(df3$feature), id = as.vector(df3$feature),  exon = as.vector(df3$feature), 
                              group = as.character(df3$symbol), genome = "hg38", 
                              name = "UCSC genes", 
                              showId=TRUE, geneSymbol=TRUE,just.group="above",
                              fill =as.vector(df3$feature),col=NULL,filled.contour=as.vector(df3$feature), 
                              #background.panel = "#FFFEDB", background.title = "darkorange", 
                              red="#f78484", black="#a8a8a8",
                              fontsize = 18, cex.group=1)
z <- ranges(ucscGenes3)
mcols(z)$exon <- df3$feature
mcols(z)$gene <- df3$transcript
mcols(z)$transcript <- df3$transcript
ucscGenes4 <- ucscGenes3 
ranges(ucscGenes4) <- z

#Gather track lists
track.list <- list(ucscGenes4)

#Highlight regions
ht1 <- HighlightTrack(trackList = track.list, start = from_zoomed, 
                      width = to_zoomed - from_zoomed, 
                      chromosome = chr, col=c("#efefef"), fill=c("#f9f9f9"))


pdf(file = paste(gene_name, "/treg_",gene_name, "_", gene, "_", snp, "_", window, "_Region.pdf", sep=""), width = 6.40, height = 3)
plotTracks(ht1,chromosome = chr,from = from, to = to, 
           sizes=c(1), fontsize = 10)
dev.off()

##############################################################################
################################################################
##############################################

###ZOOMED PLOTS

#Load K27ac peaks
if(k27acpeak != "chrNA"){
  k27acpeaks <- import("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/Merged_Day0_treg_ChM_K27ac_sorted-broad_peaks_qval_0.001_FE_2.broadPeak", format = "BED",extraCols = extraCols_broadPeak)
  k27ac_selected_peaks <- k27acpeaks[k27acpeaks@seqnames == chr & start(k27acpeaks@ranges) > from_zoomed & end(k27acpeaks@ranges) < to_zoomed, ]
  k27ac_peak_list <- split(k27ac_selected_peaks, as.factor(k27ac_selected_peaks))

#get SNP info
  selected_ids_k27ac <- c(row.names(Major_actQTL),row.names(Minor_actQTL))
  s = subset(snp_table,rownames(snp_table)==snp)
  s = s[ , which(names(s) %in% selected_ids_k27ac)]
  s[s == "0|0"] <- c("1_Hom_Maj")
  s[s == "1|0"] <- c("2_Het")
  s[s == "0|1"] <- c("2_Het")
  s[s == "1|1"] <- c("3_Hom_Min")

  ids <- as.vector(s)
  colnames(ids) <- NULL
  row.names(ids) <- NULL 
  sample_data = dplyr::data_frame(sample_id = names(s), condition = as.character(ids),scaling_factor = 1)
  sample_data = sample_data %>%
    dplyr::mutate(bigWig = paste0("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/BigWig/",sample_id, "_treg_ChM_K27ac_Merged.trimmed.ReAligned_GRCh38.sorted2_uniq_nodup_Autosomes_SPMR.bw",sep=""))
  as.data.frame(sample_data)
  track_data = dplyr::mutate(sample_data, track_id = condition, colour_group = condition)
  track_data2 = dplyr::mutate(sample_data, track_id = "H3K27ac", colour_group = condition)

### and then to make the plot
#line
  q3 <- plotCoverage2(track_data = track_data2,
                     exons = k27ac_peak_list, 
                     cdss = k27ac_peak_list,
                     #fill_palette = c("#C80D04","#F73930","#F96C66"),
                     fill_palette = rev(c("#9df4e9","#78bbb5")),
                     mean_only = TRUE, region_coords = c(from_zoomed,to_zoomed),
                     connect_exons = FALSE,  rescale_introns = FALSE,
                     transcript_label = FALSE, 
                     coverage_type="both", return_subplots_list = TRUE,
                     plot_fraction = 0.005)

  k27acpeaks <- read.table("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/Merged_Day0_treg_ChM_K27ac_sorted-broad_peaks_qval_0.001_FE_2.broadPeak", sep="\t")
  k27acpeaks <- k27acpeaks[ which(k27acpeaks$V1==chr),]
  k27acpeaks <- k27acpeaks[ which(k27acpeaks$V1==chr),] 
  k27acpeaks <- k27acpeaks[ which(k27acpeaks$V3 >= from_zoomed),] 
  k27acpeaks <- k27acpeaks[ which(k27acpeaks$V2 <= to_zoomed),] 

  t2.rect1 <- data.frame (name= paste(k27acpeaks$V1, k27acpeaks$V2, k27acpeaks$V3, k27acpeaks$V4, sep="_"),xmin=k27acpeaks$V2, xmax=k27acpeaks$V3, ymin=rep(0, length(k27acpeaks$V1)), ymax=rep(Inf, length(k27acpeaks$V1)))
  t2.rect1$color <- ifelse(t2.rect1$name == k27acpeak, "#f9f973", "#cccccc")

  q3_2 <- q3$coverage_plot + 
          geom_rect(data=t2.rect1, aes(xmin=as.numeric(xmin), xmax=as.numeric(xmax), 
                               ymin=as.numeric(ymin), ymax=as.numeric(ymax)), 
                               fill = t2.rect1$color, alpha=0.3, inherit.aes = FALSE)
#  ggsave(file = paste(gene_name, "/treg_actQTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_coverage.pdf", sep=""), plot = q3_2 , width = 15, height = 4, units = "cm",
#           dpi = 300,device = "pdf", useDingbats=FALSE)

###Smoothed zoomed coverage plots
coverage <- as.data.frame(q3$bins)
coverage$coverage <- q3$coverage
coverage$colour_group <- q3$colour_group
coverage$sample_id <- q3$sample_id
names(coverage) <- c("bins","coverage", "colour_group", "sample_id")
coverage <- coverage[complete.cases(coverage), ]
coverage_Hom_Maj <- coverage[ which(coverage$colour_group == c("1_Hom_Maj")), ] 
coverage_Hom_Min <- coverage[ which(coverage$colour_group == "3_Hom_Min"), ] 

coverage2 <- as.data.frame(tapply(coverage_Hom_Maj$coverage, cut(coverage_Hom_Maj$bins, seq(min(coverage_Hom_Maj$bins), max(coverage_Hom_Maj$bins), by=1000)), max))
names(coverage2) <- c("coverage")
coverage3 <- as.data.frame(tapply(coverage_Hom_Maj$bins, cut(coverage_Hom_Maj$bins, seq(min(coverage_Hom_Maj$bins), max(coverage_Hom_Maj$bins), by=1000)), mean))
names(coverage3) <- c("new_bin")
coverage4 <- as.data.frame(tapply(coverage_Hom_Maj$colour_group, cut(coverage_Hom_Maj$bins, seq(min(coverage_Hom_Maj$bins), max(coverage_Hom_Maj$bins), by=1000)), unique))
names(coverage4) <- c("colour_group")
coverage5 <- as.data.frame(tapply(coverage_Hom_Maj$sample_id, cut(coverage_Hom_Maj$bins, seq(min(coverage_Hom_Maj$bins), max(coverage_Hom_Maj$bins), by=1000)), unique))
names(coverage5) <- c("sample_id")
coverage6 <- merge(coverage3,coverage2, by="row.names")
coverage7 <- merge(coverage6,coverage4, by.x = "Row.names", by.y="row.names")
coverage8_Hom_Maj_k27ac <- merge(coverage7,coverage5, by.x = "Row.names", by.y="row.names")
names(coverage8_Hom_Maj_k27ac) <- c("bin","new_bin","coverage", "colour_group", "sample_id")

coverage2 <- as.data.frame(tapply(coverage_Hom_Min$coverage, cut(coverage_Hom_Min$bins, seq(min(coverage_Hom_Min$bins), max(coverage_Hom_Min$bins), by=1000)), max))
names(coverage2) <- c("coverage")
coverage3 <- as.data.frame(tapply(coverage_Hom_Min$bins, cut(coverage_Hom_Min$bins, seq(min(coverage_Hom_Min$bins), max(coverage_Hom_Min$bins), by=1000)), mean))
names(coverage3) <- c("new_bin")
coverage4 <- as.data.frame(tapply(coverage_Hom_Min$colour_group, cut(coverage_Hom_Min$bins, seq(min(coverage_Hom_Min$bins), max(coverage_Hom_Min$bins), by=1000)), unique))
names(coverage4) <- c("colour_group")
coverage5 <- as.data.frame(tapply(coverage_Hom_Min$sample_id, cut(coverage_Hom_Min$bins, seq(min(coverage_Hom_Min$bins), max(coverage_Hom_Min$bins), by=1000)), unique))
names(coverage5) <- c("sample_id")
coverage6 <- merge(coverage3,coverage2, by="row.names")
coverage7 <- merge(coverage6,coverage4, by.x = "Row.names", by.y="row.names")
coverage8_Hom_Min_k27ac <- merge(coverage7,coverage5, by.x = "Row.names", by.y="row.names")
names(coverage8_Hom_Min_k27ac) <- c("bin","new_bin","coverage", "colour_group", "sample_id")
coverage8_k27ac_zoomed <- rbind(coverage8_Hom_Maj_k27ac,coverage8_Hom_Min_k27ac)



spline_int_1_k27ac <- as.data.frame(spline(coverage8_Hom_Maj_k27ac$new_bin, coverage8_Hom_Maj_k27ac$coverage))
spline_int_1_k27ac$y <- ifelse(spline_int_1_k27ac$y < 0, 0, spline_int_1_k27ac$y)
spline_int_2_k27ac <- as.data.frame(spline(coverage8_Hom_Min_k27ac$new_bin, coverage8_Hom_Min_k27ac$coverage))
spline_int_2_k27ac$y <- ifelse(spline_int_2_k27ac$y < 0, 0, spline_int_2_k27ac$y)

if(slope < 0){
  k27acpeaks <- read.table("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/Merged_Day0_treg_ChM_K27ac_sorted-broad_peaks_qval_0.001_FE_2.broadPeak", sep="\t")
  k27acpeaks <- k27acpeaks[ which(k27acpeaks$V1==chr),]
  k27acpeaks <- k27acpeaks[ which(k27acpeaks$V1==chr),]
  k27acpeaks <- k27acpeaks[ which(k27acpeaks$V3 >= from_zoomed),]
  k27acpeaks <- k27acpeaks[ which(k27acpeaks$V2 <= to_zoomed),]

  t2.rect1_k27ac <- data.frame (name= paste(k27acpeaks$V1, k27acpeaks$V2, k27acpeaks$V3, k27acpeaks$V4, sep="_"),xmin=k27acpeaks$V2, xmax=k27acpeaks$V3, ymin=rep(0, length(k27acpeaks$V1)), ymax=rep(Inf, length(k27acpeaks$V1)))
  t2.rect1_k27ac$color <- ifelse(t2.rect1_k27ac$name == k27acpeak, "#F4D7C3", "#d9d9d9" )

  q3_2 = ggplot(coverage8_k27ac_zoomed, aes(coverage8_k27ac_zoomed$new_bin, coverage8_k27ac_zoomed$coverage)) + 
    geom_blank() +
    theme_light()+
    geom_rect(data=t2.rect1_k27ac, aes(xmin=as.numeric(xmin), xmax=as.numeric(xmax),
                                 ymin=as.numeric(ymin), ymax=as.numeric(ymax)),
              fill = t2.rect1_k27ac$color, alpha=0.5, inherit.aes = FALSE)+
    geom_line(data = spline_int_1_k27ac, aes(x = x, y = y), colour=c("#78bbb5"))+ 
    geom_area(data = spline_int_1_k27ac, aes(x = x, y = y),fill=c("#78bbb5"), alpha=0.8) +
    geom_line(data = spline_int_2_k27ac, aes(x = x, y = y), colour=c("#9df4e9"))+ 
    geom_area(data = spline_int_2_k27ac, aes(x = x, y = y),fill=c("#9df4e9"), alpha=0.8) +  
    ylab("SPMR")+xlab(paste(chr, "bp"))+
    scale_y_continuous(labels = function(x) round(as.numeric(x), digits=0))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), 
                       axis.line = element_line(colour = "#d9d9d9"),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())+
                       ylim(0, max(as.numeric(spline_int_1_k27ac$y))+0.2)

#  ggsave(file = paste(gene_name, "/treg_actQTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_coverage.pdf", sep=""), plot = q3_2 , width = 15, height = 4, units = "cm",
#           dpi = 300,device = "pdf", useDingbats=FALSE)

  }
if(slope > 0){
  k27acpeaks <- read.table("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/Merged_Day0_treg_ChM_K27ac_sorted-broad_peaks_qval_0.001_FE_2.broadPeak", sep="\t")
  k27acpeaks <- k27acpeaks[ which(k27acpeaks$V1==chr),]
  k27acpeaks <- k27acpeaks[ which(k27acpeaks$V1==chr),]
  k27acpeaks <- k27acpeaks[ which(k27acpeaks$V3 >= from_zoomed),]
  k27acpeaks <- k27acpeaks[ which(k27acpeaks$V2 <= to_zoomed),]

  t2.rect1_k27ac <- data.frame (name= paste(k27acpeaks$V1, k27acpeaks$V2, k27acpeaks$V3, k27acpeaks$V4, sep="_"),xmin=k27acpeaks$V2, xmax=k27acpeaks$V3, ymin=rep(0, length(k27acpeaks$V1)), ymax=rep(Inf, length(k27acpeaks$V1)))
  t2.rect1_k27ac$color <- ifelse(t2.rect1_k27ac$name == k27acpeak, "#F4D7C3", "#d9d9d9" )

  q3_2 = ggplot(coverage8_k27ac_zoomed, aes(coverage8_k27ac_zoomed$new_bin, coverage8_k27ac_zoomed$coverage)) +
    geom_blank() +
    theme_light()+
    geom_rect(data=t2.rect1_k27ac, aes(xmin=as.numeric(xmin), xmax=as.numeric(xmax),
                                 ymin=as.numeric(ymin), ymax=as.numeric(ymax)),
              fill = t2.rect1_k27ac$color, alpha=0.5, inherit.aes = FALSE)+
    geom_line(data = spline_int_2_k27ac, aes(x = x, y = y), colour=c("#9df4e9"))+
    geom_area(data = spline_int_2_k27ac, aes(x = x, y = y),fill=c("#9df4e9"), alpha=0.8) +
    geom_line(data = spline_int_1_k27ac, aes(x = x, y = y), colour=c("#78bbb5"))+
    geom_area(data = spline_int_1_k27ac, aes(x = x, y = y),fill=c("#78bbb5"), alpha=0.8) +
    ylab("SPMR")+xlab(paste(chr, "bp"))+
    scale_y_continuous(labels = function(x) round(as.numeric(x), digits=0))+ 
   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "#d9d9d9"),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())+
                       ylim(0, max(as.numeric(spline_int_2_k27ac$y))+0.2)
  #ggsave(file = paste(gene_name, "/treg_actQTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_coverage.pdf", sep=""), plot = q3_2 , width = 15, height = 4, units = "cm",
  #         dpi = 300,device = "pdf", useDingbats=FALSE)
  }
}
#Load K4me3 peaks
if(k4me3peak != "chrNA"){
  k4me3peaks <- import("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/Merged_Day0_treg_ChM_K4me3_sorted-sharp_peaks_qval_0.001_FE_2.narrowPeak", format = "BED",extraCols = extraCols_narrowPeak)
  k4me3_selected_peaks <- k4me3peaks[k4me3peaks@seqnames == chr & start(k4me3peaks@ranges) > from_zoomed & end(k4me3peaks@ranges) < to_zoomed, ]
  k4me3_peak_list <- split(k4me3_selected_peaks, as.factor(k4me3_selected_peaks))
#get SNP info
  selected_ids_k4me3 <- c(row.names(Major_k4me3QTL),row.names(Minor_k4me3QTL))

  s = subset(snp_table,rownames(snp_table)==snp)
  s = s[ , which(names(s) %in% selected_ids_k4me3)]
  s[s == "0|0"] <- c("1_Hom_Maj")
  s[s == "1|0"] <- c("2_Het")
  s[s == "0|1"] <- c("2_Het")
  s[s == "1|1"] <- c("3_Hom_Min")

  ids <- as.vector(s)
  colnames(ids) <- NULL
  row.names(ids) <- NULL 
  sample_data = dplyr::data_frame(sample_id = names(s), condition = as.character(ids),scaling_factor = 1)
  sample_data = sample_data %>%
    dplyr::mutate(bigWig = paste0("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/BigWig/",sample_id, "_treg_ChM_K4me3_Merged.trimmed.ReAligned_GRCh38.sorted2_uniq_nodup_Autosomes-sharp_SPMR.bw",sep=""))
  as.data.frame(sample_data)
  track_data = dplyr::mutate(sample_data, track_id = condition, colour_group = condition)
  track_data2 = dplyr::mutate(sample_data, track_id = "H3K4me3", colour_group = condition)

### and then to make the plot
#line
  q4 <- plotCoverage2(track_data = track_data2,
                     exons = k4me3_peak_list, 
                     cdss = k4me3_peak_list,
                     #fill_palette = c("#03A461", "#24CB86", "#2BFBA5"),
                     fill_palette = c("#f47362", "#f4adab"),
                     mean_only = TRUE, region_coords = c(from_zoomed,to_zoomed),
                     connect_exons = FALSE,  rescale_introns = FALSE,
                     transcript_label = FALSE, 
                     coverage_type="both", return_subplots_list = TRUE,
                     plot_fraction = 0.005)
  k4me3peaks <- read.table("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/Merged_Day0_treg_ChM_K4me3_sorted-sharp_peaks_qval_0.001_FE_2.narrowPeak", sep="\t")
  k4me3peaks <- k4me3peaks[ which(k4me3peaks$V1==chr),]
  k4me3peaks <- k4me3peaks[ which(k4me3peaks$V1==chr),] 
  k4me3peaks <- k4me3peaks[ which(k4me3peaks$V3 >= from_zoomed),] 
  k4me3peaks <- k4me3peaks[ which(k4me3peaks$V2 <= to_zoomed),] 

  t2.rect2 <- data.frame (name= k4me3peaks$V4,xmin=k4me3peaks$V2, xmax=k4me3peaks$V3, ymin=rep(0, length(k4me3peaks$V1)), ymax=rep(Inf, length(k4me3peaks$V1)))
  t2.rect2$color <- ifelse(t2.rect2$name == k4me3peak, "#F4D7C3", "#d9d9d9")

  q4_2 <- q4$coverage_plot + 
          geom_rect(data=t2.rect2, aes(xmin=as.numeric(xmin), xmax=as.numeric(xmax), 
                                 ymin=as.numeric(ymin), ymax=as.numeric(ymax)), 
                                 fill = t2.rect2$color, alpha=0.5, inherit.aes = FALSE)
    #ggsave(file = paste(gene_name, "/treg_K4me3QTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_coverage.pdf", sep=""), plot = q4_2 , width = 15, height = 4, units = "cm",
    #       dpi = 300,device = "pdf", useDingbats=FALSE)

###Smoothed zoomed coverage plots
coverage <- as.data.frame(q4$bins)
coverage$coverage <- q4$coverage
coverage$colour_group <- q4$colour_group
coverage$sample_id <- q4$sample_id
names(coverage) <- c("bins","coverage", "colour_group", "sample_id")
coverage <- coverage[complete.cases(coverage), ]
coverage_Hom_Maj <- coverage[ which(coverage$colour_group == c("1_Hom_Maj")), ]
coverage_Hom_Min <- coverage[ which(coverage$colour_group == "3_Hom_Min"), ]

coverage2 <- as.data.frame(tapply(coverage_Hom_Maj$coverage, cut(coverage_Hom_Maj$bins, seq(min(coverage_Hom_Maj$bins), max(coverage_Hom_Maj$bins), by=1000)), max))
names(coverage2) <- c("coverage")
coverage3 <- as.data.frame(tapply(coverage_Hom_Maj$bins, cut(coverage_Hom_Maj$bins, seq(min(coverage_Hom_Maj$bins), max(coverage_Hom_Maj$bins), by=1000)), mean))
names(coverage3) <- c("new_bin")
coverage4 <- as.data.frame(tapply(coverage_Hom_Maj$colour_group, cut(coverage_Hom_Maj$bins, seq(min(coverage_Hom_Maj$bins), max(coverage_Hom_Maj$bins), by=1000)), unique))
names(coverage4) <- c("colour_group")
coverage5 <- as.data.frame(tapply(coverage_Hom_Maj$sample_id, cut(coverage_Hom_Maj$bins, seq(min(coverage_Hom_Maj$bins), max(coverage_Hom_Maj$bins), by=1000)), unique))
names(coverage5) <- c("sample_id")
coverage6 <- merge(coverage3,coverage2, by="row.names")
coverage7 <- merge(coverage6,coverage4, by.x = "Row.names", by.y="row.names")
coverage8_Hom_Maj_k4me3 <- merge(coverage7,coverage5, by.x = "Row.names", by.y="row.names")
names(coverage8_Hom_Maj_k4me3) <- c("bin","new_bin","coverage", "colour_group", "sample_id")

coverage2 <- as.data.frame(tapply(coverage_Hom_Min$coverage, cut(coverage_Hom_Min$bins, seq(min(coverage_Hom_Min$bins), max(coverage_Hom_Min$bins), by=1000)), max))
names(coverage2) <- c("coverage")
coverage3 <- as.data.frame(tapply(coverage_Hom_Min$bins, cut(coverage_Hom_Min$bins, seq(min(coverage_Hom_Min$bins), max(coverage_Hom_Min$bins), by=1000)), mean))
names(coverage3) <- c("new_bin")
coverage4 <- as.data.frame(tapply(coverage_Hom_Min$colour_group, cut(coverage_Hom_Min$bins, seq(min(coverage_Hom_Min$bins), max(coverage_Hom_Min$bins), by=1000)), unique))
names(coverage4) <- c("colour_group")
coverage5 <- as.data.frame(tapply(coverage_Hom_Min$sample_id, cut(coverage_Hom_Min$bins, seq(min(coverage_Hom_Min$bins), max(coverage_Hom_Min$bins), by=1000)), unique))
names(coverage5) <- c("sample_id")
coverage6 <- merge(coverage3,coverage2, by="row.names")
coverage7 <- merge(coverage6,coverage4, by.x = "Row.names", by.y="row.names")
coverage8_Hom_Min_k4me3 <- merge(coverage7,coverage5, by.x = "Row.names", by.y="row.names")
names(coverage8_Hom_Min_k4me3) <- c("bin","new_bin","coverage", "colour_group", "sample_id")
coverage8_k4me3_zoomed <- rbind(coverage8_Hom_Maj_k4me3,coverage8_Hom_Min_k4me3)

spline_int_1_k4me3 <- as.data.frame(spline(coverage8_Hom_Maj_k4me3$new_bin, coverage8_Hom_Maj_k4me3$coverage))
spline_int_1_k4me3$y <- ifelse(spline_int_1_k4me3$y < 0, 0, spline_int_1_k4me3$y)
spline_int_2_k4me3 <- as.data.frame(spline(coverage8_Hom_Min_k4me3$new_bin, coverage8_Hom_Min_k4me3$coverage))
spline_int_2_k4me3$y <- ifelse(spline_int_2_k4me3$y < 0, 0, spline_int_2_k4me3$y)
if(slope < 0){
  k4me3peaks <- read.table("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/Merged_Day0_treg_ChM_K4me3_sorted-sharp_peaks_qval_0.001_FE_2.narrowPeak", sep="\t")
  k4me3peaks <- k4me3peaks[ which(k4me3peaks$V1==chr),]
  k4me3peaks <- k4me3peaks[ which(k4me3peaks$V1==chr),]
  k4me3peaks <- k4me3peaks[ which(k4me3peaks$V3 >= from_zoomed),]
  k4me3peaks <- k4me3peaks[ which(k4me3peaks$V2 <= to_zoomed),]
  t2.rect1_k4me3 <- data.frame (name= k4me3peaks$V4,xmin=k4me3peaks$V2, xmax=k4me3peaks$V3, ymin=rep(0, length(k4me3peaks$V1)), ymax=rep(Inf, length(k4me3peaks$V1)))
  t2.rect1_k4me3$color <- ifelse(t2.rect1_k4me3$name == k4me3peak, "#F4D7C3", "#d9d9d9" )

  q4_2 = ggplot(coverage8_k4me3_zoomed, aes(coverage8_k4me3_zoomed$new_bin, coverage8_k4me3_zoomed$coverage)) +
    geom_blank() +
    theme_light()+
    geom_rect(data=t2.rect1_k4me3, aes(xmin=as.numeric(xmin), xmax=as.numeric(xmax),
                                 ymin=as.numeric(ymin), ymax=as.numeric(ymax)),
              fill = t2.rect1_k4me3$color, alpha=0.5, inherit.aes = FALSE)+
    geom_line(data = spline_int_1_k4me3, aes(x = x, y = y), colour=c("#f47362"))+
    geom_area(data = spline_int_1_k4me3, aes(x = x, y = y),fill=c("#f47362"), alpha=0.8) +
    geom_line(data = spline_int_2_k4me3, aes(x = x, y = y), colour=c("#f4adab"))+
    geom_area(data = spline_int_2_k4me3, aes(x = x, y = y),fill=c("#f4adab"), alpha=0.8) +
    ylab("SPMR")+xlab(paste(chr, "bp"))+
    scale_y_continuous(labels = function(x) round(as.numeric(x), digits=0))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "#d9d9d9"),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())+
                       ylim(0, max(as.numeric(spline_int_1_k4me3$y))+0.2)

 #   ggsave(file = paste(gene_name, "/treg_K4me3QTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_coverage.pdf", sep=""), plot = q4_2 , width = 15, height = 4, units = "cm",
 #          dpi = 300,device = "pdf", useDingbats=FALSE)
  
}
if(slope > 0){
  k4me3peaks <- read.table("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/Merged_Day0_treg_ChM_K4me3_sorted-sharp_peaks_qval_0.001_FE_2.narrowPeak", sep="\t")
  k4me3peaks <- k4me3peaks[ which(k4me3peaks$V1==chr),]
  k4me3peaks <- k4me3peaks[ which(k4me3peaks$V1==chr),]
  k4me3peaks <- k4me3peaks[ which(k4me3peaks$V3 >= from_zoomed),]
  k4me3peaks <- k4me3peaks[ which(k4me3peaks$V2 <= to_zoomed),]

  t2.rect1_k4me3 <- data.frame (name= k4me3peaks$V4, xmin=k4me3peaks$V2, xmax=k4me3peaks$V3, ymin=rep(0, length(k4me3peaks$V1)), ymax=rep(Inf, length(k4me3peaks$V1)))
  t2.rect1_k4me3$color <- ifelse(t2.rect1_k4me3$name == k4me3peak, "#F4D7C3", "#d9d9d9" )

  q4_2 = ggplot(coverage8_k4me3_zoomed, aes(coverage8_k4me3_zoomed$new_bin, coverage8_k4me3_zoomed$coverage)) +
    geom_blank() +
    theme_light()+
    geom_rect(data=t2.rect1_k4me3, aes(xmin=as.numeric(xmin), xmax=as.numeric(xmax),
                                 ymin=as.numeric(ymin), ymax=as.numeric(ymax)),
              fill = t2.rect1_k4me3$color, alpha=0.5, inherit.aes = FALSE)+
    geom_line(data = spline_int_2_k4me3, aes(x = x, y = y), colour=c("#f4adab"))+
    geom_area(data = spline_int_2_k4me3, aes(x = x, y = y),fill=c("#f4adab"), alpha=0.8) +
    geom_line(data = spline_int_1_k4me3, aes(x = x, y = y), colour=c("#f47362"))+
    geom_area(data = spline_int_1_k4me3, aes(x = x, y = y),fill=c("#f47362"), alpha=0.8) +
    ylab("SPMR")+xlab(paste(chr, "bp"))+
    scale_y_continuous(labels = function(x) round(as.numeric(x), digits=0))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "#d9d9d9"),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())+
                       ylim(0, max(as.numeric(spline_int_2_k4me3$y))+0.2)

#    ggsave(file = paste(gene_name, "/treg_K4me3QTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_coverage.pdf", sep=""), plot = q4_2 , width = 15, height = 4, units = "cm",
#           dpi = 300,device = "pdf", useDingbats=FALSE)  
  }
}
#Load ATAC peaks
if(atacpeak != "chrNA"){
  atacpeaks <- import("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ATACseq/QTLtools/featureCounts/2MRandomReads_Autosomes_sorted_10Reads58Samples.narrowPeak", format = "BED",extraCols = extraCols_narrowPeak)
  atac_selected_peaks <- atacpeaks[atacpeaks@seqnames == chr & start(atacpeaks@ranges) > from_zoomed & end(atacpeaks@ranges) < to_zoomed, ]
  atac_peak_list <- split(atac_selected_peaks, as.factor(atac_selected_peaks))
#get SNP info
  selected_ids_atac <- c(row.names(Major_atacQTL),row.names(Minor_atacQTL))

  s = subset(snp_table,rownames(snp_table)==snp)
  s = s[ , which(names(s) %in% selected_ids_atac)]
  s[s == "0|0"] <- c("1_Hom_Maj")
  s[s == "1|0"] <- c("2_Het")
  s[s == "0|1"] <- c("2_Het")
  s[s == "1|1"] <- c("3_Hom_Min")

  ids <- as.vector(s)
  colnames(ids) <- NULL
  row.names(ids) <- NULL
  sample_data = dplyr::data_frame(sample_id = names(s), condition = as.character(ids),scaling_factor = 1)
  sample_data = sample_data %>%
    dplyr::mutate(bigWig = paste0("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ATACseq/BigWigs/",sample_id, "_ATAC.bw",sep=""))
  as.data.frame(sample_data)
  track_data = dplyr::mutate(sample_data, track_id = condition, colour_group = condition)
  track_data2 = dplyr::mutate(sample_data, track_id = "ATAC", colour_group = condition)

### and then to make the plot
#line
  q5 <- plotCoverage2(track_data = track_data2,
                     exons = atac_peak_list,
                     cdss = atac_peak_list,
                     #fill_palette = c("#06658C", "#0A9EDB", "#68D1FC"),
                     fill_palette = c("#bb804f", "#efd6c5"),
                     mean_only = TRUE, region_coords = c(from_zoomed,to_zoomed),
                     connect_exons = FALSE,  rescale_introns = FALSE,
                     transcript_label = FALSE,
                     coverage_type="both", return_subplots_list = TRUE,
                     plot_fraction = 0.005)
  atacpeaks <- read.table("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ATACseq/QTLtools/featureCounts/2MRandomReads_Autosomes_sorted_10Reads58Samples.narrowPeak", sep="\t")
  atacpeaks <- atacpeaks[ which(atacpeaks$V1==chr),]
  atacpeaks <- atacpeaks[ which(atacpeaks$V1==chr),]
  atacpeaks <- atacpeaks[ which(atacpeaks$V3 >= from_zoomed),]
  atacpeaks <- atacpeaks[ which(atacpeaks$V2 <= to_zoomed),]

  t2.rect2 <- data.frame (name= atacpeaks$V4,xmin=atacpeaks$V2, xmax=atacpeaks$V3, ymin=rep(0, length(atacpeaks$V1)), ymax=rep(Inf, length(atacpeaks$V1)))
  t2.rect2$color <- ifelse(t2.rect2$name == atacpeak, "#f9f973", "#cccccc")

  q5_2 <- q5$coverage_plot +
          geom_rect(data=t2.rect2, aes(xmin=as.numeric(xmin), xmax=as.numeric(xmax),
                                 ymin=as.numeric(ymin), ymax=as.numeric(ymax)),
                                 fill = t2.rect2$color, alpha=0.3, inherit.aes = FALSE)
#    ggsave(file = paste(gene_name, "/treg_atacQTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_coverage.pdf", sep=""), plot = q5_2 , width = 15, height = 4, units = "cm",
#           dpi = 300,device = "pdf", useDingbats=FALSE)

###Smoothed zoomed coverage plots
coverage <- as.data.frame(q5$bins)
coverage$coverage <- q5$coverage
coverage$colour_group <- q5$colour_group
coverage$sample_id <- q5$sample_id
names(coverage) <- c("bins","coverage", "colour_group", "sample_id")
coverage <- coverage[complete.cases(coverage), ]
coverage_Hom_Maj <- coverage[ which(coverage$colour_group == c("1_Hom_Maj")), ]
coverage_Hom_Min <- coverage[ which(coverage$colour_group == "3_Hom_Min"), ]

coverage2 <- as.data.frame(tapply(coverage_Hom_Maj$coverage, cut(coverage_Hom_Maj$bins, seq(min(coverage_Hom_Maj$bins), max(coverage_Hom_Maj$bins), by=1000)), max))
names(coverage2) <- c("coverage")
coverage3 <- as.data.frame(tapply(coverage_Hom_Maj$bins, cut(coverage_Hom_Maj$bins, seq(min(coverage_Hom_Maj$bins), max(coverage_Hom_Maj$bins), by=1000)), mean))
names(coverage3) <- c("new_bin")
coverage4 <- as.data.frame(tapply(coverage_Hom_Maj$colour_group, cut(coverage_Hom_Maj$bins, seq(min(coverage_Hom_Maj$bins), max(coverage_Hom_Maj$bins), by=1000)), unique))
names(coverage4) <- c("colour_group")
coverage5 <- as.data.frame(tapply(coverage_Hom_Maj$sample_id, cut(coverage_Hom_Maj$bins, seq(min(coverage_Hom_Maj$bins), max(coverage_Hom_Maj$bins), by=1000)), unique))
names(coverage5) <- c("sample_id")
coverage6 <- merge(coverage3,coverage2, by="row.names")
coverage7 <- merge(coverage6,coverage4, by.x = "Row.names", by.y="row.names")
coverage8_Hom_Maj_atac <- merge(coverage7,coverage5, by.x = "Row.names", by.y="row.names")
names(coverage8_Hom_Maj_atac) <- c("bin","new_bin","coverage", "colour_group", "sample_id")

coverage2 <- as.data.frame(tapply(coverage_Hom_Min$coverage, cut(coverage_Hom_Min$bins, seq(min(coverage_Hom_Min$bins), max(coverage_Hom_Min$bins), by=1000)), max))
names(coverage2) <- c("coverage")
coverage3 <- as.data.frame(tapply(coverage_Hom_Min$bins, cut(coverage_Hom_Min$bins, seq(min(coverage_Hom_Min$bins), max(coverage_Hom_Min$bins), by=1000)), mean))
names(coverage3) <- c("new_bin")
coverage4 <- as.data.frame(tapply(coverage_Hom_Min$colour_group, cut(coverage_Hom_Min$bins, seq(min(coverage_Hom_Min$bins), max(coverage_Hom_Min$bins), by=1000)), unique))
names(coverage4) <- c("colour_group")
coverage5 <- as.data.frame(tapply(coverage_Hom_Min$sample_id, cut(coverage_Hom_Min$bins, seq(min(coverage_Hom_Min$bins), max(coverage_Hom_Min$bins), by=1000)), unique))
names(coverage5) <- c("sample_id")
coverage6 <- merge(coverage3,coverage2, by="row.names")
coverage7 <- merge(coverage6,coverage4, by.x = "Row.names", by.y="row.names")
coverage8_Hom_Min_atac <- merge(coverage7,coverage5, by.x = "Row.names", by.y="row.names")
names(coverage8_Hom_Min_atac) <- c("bin","new_bin","coverage", "colour_group", "sample_id")
coverage8_atac_zoomed <- rbind(coverage8_Hom_Maj_atac,coverage8_Hom_Min_atac)


spline_int_1_atac <- as.data.frame(spline(coverage8_Hom_Maj_atac$new_bin, coverage8_Hom_Maj_atac$coverage))
spline_int_1_atac$y <- ifelse(spline_int_1_atac$y < 0, 0, spline_int_1_atac$y)
spline_int_2_atac <- as.data.frame(spline(coverage8_Hom_Min_atac$new_bin, coverage8_Hom_Min_atac$coverage))
spline_int_2_atac$y <- ifelse(spline_int_2_atac$y < 0, 0, spline_int_2_atac$y)
if(slope < 0){
  atacpeaks <- read.table("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ATACseq/QTLtools/featureCounts/2MRandomReads_Autosomes_sorted_10Reads58Samples.narrowPeak", sep="\t")
  atacpeaks <- atacpeaks[ which(atacpeaks$V1==chr),]
  atacpeaks <- atacpeaks[ which(atacpeaks$V1==chr),]
  atacpeaks <- atacpeaks[ which(atacpeaks$V3 >= from_zoomed),]
  atacpeaks <- atacpeaks[ which(atacpeaks$V2 <= to_zoomed),]

  t2.rect1_atac <- data.frame (name= atacpeaks$V4 ,xmin=atacpeaks$V2, xmax=atacpeaks$V3, ymin=rep(0, length(atacpeaks$V1)), ymax=rep(Inf, length(atacpeaks$V1)))
  t2.rect1_atac$color <- ifelse(t2.rect1_atac$name == atacpeak, "#F4D7C3", "#d9d9d9" )


  q5_2 = ggplot(coverage8_atac_zoomed, aes(coverage8_atac_zoomed$new_bin, coverage8_atac_zoomed$coverage)) +
    geom_blank() +
    theme_light()+
    geom_rect(data=t2.rect1_atac, aes(xmin=as.numeric(xmin), xmax=as.numeric(xmax),
                                 ymin=as.numeric(ymin), ymax=as.numeric(ymax)),
              fill = t2.rect1_atac$color, alpha=0.5, inherit.aes = FALSE)+
    geom_line(data = spline_int_1_atac, aes(x = x, y = y), colour=c("#bb804f"))+
    geom_area(data = spline_int_1_atac, aes(x = x, y = y),fill=c("#bb804f"), alpha=0.8) +
    geom_line(data = spline_int_2_atac, aes(x = x, y = y), colour=c("#efd6c5"))+
    geom_area(data = spline_int_2_atac, aes(x = x, y = y),fill=c("#efd6c5"), alpha=0.8) +
    ylab("SPMR")+xlab(paste(chr, "bp"))+
    scale_y_continuous(labels = function(x) round(as.numeric(x), digits=0))+
    scale_y_continuous(labels = function(x) round(as.numeric(x), digits=0))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "#d9d9d9"),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())+
                       ylim(0, max(as.numeric(spline_int_1_atac$y))+0.2)

    #ggsave(file = paste(gene_name, "/treg_atacQTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_coverage.pdf", sep=""), plot = q5_2 , width = 15, height = 4, units = "cm",
    #       dpi = 300,device = "pdf", useDingbats=FALSE)
  }
if(slope > 0){
  atacpeaks <- read.table("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ATACseq/QTLtools/featureCounts/2MRandomReads_Autosomes_sorted_10Reads58Samples.narrowPeak", sep="\t")
  atacpeaks <- atacpeaks[ which(atacpeaks$V1==chr),]
  atacpeaks <- atacpeaks[ which(atacpeaks$V1==chr),]
  atacpeaks <- atacpeaks[ which(atacpeaks$V3 >= from_zoomed),]
  atacpeaks <- atacpeaks[ which(atacpeaks$V2 <= to_zoomed),]

  t2.rect1_atac <- data.frame (name= atacpeaks$V4, xmin=atacpeaks$V2, xmax=atacpeaks$V3, ymin=rep(0, length(atacpeaks$V1)), ymax=rep(Inf, length(atacpeaks$V1)))
  t2.rect1_atac$color <- ifelse(t2.rect1_atac$name == atacpeak, "#F4D7C3", "#d9d9d9" )

  q5_2 = ggplot(coverage8_atac_zoomed, aes(coverage8_atac_zoomed$new_bin, coverage8_atac_zoomed$coverage)) +
    geom_blank() +
    theme_light()+
    geom_rect(data=t2.rect1_atac, aes(xmin=as.numeric(xmin), xmax=as.numeric(xmax),
                                 ymin=as.numeric(ymin), ymax=as.numeric(ymax)),
              fill = t2.rect1_atac$color, alpha=0.5, inherit.aes = FALSE)+
    geom_line(data = spline_int_2_atac, aes(x = x, y = y), colour=c("#efd6c5"))+
    geom_area(data = spline_int_2_atac, aes(x = x, y = y),fill=c("#efd6c5"), alpha=0.8) +
    geom_line(data = spline_int_1_atac, aes(x = x, y = y), colour=c("#bb804f"))+
    geom_area(data = spline_int_1_atac, aes(x = x, y = y),fill=c("#bb804f"), alpha=0.8) +
    ylab("SPMR")+xlab(paste(chr, "bp"))+
    scale_y_continuous(labels = function(x) round(as.numeric(x), digits=0))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "#d9d9d9"),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())+
                       ylim(0, max(as.numeric(spline_int_2_atac$y))+0.2)

    #ggsave(file = paste(gene_name, "/treg_atacQTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_coverage.pdf", sep=""), plot = q5_2 , width = 15, height = 4, units = "cm",
    #       dpi = 300,device = "pdf", useDingbats=FALSE)
  }
}
##combined zoomed coverage plots
if(atacpeak != "chrNA" & k4me3peak != "chrNA" & k27acpeak != "chrNA"){
  allcor <- arrangeGrob(arrangeGrob(q3_2), arrangeGrob(q4_2), arrangeGrob(q5_2),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_coverage.pdf", sep =""), plot = allcor, width = 15, height = 4*3, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(atacpeak == "chrNA" & k4me3peak != "chrNA" & k27acpeak != "chrNA"){
  allcor <- arrangeGrob(arrangeGrob(q3_2), arrangeGrob(q4_2),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_coverage.pdf", sep =""), plot = allcor, width = 15, height = 4*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(atacpeak != "chrNA" & k4me3peak == "chrNA" & k27acpeak != "chrNA"){
  allcor <- arrangeGrob(arrangeGrob(q3_2), arrangeGrob(q5_2),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_coverage.pdf", sep =""), plot = allcor, width = 15, height = 4*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(atacpeak != "chrNA" & k4me3peak != "chrNA" & k27acpeak == "chrNA"){
  allcor <- arrangeGrob(arrangeGrob(q4_2), arrangeGrob(q5_2),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_coverage.pdf", sep =""), plot = allcor, width = 15, height = 4*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(atacpeak != "chrNA" & k4me3peak == "chrNA" & k27acpeak == "chrNA"){
  allcor <- arrangeGrob(arrangeGrob(q5_2),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_coverage.pdf", sep =""), plot = allcor, width = 15, height = 4, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(atacpeak == "chrNA" & k4me3peak != "chrNA" & k27acpeak == "chrNA"){
  allcor <- arrangeGrob(arrangeGrob(q4_2),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_coverage.pdf", sep =""), plot = allcor, width = 15, height = 4, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(atacpeak == "chrNA" & k4me3peak == "chrNA" & k27acpeak != "chrNA"){
  allcor <- arrangeGrob(arrangeGrob(q3_2),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_coverage.pdf", sep =""), plot = allcor, width = 15, height = 4, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

#Create Manhattan plot with for p-values disease and different QTLs
if(k27acpeak != "chrNA"){
  k27acQTL <- read.table(file = paste("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/QTLtools/Ultimate/K27ac/featureCounts/500kb/output/ChM_K27ac_Peak_featureCounts_NoBatchCorrection_CQN.500kb.nominal.sorted_",chr, ".txt", sep =""), sep ="\t", header = FALSE)
  k27acQTL <- k27acQTL[ which(k27acQTL$V1==k27acpeak),] 
  k27acQTL <- k27acQTL[ which(k27acQTL$V10 >= from_zoomed),] 
  k27acQTL <- k27acQTL[ which(k27acQTL$V10 <= to_zoomed),]
  k27acQTL$V12 <- -log10(k27acQTL$V12)
  k27acQTL$GWAS_lead <- ifelse(k27acQTL$V8 == snp, 1 , 0)
  names(k27acQTL) <- c("k27ac_peak", "chr","start", "end", "strand", "variants", "distance", "snp", "chr2","snp_bp", "snp_bp2", "nominal", "slope", "lead_k27acQTL", "GWAS_lead")

##Add LD info
  missingLDblock <- LD4[-which(LD4$snp2 %in% k27acQTL$snp),]
  missingLDblock <- missingLDblock[, names(missingLDblock) %in% c("snp2", "bp2")]
  names(missingLDblock) <- c("snp_bp", "snp")
  missingLDblock$k27ac_peak <- rep(k27acpeak,length(missingLDblock$snp))
  missingLDblock$chr <- rep(chr,length(missingLDblock$snp))
  missingLDblock$start <- rep(0,length(missingLDblock$snp))
  missingLDblock$end <- rep(0,length(missingLDblock$snp))
  missingLDblock$strand <- rep("+",length(missingLDblock$snp))
  missingLDblock$variants <- rep(0,length(missingLDblock$snp))
  missingLDblock$distance <- rep(0,length(missingLDblock$snp))
  missingLDblock$chr2 <- rep(chr,length(missingLDblock$snp))
  missingLDblock$snp_bp2 <- rep(0,length(missingLDblock$snp))
  missingLDblock$nominal <- rep(0,length(missingLDblock$snp))
  missingLDblock$slope <- rep(0,length(missingLDblock$snp))
  missingLDblock$lead_k27acQTL <- rep(0,length(missingLDblock$snp))
  missingLDblock$GWAS_lead <- rep(0,length(missingLDblock$snp))
  missingLDblock <- missingLDblock[,c("k27ac_peak", "chr","start", "end", "strand", "variants", "distance", "snp", "chr2","snp_bp", "snp_bp2", "nominal", "slope", "lead_k27acQTL", "GWAS_lead")]
  k27acQTL <- rbind(k27acQTL, missingLDblock)

  LD <- read.table(paste(snp, ".ld", sep=""), sep="\t")
  names(LD) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
  Add_autoLD <- c(chr, snp_bp, snp, chr, snp_bp, snp, "1", "NA")
  LD <- rbind(LD, Add_autoLD)
  LD2 <- LD[ which(LD$bp1==snp_bp),]
  LD3 <- LD[ which(LD$bp2==snp_bp),]
  LD3 <- LD3[,c(1,5,6,4,2,3,7,8)]
  names(LD3) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
  LD5 <- rbind(LD2,LD3)
  rm(LD2,LD3)
  LD2 <- merge(k27acQTL,LD5, by.x = "snp_bp", by.y = "bp2", all.x = TRUE) 
  LD2 <- unique(LD2)
  LD2$r2 <- as.numeric(LD2$r2)
  LD2$GWAS_lead <- ifelse(LD2$snp_bp == snp_bp, 1, 0)
  LD2_k27ac <- LD2
  LD2_k27ac$r2[is.na(LD2_k27ac$r2)] <- 0


  pval1_2 <- ggplot(LD2_k27ac, aes(x=as.numeric(LD2_k27ac$snp_bp), y=as.numeric(LD2_k27ac$nominal), 
                           shape=as.factor(LD2_k27ac$GWAS_lead)))+
    annotate("rect", xmin=as.numeric(k27acpeak_region[2]), xmax=as.numeric(k27acpeak_region[3]), ymin=0, ymax=Inf, alpha=0.5, fill="#F4D7C3")+
    geom_point(aes(fill = cut(as.numeric(LD2_k27ac$r2), c(0, 0.2, 0.5, 0.8, 1), include.lowest = TRUE)),
               size = LD2_k27ac$lead_k27acQTL+3, colour =  "#bebebe") +
    scale_shape_manual(values=c(21, 23))+
    scale_fill_manual(name = paste("r2 with ",snp, sep=""),
                      values = c("[0,0.2]" = "#f9f9f9",
                                 "(0.2,0.5]" = "#34495e",
                                 "(0.5,0.8]" = "#3498db",
                                 "(0.8,1]" = "#9b59b6"),
                      labels = c("<= 0.2", "0.2 < r2 <= 0.5", "0.5 < r2 <= 0.8","> 0.8"))+
    coord_cartesian(xlim = c(from_zoomed, to_zoomed))+
    labs(fill = paste("R2 with ", snp), sep = "") +
    xlab(paste(chr, " bp", sep = ""))+ ylab("-LOG10(actQTL p-value)")+
    theme_classic()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
           legend.position = "none")

#  ggsave(paste(gene_name, "/treg_K27acQTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_pvalue.pdf", sep =""), plot = pval1_2, width = 15, height = 5, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(k4me3peak != "chrNA"){
  k4me3QTL <- read.table(file = paste("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/QTLtools/Ultimate/K4me3/featureCounts/500kb/output/ChM_K4me3_Peak_featureCounts_NoBatchCorrection_CQN.nominal.sorted_",gsub("chr","",chr), ".txt", sep =""), sep ="\t", header = FALSE)
  k4me3QTL <- k4me3QTL[ which(k4me3QTL$V1==k4me3peak),]
  k4me3QTL <- k4me3QTL[ which(k4me3QTL$V10 >= from_zoomed),]
  k4me3QTL <- k4me3QTL[ which(k4me3QTL$V10 <= to_zoomed),]
  k4me3QTL$V12 <- -log10(k4me3QTL$V12)
  k4me3QTL$GWAS_lead <- ifelse(k4me3QTL$V8 == snp, 1 , 0)
  names(k4me3QTL) <- c("k4me3_peak", "chr","start", "end", "strand", "variants", "distance", "snp", "chr2","snp_bp", "snp_bp2", "nominal", "slope", "lead_k4me3QTL", "GWAS_lead")

##Add LD info
  missingLDblock <- LD4[-which(LD4$snp2 %in% k4me3QTL$snp),]
  missingLDblock <- missingLDblock[, names(missingLDblock) %in% c("snp2", "bp2")]
  names(missingLDblock) <- c("snp_bp", "snp")
  missingLDblock$k4me3_peak <- rep(k4me3peak,length(missingLDblock$snp))
  missingLDblock$chr <- rep(chr,length(missingLDblock$snp))
  missingLDblock$start <- rep(0,length(missingLDblock$snp))
  missingLDblock$end <- rep(0,length(missingLDblock$snp))
  missingLDblock$strand <- rep("+",length(missingLDblock$snp))
  missingLDblock$variants <- rep(0,length(missingLDblock$snp))
  missingLDblock$distance <- rep(0,length(missingLDblock$snp))
  missingLDblock$chr2 <- rep(chr,length(missingLDblock$snp))
  missingLDblock$snp_bp2 <- rep(0,length(missingLDblock$snp))
  missingLDblock$nominal <- rep(0,length(missingLDblock$snp))
  missingLDblock$slope <- rep(0,length(missingLDblock$snp))
  missingLDblock$lead_k4me3QTL <- rep(0,length(missingLDblock$snp))
  missingLDblock$GWAS_lead <- rep(0,length(missingLDblock$snp))
  missingLDblock <- missingLDblock[,c("k4me3_peak", "chr","start", "end", "strand", "variants", "distance", "snp", "chr2","snp_bp", "snp_bp2", "nominal", "slope", "lead_k4me3QTL", "GWAS_lead")]
  k4me3QTL <- rbind(k4me3QTL, missingLDblock)

  LD <- read.table(paste(snp, ".ld", sep=""), sep="\t")
  names(LD) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
  Add_autoLD <- c(chr, snp_bp, snp, chr, snp_bp, snp, "1", "NA")
  LD <- rbind(LD, Add_autoLD)
  LD2 <- LD[ which(LD$bp1==snp_bp),]
  LD3 <- LD[ which(LD$bp2==snp_bp),]
  LD3 <- LD3[,c(1,5,6,4,2,3,7,8)]
  names(LD3) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
  LD5 <- rbind(LD2,LD3)
  rm(LD2,LD3)
  LD2 <- merge(k4me3QTL,LD5, by.x = "snp_bp", by.y = "bp2", all.x = TRUE)
  LD2 <- unique(LD2)
  LD2$r2 <- as.numeric(LD2$r2)
  LD2$GWAS_lead <- ifelse(LD2$snp_bp == snp_bp, 1, 0)
  LD2_k4me3 <- LD2
  LD2_k4me3$r2[is.na(LD2_k4me3$r2)] <- 0

  pval1_4 <- ggplot(LD2_k4me3, aes(x=as.numeric(LD2_k4me3$snp_bp), y=as.numeric(LD2_k4me3$nominal),
                           shape=as.factor(LD2_k4me3$GWAS_lead)))+
    annotate("rect", xmin=as.numeric(k4me3peak_region[2]), xmax=as.numeric(k4me3peak_region[3]), ymin=0, ymax=Inf, alpha=0.5, fill="#F4D7C3")+
    geom_point(aes(fill = cut(as.numeric(LD2_k4me3$r2), c(0, 0.2, 0.5, 0.8, 1), include.lowest = TRUE)),
               size = LD2_k4me3$lead_k4me3QTL+3, colour =  "#bebebe") +
    scale_shape_manual(values=c(21, 23))+
    scale_fill_manual(name = paste("r2 with ",snp, sep=""),
                      values = c("[0,0.2]" = "#f9f9f9",
                                 "(0.2,0.5]" = "#34495e",
                                 "(0.5,0.8]" = "#3498db",
                                 "(0.8,1]" = "#9b59b6"),
                      labels = c("<= 0.2", "0.2 < r2 <= 0.5", "0.5 < r2 <= 0.8","> 0.8"))+
    coord_cartesian(xlim = c(from_zoomed, to_zoomed))+
    labs(fill = paste("R2 with ", snp), sep = "") +
    xlab(paste(chr, " bp", sep = ""))+ ylab("-LOG10(H3K4me3QTL p-value)")+
    theme_classic()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
           legend.position = "none")

#  ggsave(paste(gene_name, "/treg_K4me3QTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_pvalue.pdf", sep =""), plot = pval1_4, width = 15, height = 5, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(atacpeak != "chrNA"){
  atacQTL <- read.table(file = paste("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ATACseq/QTLtools/featureCounts/500kb/output/ATAC_Counts_featureCounts_NoBatchCorrection_CQN.nominal.sorted_",gsub("chr","",chr), ".txt", sep =""), sep ="\t", header = FALSE)
  atacQTL <- atacQTL[ which(atacQTL$V1==atacpeak),]
  atacQTL <- atacQTL[ which(atacQTL$V10 >= from_zoomed),]
  atacQTL <- atacQTL[ which(atacQTL$V10 <= to_zoomed),]
  atacQTL$V12 <- -log10(atacQTL$V12)
  atacQTL$GWAS_lead <- ifelse(atacQTL$V8 == snp, 1 , 0)
  names(atacQTL) <- c("atac_peak", "chr","start", "end", "strand", "variants", "distance", "snp", "chr2","snp_bp", "snp_bp2", "nominal", "slope", "lead_atacQTL", "GWAS_lead")

##Add LD info
  missingLDblock <- LD4[-which(LD4$snp2 %in% atacQTL$snp),]
  missingLDblock <- missingLDblock[, names(missingLDblock) %in% c("snp2", "bp2")]
  names(missingLDblock) <- c("snp_bp", "snp")
  missingLDblock$atac_peak <- rep(atacpeak,length(missingLDblock$snp))
  missingLDblock$chr <- rep(chr,length(missingLDblock$snp))
  missingLDblock$start <- rep(0,length(missingLDblock$snp))
  missingLDblock$end <- rep(0,length(missingLDblock$snp))
  missingLDblock$strand <- rep("+",length(missingLDblock$snp))
  missingLDblock$variants <- rep(0,length(missingLDblock$snp))
  missingLDblock$distance <- rep(0,length(missingLDblock$snp))
  missingLDblock$chr2 <- rep(chr,length(missingLDblock$snp))
  missingLDblock$snp_bp2 <- rep(0,length(missingLDblock$snp))
  missingLDblock$nominal <- rep(0,length(missingLDblock$snp))
  missingLDblock$slope <- rep(0,length(missingLDblock$snp))
  missingLDblock$lead_atacQTL <- rep(0,length(missingLDblock$snp))
  missingLDblock$GWAS_lead <- rep(0,length(missingLDblock$snp))
  missingLDblock <- missingLDblock[,c("atac_peak", "chr","start", "end", "strand", "variants", "distance", "snp", "chr2","snp_bp", "snp_bp2", "nominal", "slope", "lead_atacQTL", "GWAS_lead")]
  atacQTL <- rbind(atacQTL, missingLDblock)

  LD <- read.table(paste(snp, ".ld", sep=""), sep="\t")
  names(LD) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
  Add_autoLD <- c(chr, snp_bp, snp, chr, snp_bp, snp, "1", "NA")
  LD <- rbind(LD, Add_autoLD)
  LD2 <- LD[ which(LD$bp1==snp_bp),]
  LD3 <- LD[ which(LD$bp2==snp_bp),]
  LD3 <- LD3[,c(1,5,6,4,2,3,7,8)]
  names(LD3) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
  LD5 <- rbind(LD2,LD3)
  rm(LD2,LD3)
  LD2 <- merge(atacQTL,LD5, by.x = "snp_bp", by.y = "bp2", all.x = TRUE)
  LD2 <- unique(LD2)
  LD2$r2 <- as.numeric(LD2$r2)
  LD2$GWAS_lead <- ifelse(LD2$snp_bp == snp_bp, 1, 0)
  LD2_atac <- LD2
  LD2_atac$r2[is.na(LD2_atac$r2)] <- 0

  pval1_5 <- ggplot(LD2_atac, aes(x=as.numeric(LD2_atac$snp_bp), y=as.numeric(LD2_atac$nominal),
                           shape=as.factor(LD2_atac$GWAS_lead)))+
    annotate("rect", xmin=as.numeric(atacpeak_region[2]), xmax=as.numeric(atacpeak_region[3]), ymin=0, ymax=Inf, alpha=0.5, fill="#F4D7C3")+
    geom_point(aes(fill = cut(as.numeric(LD2_atac$r2), c(0, 0.2, 0.5, 0.8, 1), include.lowest = TRUE)),
               size = LD2_atac$lead_atacQTL+3, colour =  "#bebebe") +
    scale_shape_manual(values=c(21, 23))+
    scale_fill_manual(name = paste("r2 with ",snp, sep=""),
                      values = c("[0,0.2]" = "#f9f9f9",
                                 "(0.2,0.5]" = "#34495e",
                                 "(0.5,0.8]" = "#3498db",
                                 "(0.8,1]" = "#9b59b6"),
                      labels = c("<= 0.2", "0.2 < r2 <= 0.5", "0.5 < r2 <= 0.8","> 0.8"))+
    coord_cartesian(xlim = c(from_zoomed, to_zoomed))+
    labs(fill = paste("R2 with ", snp), sep = "") +
    xlab(paste(chr, " bp", sep = ""))+ ylab("-LOG10(atacQTL p-value)")+
    theme_classic()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
           legend.position = "none")

#  ggsave(paste(gene_name, "/treg_atacQTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_pvalue.pdf", sep =""), plot = pval1_5, width = 15, height = 5, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)
}

#Create Manhattan plot with for eQTLs
if(gene != "Peak"){
  eQTL <- read.table(file = paste("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_RNAseq/Dafnis_Results/day0.nominal.sorted_", chr, ".txt", sep = ""), header = FALSE)
  eQTL <- eQTL[ which(eQTL$V1==gene),] 
  eQTL <- eQTL[ which(eQTL$V10 >= from_zoomed),] 
  eQTL <- eQTL[ which(eQTL$V10 <= to_zoomed),]
  eQTL$V12 <- -log10(eQTL$V12)
  eQTL$GWAS_lead <- ifelse(eQTL$V8 == snp, 1 , 0)
  names(eQTL) <- c("gene", "chr","start", "end", "strand", "variants", "distance", "snp", "chr2","snp_bp", "snp_bp2", "nominal", "slope", "lead_eQTL", "GWAS_lead")

##Add LD info
  missingLDblock <- LD4[-which(LD4$snp2 %in% eQTL$snp),]
  missingLDblock <- missingLDblock[, names(missingLDblock) %in% c("snp2", "bp2")]
  names(missingLDblock) <- c("snp_bp", "snp")
  missingLDblock$gene <- rep(gene,length(missingLDblock$snp))
  missingLDblock$chr <- rep(chr,length(missingLDblock$snp))
  missingLDblock$start <- rep(0,length(missingLDblock$snp))
  missingLDblock$end <- rep(0,length(missingLDblock$snp))
  missingLDblock$strand <- rep("+",length(missingLDblock$snp))
  missingLDblock$variants <- rep(0,length(missingLDblock$snp))
  missingLDblock$distance <- rep(0,length(missingLDblock$snp))
  missingLDblock$chr2 <- rep(chr,length(missingLDblock$snp))
  missingLDblock$snp_bp2 <- rep(0,length(missingLDblock$snp))
  missingLDblock$nominal <- rep(0,length(missingLDblock$snp))
  missingLDblock$slope <- rep(0,length(missingLDblock$snp))
  missingLDblock$lead_eQTL <- rep(0,length(missingLDblock$snp))
  missingLDblock$GWAS_lead <- rep(0,length(missingLDblock$snp))
  missingLDblock <- missingLDblock[,c("gene", "chr","start", "end", "strand", "variants", "distance", "snp", "chr2","snp_bp", "snp_bp2", "nominal", "slope", "lead_eQTL", "GWAS_lead")]
  eQTL <- rbind(eQTL, missingLDblock)

  LD <- read.table(paste(snp, ".ld", sep=""), sep="\t")
  names(LD) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
  Add_autoLD <- c(chr, snp_bp, snp, chr, snp_bp, snp, "1", "NA")
  LD <- rbind(LD, Add_autoLD)
  LD2 <- LD[ which(LD$bp1==snp_bp),]
  LD3 <- LD[ which(LD$bp2==snp_bp),]
  LD3 <- LD3[,c(1,5,6,4,2,3,7,8)]
  names(LD3) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
  LD5 <- rbind(LD2,LD3)
  rm(LD2,LD3)
  LD2 <- merge(eQTL,LD5, by.x = "snp_bp", by.y = "bp2", all.x = TRUE) 
  LD2 <- unique(LD2)
  LD2$r2 <- as.numeric(LD2$r2)
  LD2$GWAS_lead <- ifelse(LD2$snp_bp == snp_bp, 1, 0)
  LD2_rna <- LD2
  LD2_rna$r2[is.na(LD2_rna$r2)] <- 0

if(k27acpeak != "chrNA"){

  pval1_3 <- ggplot(LD2_rna, aes(x=as.numeric(LD2_rna$snp_bp), y=as.numeric(LD2_rna$nominal), 
                             shape=as.factor(LD2_rna$GWAS_lead)))+
    annotate("rect", xmin=as.numeric(k27acpeak_region[2]), xmax=as.numeric(k27acpeak_region[3]), ymin=0, ymax=Inf, alpha=0.5, fill="#F4D7C3")+
    geom_point(aes(fill = cut(as.numeric(LD2_rna$r2), c(0, 0.2, 0.5, 0.8, 1), include.lowest = TRUE)),
               size = LD2_rna$lead_eQTL+3, colour =  "#bebebe") +
    scale_shape_manual(values=c(21, 23))+
    scale_fill_manual(name = paste("r2 with ",snp, sep=""),
                      values = c("[0,0.2]" = "#f9f9f9",
                                 "(0.2,0.5]" = "#34495e",
                                 "(0.5,0.8]" = "#3498db",
                                 "(0.8,1]" = "#9b59b6"),
                      labels = c("<= 0.2", "0.2 < r2 <= 0.5", "0.5 < r2 <= 0.8","> 0.8"))+
    coord_cartesian(xlim = c(from_zoomed, to_zoomed))+
    labs(fill = paste("R2 with ", snp), sep = "") +
    xlab(paste(chr, " bp", sep = ""))+ ylab("-LOG10(eQTL p-value)")+
    theme_classic()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          legend.position = "none")

#  ggsave(paste(gene_name, "/treg_eQTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_pvalue.pdf", sep =""), plot = pval1_3, width = 15, height = 5, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)
}
if(k4me3peak != "chrNA" & k27acpeak == "chrNA"){

  pval1_3 <- ggplot(LD2_rna, aes(x=LD2_rna$snp_bp, y=LD2_rna$nominal, 
                             shape=as.factor(LD2_rna$GWAS_lead)))+
    annotate("rect", xmin=as.numeric(k4me3peak_region[2]), xmax=as.numeric(k4me3peak_region[3]), ymin=0, ymax=Inf, alpha=0.5, fill="#F4D7C3")+
    geom_point(aes(fill = cut(as.numeric(LD2_rna$r2), c(0, 0.2, 0.5, 0.8, 1), include.lowest = TRUE)),
               size = LD2_rna$lead_eQTL+3, colour =  "#bebebe") +
    scale_shape_manual(values=c(21, 23))+
    scale_fill_manual(name = paste("r2 with ",snp, sep=""),
                      values = c("[0,0.2]" = "#f9f9f9",
                                 "(0.2,0.5]" = "#34495e",
                                 "(0.5,0.8]" = "#3498db",
                                 "(0.8,1]" = "#9b59b6"),
                      labels = c("<= 0.2", "0.2 < r2 <= 0.5", "0.5 < r2 <= 0.8","> 0.8"))+
    coord_cartesian(xlim = c(from_zoomed, to_zoomed))+
    labs(fill = paste("R2 with ", snp), sep = "") +
    xlab(paste(chr, " bp", sep = ""))+ ylab("-LOG10(eQTL p-value)")+
    theme_classic()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          legend.position = "none")

#  ggsave(paste(gene_name, "/treg_eQTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_pvalue.pdf", sep =""), plot = pval1_3, width = 15, height = 5, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)
}
if(atacpeak != "chrNA"  & k27acpeak == "chrNA" & k4me3peak == "chrNA"){

  pval1_3 <- ggplot(LD2_rna, aes(x=as.numeric(LD2_rna$snp_bp), y=as.numeric(LD2_rna$nominal), 
                             shape=as.factor(LD2_rna$GWAS_lead)))+
    annotate("rect", xmin=as.numeric(atacpeak_region[2]), xmax=as.numeric(atacpeak_region[3]), ymin=0, ymax=Inf, alpha=0.5, fill="#F4D7C3")+
    geom_point(aes(fill = cut(as.numeric(LD2_rna$r2), c(0, 0.2, 0.5, 0.8, 1), include.lowest = TRUE)),
               size = LD2_rna$lead_eQTL+3, colour =  "#bebebe") +
    scale_shape_manual(values=c(21, 23))+
    scale_fill_manual(name = paste("r2 with ",snp, sep=""),
                      values = c("[0,0.2]" = "#f9f9f9",
                                 "(0.2,0.5]" = "#34495e",
                                 "(0.5,0.8]" = "#3498db",
                                 "(0.8,1]" = "#9b59b6"),
                      labels = c("<= 0.2", "0.2 < r2 <= 0.5", "0.5 < r2 <= 0.8","> 0.8"))+
    coord_cartesian(xlim = c(from_zoomed, to_zoomed))+
    labs(fill = paste("R2 with ", snp), sep = "") +
    xlab(paste(chr, " bp", sep = ""))+ ylab("-LOG10(eQTL p-value)")+
    theme_classic()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          legend.position = "none")

#  ggsave(paste(gene_name, "/treg_eQTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_pvalue.pdf", sep =""), plot = pval1_3, width = 15, height = 5, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)
}
if(k27acpeak == "chrNA"  & k4me3peak == "chrNA" &  atacpeak == "chrNA"){

  pval1_3 <- ggplot(LD2_rna, aes(x=as.numeric(LD2_rna$snp_bp), y=as.numeric(LD2_rna$nominal), 
                             shape=as.factor(LD2_rna$GWAS_lead)))+
    geom_point(aes(fill = cut(as.numeric(LD2_rna$r2), c(0, 0.2, 0.5, 0.8, 1), include.lowest = TRUE)),
               size = LD2_rna$lead_eQTL+3, colour =  "#bebebe") +
    scale_shape_manual(values=c(21, 23))+
    scale_fill_manual(name = paste("r2 with ",snp, sep=""),
                      values = c("[0,0.2]" = "#f9f9f9",
                                 "(0.2,0.5]" = "#34495e",
                                 "(0.5,0.8]" = "#3498db",
                                 "(0.8,1]" = "#9b59b6"),
                      labels = c("<= 0.2", "0.2 < r2 <= 0.5", "0.5 < r2 <= 0.8","> 0.8"))+
    coord_cartesian(xlim = c(from_zoomed, to_zoomed))+
    labs(fill = paste("R2 with ", snp), sep = "") +
    xlab(paste(chr, " bp", sep = ""))+ ylab("-LOG10(eQTL p-value)")+
    theme_classic()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          legend.position = "none")

#  ggsave(paste(gene_name, "/treg_eQTL_",gene_name, "_", gene, "_", snp, "_","_Zoomed_pvalue.pdf", sep =""), plot = pval1_3, width = 15, height = 5, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)
}

}

#Plot GWAS association p-values

GWAS1 <- read.table(GWAS, sep="\t", header = TRUE)
GWAS1$Chr <- paste("chr", GWAS1$Chr, sep="")
GWAS1 <- GWAS1[which(GWAS1$Chr == chr),]
GWAS1 <- GWAS1[which(GWAS1$Pos >= snp_bp - 5000000),]
GWAS1 <- GWAS1[which(GWAS1$Pos <= snp_bp + 5000000),]
GWAS1$RSid <- paste("chr", GWAS1$RSid, sep="")
GWAS1$RSid <- gsub(':', '_', GWAS1$RSid)
GWAS1$pval <- -log10(GWAS1$pval)
GWAS1_1 <- GWAS1
for(i in 1:length(GWAS1_1$RSid)){
  GWAS1_1$kk[i] <- as.character(as.vector(unlist(strsplit(GWAS1_1$RSid[i], "_")))[1])
  GWAS1_1$kk2[i] <- as.character(as.vector(unlist(strsplit(GWAS1_1$RSid[i], "_")))[2])
  GWAS1_1$RSid[i] <- paste(GWAS1_1$kk[i], GWAS1_1$kk2[i],sep="_")
} 
##Add GRCh38 position

GRCh38_bp <- read.table(paste("/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_Genotypes/BLUEPRINT_mimic_Imputation/1KGenomes_10kUK_Reference/UK10K_",  gsub("chr", "",chr), "_MAF0.1_GRCh38_sorted.recode.positions", sep=""), sep ="\t", header = TRUE, stringsAsFactors=FALSE)
GRCh38_bp <- GRCh38_bp[which(GRCh38_bp$POS >= snp_bp - 5000000),]
GRCh38_bp <- GRCh38_bp[which(GRCh38_bp$POS <= snp_bp + 5000000),]
for(i in 1:length(GRCh38_bp$CHROM)){
  GRCh38_bp$kk[i] <- as.character(as.vector(unlist(strsplit(GRCh38_bp$ID[i], "_")))[1])
  GRCh38_bp$kk2[i] <- as.character(as.vector(unlist(strsplit(GRCh38_bp$ID[i], "_")))[2])
  GRCh38_bp$ID[i] <- paste(GRCh38_bp$kk[i], GRCh38_bp$kk2[i],sep="_")
}
GWAS1_GRCh38 <- merge(GWAS1_1, GRCh38_bp, by.x  = "RSid", by.y = "ID", all.x = TRUE)
GWAS1_GRCh38_2 <- subset(GWAS1_GRCh38,!GWAS1_GRCh38$POS %in% GWAS1_GRCh38$POS[is.na(GWAS1_GRCh38$POS)])
###
GWAS1_GRCh38_2$GWAS_lead <- ifelse(GWAS1_GRCh38_2$POS == snp_bp, 1 , 0)

##Add LD info
LD <- read.table(paste(snp, ".ld", sep=""), sep="\t", stringsAsFactors=FALSE)
names(LD) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
Add_autoLD <- c(chr, snp_bp, snp, chr, snp_bp, snp, "1", "NA")
LD <- rbind(LD, Add_autoLD)
LD2 <- LD[ which(LD$bp1==snp_bp),]
LD3 <- LD[ which(LD$bp2==snp_bp),]
LD3 <- LD3[,c(1,5,6,4,2,3,7,8)]
names(LD3) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
LD4 <- rbind(LD2,LD3)
LD4 <- LD4[ which(LD4$r2 >= 0.8),]
#LD4$color <- ifelse((LD4$bp2 >= as.numeric(k27acpeak_region[2])) & (LD4$bp2 <= as.numeric(k27acpeak_region[3])), "red", "black")
LD4 <- LD4[order(LD4$bp2),]
from_zoomed = min(as.numeric(LD4$bp2))-55000
to_zoomed = max(as.numeric(LD4$bp2))+55000

for(i in 1:length(LD4$snp2)){
  LD4$kk[i] <- as.character(as.vector(unlist(strsplit(LD4$snp2[i], "_")))[1])
  LD4$kk2[i] <- as.character(as.vector(unlist(strsplit(LD4$snp2[i], "_")))[2])
  LD4$snp2[i] <- paste(LD4$kk[i], LD4$kk2[i],sep="_")
}

missingLDblock2 <- LD4[-which(LD4$snp2 %in% GWAS1_GRCh38$RSid),]
missingLDblock2 <- missingLDblock2[, names(missingLDblock2) %in% c("snp2", "bp2")]
names(missingLDblock2) <- c("POS","RSid")
missingLDblock2$Chr <- rep(chr,length(missingLDblock2$RSid))
missingLDblock2$Pos <- rep(0,length(missingLDblock2$RSid))
missingLDblock2$Eff_allele <- rep(0,length(missingLDblock2$RSid))
missingLDblock2$MAF <- rep(0,length(missingLDblock2$RSid))
missingLDblock2$pval <- rep(0,length(missingLDblock2$RSid))
missingLDblock2$beta <- rep(0,length(missingLDblock2$RSid))
missingLDblock2$OR <- rep(0,length(missingLDblock2$RSid))
missingLDblock2$log_OR <- rep(0,length(missingLDblock2$RSid))
missingLDblock2$se <- rep(0,length(missingLDblock2$RSid))
missingLDblock2$z.score <- rep(0,length(missingLDblock2$RSid))
missingLDblock2$Disease <- rep(0,length(missingLDblock2$RSid))
missingLDblock2$PubmedID <- rep(0,length(missingLDblock2$RSid))
missingLDblock2$used_file <- rep(0,length(missingLDblock2$RSid))
missingLDblock2$GWAS_lead <- rep(0,length(missingLDblock2$RSid))
missingLDblock2$CHROM <- rep(0,length(missingLDblock2$RSid))
GWAS1_GRCh38_2 <- GWAS1_GRCh38_2[ , -which(names(GWAS1_GRCh38_2) %in% c("kk.x","kk2.x","kk.y", "kk2.y", "X", "NA."))]
missingLDblock2 <- missingLDblock2[,names(GWAS1_GRCh38_2)]
GWAS1_GRCh38_2 <- rbind(GWAS1_GRCh38_2, missingLDblock2)

LD <- read.table(paste(snp, ".ld", sep=""), sep="\t")
names(LD) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
Add_autoLD <- c(chr, snp_bp, snp, chr, snp_bp, snp, "1", "NA")
LD <- rbind(LD, Add_autoLD)
LD2 <- LD[ which(LD$bp1==snp_bp),]
LD3 <- LD[ which(LD$bp2==snp_bp),]
LD3 <- LD3[,c(1,5,6,4,2,3,7,8)]
names(LD3) <- c("chr_snp", "bp1", "snp1", "chr_snp2","bp2", "snp2", "r2", "null")
LD5 <- rbind(LD2,LD3)
rm(LD2,LD3)
LD2 <- merge(GWAS1_GRCh38,LD5, by.x = "POS", by.y = "bp2", all.x = TRUE) 
LD2 <- unique(LD2)
LD2$r2 <- as.numeric(LD2$r2)
LD2 <- LD2[which(LD2$POS >= from_zoomed),]
LD2 <- LD2[which(LD2$POS <= to_zoomed),]
LD2$GWAS_lead <- ifelse(LD2$POS == snp_bp, 1 ,0)
LD2_gwas <- LD2
LD2_gwas$r2[is.na(LD2_gwas$r2)] <- 0

if(k27acpeak != "chrNA"){
  pvalGWAS1_2 <- ggplot(LD2_gwas, aes(x=as.numeric(LD2_gwas$POS), y=as.numeric(LD2_gwas$pval), 
                              shape=as.factor(LD2_gwas$GWAS_lead)))+
    annotate("rect", xmin=as.numeric(k27acpeak_region[2]), xmax=as.numeric(k27acpeak_region[3]), ymin=0, ymax=Inf, alpha=0.5, fill="#F4D7C3")+
    geom_point(aes(fill = cut(as.numeric(LD2_gwas$r2), c(0, 0.2, 0.5, 0.8, 1), include.lowest = TRUE)),
               size = LD2_gwas$GWAS_lead+3, colour =  "#bebebe") +
    scale_shape_manual(values=c(21, 23))+
    scale_fill_manual(name = paste("r2 with ",snp, sep=""),
                      values = c("[0,0.2]" = "#f9f9f9",
                                 "(0.2,0.5]" = "#34495e",
                                 "(0.5,0.8]" = "#3498db",
                                 "(0.8,1]" = "#9b59b6"),
                      labels = c("<= 0.2", "0.2 < r2 <= 0.5", "0.5 < r2 <= 0.8","> 0.8"))+
    coord_cartesian(xlim = c(from_zoomed, to_zoomed))+
    labs(fill = paste("R2 with ", snp), sep = "") +
    xlab(paste(chr, " bp", sep = ""))+ ylab(paste("-LOG10(",GWAS, " p-value)", sep=""))+
    theme_classic()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          legend.position = "none")

#  ggsave(paste(gene_name, "/treg_",GWAS[1],"_",gene_name, "_", gene, "_", snp, "_","_Zoomed_pvalue.pdf", sep =""), plot = pvalGWAS1_2, width = 15, height = 5, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(k27acpeak == "chrNA" & k4me3peak != "chrNA"){
  pvalGWAS1_2 <- ggplot(LD2_gwas, aes(x=as.numeric(LD2_gwas$POS), y=as.numeric(LD2_gwas$pval),
                              shape=as.factor(LD2_gwas$GWAS_lead)))+
    annotate("rect", xmin=as.numeric(k4me3peak_region[2]), xmax=as.numeric(k4me3peak_region[3]), ymin=0, ymax=Inf, alpha=0.3, fill="#F4D7C3")+
    geom_point(aes(fill = cut(as.numeric(LD2_gwas$r2), c(0, 0.2, 0.5, 0.8, 1), include.lowest = TRUE)),
               size = LD2_gwas$GWAS_lead+3, colour =  "#bebebe") +
    scale_shape_manual(values=c(21, 23))+
    scale_fill_manual(name = paste("r2 with ",snp, sep=""),
                      values = c("[0,0.2]" = "#f9f9f9",
                                 "(0.2,0.5]" = "#34495e",
                                 "(0.5,0.8]" = "#3498db",
                                 "(0.8,1]" = "#9b59b6"),
                      labels = c("<= 0.2", "0.2 < r2 <= 0.5", "0.5 < r2 <= 0.8","> 0.8"))+
    coord_cartesian(xlim = c(from_zoomed, to_zoomed))+
    labs(fill = paste("R2 with ", snp), sep = "") +
    xlab(paste(chr, " bp", sep = ""))+ ylab(paste("-LOG10(",GWAS, " p-value)", sep=""))+
    theme_classic()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          legend.position = "none")

#  ggsave(paste(gene_name, "/treg_",GWAS[1],"_",gene_name, "_", gene, "_", snp, "_","_Zoomed_pvalue.pdf", sep =""), plot = pvalGWAS1_2, width = 15, height = 5, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(k27acpeak == "chrNA" & k4me3peak == "chrNA" & atacpeak != "chrNA"){
  pvalGWAS1_2 <- ggplot(LD2_gwas, aes(x=as.numeric(LD2_gwas$POS), y=as.numeric(LD2_gwas$pval),
                              shape=as.factor(LD2_gwas$GWAS_lead)))+
    annotate("rect", xmin=as.numeric(atacpeak_region[2]), xmax=as.numeric(atacpeak_region[3]), ymin=0, ymax=Inf, alpha=0.5, fill="#F4D7C3")+
    geom_point(aes(fill = cut(as.numeric(LD2_gwas$r2), c(0, 0.2, 0.5, 0.8, 1), include.lowest = TRUE)),
               size = LD2_gwas$GWAS_lead+3, colour =  "#bebebe") +
    scale_shape_manual(values=c(21, 23))+
    scale_fill_manual(name = paste("r2 with ",snp, sep=""),
                      values = c("[0,0.2]" = "#f9f9f9",
                                 "(0.2,0.5]" = "#34495e",
                                 "(0.5,0.8]" = "#3498db",
                                 "(0.8,1]" = "#9b59b6"),
                      labels = c("<= 0.2", "0.2 < r2 <= 0.5", "0.5 < r2 <= 0.8","> 0.8"))+
    coord_cartesian(xlim = c(from_zoomed, to_zoomed))+
    labs(fill = paste("R2 with ", snp), sep = "") +
    xlab(paste(chr, " bp", sep = ""))+ ylab(paste("-LOG10(",GWAS, " p-value)", sep=""))+
    theme_classic()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          legend.position = "none")

#  ggsave(paste(gene_name, "/treg_",GWAS[1],"_",gene_name, "_", gene, "_", snp, "_","_Zoomed_pvalue.pdf", sep =""), plot = pvalGWAS1_2, width = 15, height = 5, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(k27acpeak == "chrNA" & k4me3peak == "chrNA" & atacpeak == "chrNA"){
  pvalGWAS1_2 <- ggplot(LD2_gwas, aes(x=as.numeric(LD2_gwas$POS), y=as.numeric(LD2_gwas$pval),
                              shape=as.factor(LD2_gwas$GWAS_lead)))+
    geom_point(aes(fill = cut(as.numeric(LD2_gwas$r2), c(0, 0.2, 0.5, 0.8, 1), include.lowest = TRUE)),
               size = LD2_gwas$GWAS_lead+3, colour =  "#bebebe") +
    scale_shape_manual(values=c(21, 23))+
    scale_fill_manual(name = paste("r2 with ",snp, sep=""),
                      values = c("[0,0.2]" = "#f9f9f9",
                                 "(0.2,0.5]" = "#34495e",
                                 "(0.5,0.8]" = "#3498db",
                                 "(0.8,1]" = "#9b59b6"),
                      labels = c("<= 0.2", "0.2 < r2 <= 0.5", "0.5 < r2 <= 0.8","> 0.8"))+
    coord_cartesian(xlim = c(from_zoomed, to_zoomed))+
    labs(fill = paste("R2 with ", snp), sep = "") +
    xlab(paste(chr, " bp", sep = ""))+ ylab(paste("-LOG10(",GWAS, " p-value)", sep=""))+
    theme_classic()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          legend.position = "none")

#  ggsave(paste(gene_name, "/treg_",GWAS[1],"_",gene_name, "_", gene, "_", snp, "_","_Zoomed_pvalue.pdf", sep =""), plot = pvalGWAS1_2, width = 15, height = 5, units = "cm",
#         dpi = 300,device = "pdf", useDingbats=FALSE)
}



##combined zoomed Manhattan plots
if(gene != "Peak" & atacpeak == "chrNA" & k4me3peak == "chrNA" & k27acpeak == "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_3),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene != "Peak" & atacpeak != "chrNA" & k4me3peak != "chrNA" & k27acpeak != "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_3), arrangeGrob(pval1_2), arrangeGrob(pval1_4),arrangeGrob(pval1_5),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*5, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene != "Peak" & atacpeak == "chrNA" & k4me3peak == "chrNA" & k27acpeak != "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_3), arrangeGrob(pval1_2),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*3, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}


if(gene != "Peak" & atacpeak == "chrNA" & k4me3peak != "chrNA" & k27acpeak == "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_3), arrangeGrob(pval1_4),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*3, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene != "Peak" & atacpeak != "chrNA" & k4me3peak == "chrNA" & k27acpeak == "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_3), arrangeGrob(pval1_5),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*3, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene != "Peak" & atacpeak == "chrNA" & k4me3peak != "chrNA" & k27acpeak != "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_3), arrangeGrob(pval1_2), arrangeGrob(pval1_4),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*4, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene != "Peak" & atacpeak != "chrNA" & k4me3peak == "chrNA" & k27acpeak != "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_3), arrangeGrob(pval1_2),arrangeGrob(pval1_5),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*4, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene != "Peak" & atacpeak != "chrNA" & k4me3peak != "chrNA" & k27acpeak == "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_3), arrangeGrob(pval1_4),arrangeGrob(pval1_5),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*4, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}
####
if(gene == "Peak" & atacpeak != "chrNA" & k4me3peak != "chrNA" & k27acpeak != "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_2), arrangeGrob(pval1_4),arrangeGrob(pval1_5),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height =6*4, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene == "Peak" & atacpeak == "chrNA" & k4me3peak == "chrNA" & k27acpeak != "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_2),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}


if(gene == "Peak" & atacpeak == "chrNA" & k4me3peak != "chrNA" & k27acpeak == "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_4),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene == "Peak" & atacpeak != "chrNA" & k4me3peak == "chrNA" & k27acpeak == "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_5),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene == "Peak" & atacpeak == "chrNA" & k4me3peak != "chrNA" & k27acpeak != "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_2), arrangeGrob(pval1_4),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene == "Peak" & atacpeak != "chrNA" & k4me3peak == "chrNA" & k27acpeak != "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_2),arrangeGrob(pval1_5),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}

if(gene == "Peak" & atacpeak != "chrNA" & k4me3peak != "chrNA" & k27acpeak == "chrNA"){
  allpval <- arrangeGrob(arrangeGrob(pvalGWAS1_2), arrangeGrob(pval1_4),arrangeGrob(pval1_5),ncol=1)
  ggsave(paste(gene_name, "/treg_AllQTL_",gene_name, "_", gene, "_", snp, "_",window, "_Zoomed_pvalues.pdf", sep =""), plot = allpval, width = 15, height = 6*2, units = "cm",
           dpi = 300,device = "pdf", useDingbats=FALSE)
}
