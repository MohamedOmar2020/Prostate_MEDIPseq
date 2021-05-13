
rm(list = ls())

library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicFeatures)
library(annotatr)

#########
## Change the plotting parameters
getOption("Gviz.scheme")
#scheme <- getScheme()
#scheme$GeneRegionTrack$fill <- "salmon"
#scheme$GeneRegionTrack$col <- NULL
#scheme$GeneRegionTrack$transcriptAnnotation <- "symbol"
# scheme$GdObject$background.title <- "white"
# scheme$GdObject$cex.axis <- 12
# addScheme(scheme, "myScheme")
# options(Gviz.scheme = "myScheme")


############################
## Build annotation tracks

txdb <- makeTxDbFromGFF("./data/ForIndex_hg19/hg19.refGene.gtf", format="gtf", organism = "Homo sapiens")

# Select annotations for intersection with regions
annots <- c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
            'hg19_genes_intronexonboundaries')

# Build the annotations (a single GRanges object)
cpgTrack <- build_annotations(genome = 'hg19', annotations = "hg19_cpgs")
cpgTrack <- AnnotationTrack(cpgTrack, name = "CpG")

# Build a genome axis track
# This indicate the genomic coordinates we are currently looking at in order to provide some reference.
gtrack <- GenomeAxisTrack()

# Build chromosome ideogram track
itrack <- IdeogramTrack(genome = "hg19")

# Build genome region track
grtrack <- GeneRegionTrack(txdb, genome = "hg19")

# Build sequence track
strack <- SequenceTrack(Hsapiens)

##########################################
## Load data and build data tracks

## First, the bigwig files, converted from bedgraph files generated during peak calling by macs3
sgC_5mC_BigWig <- import.bw("./MACS/bigwig/sgC_treat_pileup.bw")
#sgC_input_BigWig <- import.bw("./MACS/bigwig/sgC_control_lambda.bw")
sgPKCI_5mC_BigWig <- import.bw("./MACS/bigwig/sgPKCI_treat_pileup.bw")
#sgPKCI_input_BigWig <- import.bw("./MACS/bigwig/sgPKCI_control_lambda.bw")


sgC_5mC_BigWig <- DataTrack(range = sgC_5mC_BigWig, genome = "hg19", name = "sgC5mC", fill.histogram = "blue", col.histogram = "blue", cex.sampleNames = 12, type = "histogram")
#sgC_input_BigWig <- DataTrack(range = sgC_input_BigWig, genome = "hg19", name = "sgCinput")
sgPKCI_5mC_BigWig <- DataTrack(range = sgPKCI_5mC_BigWig, genome = "hg19", name = "sgPKCI5mC", fill.histogram = "red", col.histogram = "red", cex.sampleNames	= 12, type = "histogram")
#sgPKCI_input_BigWig <- DataTrack(range = sgPKCI_input_BigWig, genome = "hg19", name = "sgPKCIinput")

#displayPars(sgPKCI_5mC_BigWig)
#displayPars(sgPKCI_5mC_BigWig)$fill.histogram <- "red"

## Then the bed files resulting from the differntail methylation analysis using macs3 bdgdiff
MRs_Common <- "./MACS/Diff_cutoff_3/diff_sgC_vs_sgPRKCI_c3.0_common.bed"
MRs_sgC <- "./MACS/Diff_cutoff_3/diff_sgC_vs_sgPRKCI_c3.0_cond1.bed"
MRs_sgPKCI <- "./MACS/Diff_cutoff_3/diff_sgC_vs_sgPRKCI_c3.0_cond2.bed"

MRs_Common <- DataTrack(range = MRs_Common, genome = "hg19", name = "MRs_Common", size	= 3, col.histogram = "black", fill.histogram = "black", type = "histogram")
MRs_sgC <- DataTrack(range = MRs_sgC, genome = "hg19", name = "MRs_sgC", size	= 3, col.histogram = "blue", fill.histogram = "blue", type = "histogram")
MRs_sgPKCI <- DataTrack(range = MRs_sgPKCI, genome = "hg19", name = "MRs_sgPKCI", size	= 3, col.histogram = "red", fill.histogram = "red", type = "histogram")


##########################################################
## Plotting

## Wnt1
png(filename = "./figs/wnt1.png", width = 2500, height = 2000, res = 300)
plotTracks(list(itrack, gtrack, grtrack, sgC_5mC_BigWig, sgPKCI_5mC_BigWig, 
                #MRs_Common, MRs_sgC, 
                MRs_sgPKCI), 
           #type = "histogram",
           showSampleNames = TRUE, 
           separator = 1, from = 49369522, to = 49378986,
           #sizes = c(25, 25, 100, 100),
           chromosome = "chr12",
           geneSymbols = T,
           transcriptAnnotation = "gene",
           #family = "gaussian",
           #evaluation = 50,
           #col.histogram	= c("blue", "red"),
           #fill.histogram	= c("blue", "red"),
           ylim = c(0, 100),
           background.title = "white",
           fontsize = 12,
           col.title = "black"
           )
dev.off()

##########
## CDKN1A
png(filename = "./figs/CDKN1A.png", width = 2500, height = 2000, res = 300)
plotTracks(list(itrack, gtrack, grtrack, sgC_5mC_BigWig, sgPKCI_5mC_BigWig
                #, MRs_sgPKCI
                ), 
           #type = "histogram",
           showSampleNames = TRUE, 
           separator = 1, from = 36643047, to = 36655146,
           #sizes = c(25, 25, 100, 100),
           chromosome = "chr6",
           geneSymbols = T,
           #collapseTranscripts = T,
           transcriptAnnotation = "gene",
           #family = "gaussian",
           #evaluation = 50,
           #col.histogram	= c("blue", "red"),
           #fill.histogram	= c("blue", "red"),
           ylim = c(0, 120),
           background.title = "white",
           fontsize = 12,
           col.title = "black"
)
dev.off()

###########
## ADAMTS1
png(filename = "./figs/ADAMTS1.png", width = 2500, height = 2000, res = 300)
plotTracks(list(itrack, gtrack, grtrack, sgC_5mC_BigWig, sgPKCI_5mC_BigWig
                #, MRs_sgPKCI
                ), 
            #type = "histogram",
            showSampleNames = TRUE, 
            separator = 1, from = 28208658, to = 28222312,
            #sizes = c(25, 25, 100, 100),
            chromosome = "chr21",
            geneSymbols = T,
            #collapseTranscripts = T,
            transcriptAnnotation = "gene",
            #family = "gaussian",
            #evaluation = 50,
            #col.histogram	= c("blue", "red"),
            #fill.histogram	= c("blue", "red"),
            ylim = c(0, 100),
            background.title = "white",
            fontsize = 12,
            col.title = "black"
            )
dev.off()

###########
## TMPRSS2
png(filename = "./figs/TMPRSS2.png", width = 2500, height = 2000, res = 300)
plotTracks(list(itrack, gtrack, grtrack, sgC_5mC_BigWig, sgPKCI_5mC_BigWig
                #, MRs_sgPKCI
), 
#type = "histogram",
showSampleNames = TRUE, 
separator = 1, from = 42832634, to = 42887088,
#sizes = c(25, 25, 100, 100),
chromosome = "chr21",
geneSymbols = T,
#collapseTranscripts = T,
transcriptAnnotation = "gene",
#family = "gaussian",
#evaluation = 50,
#col.histogram	= c("blue", "red"),
#fill.histogram	= c("blue", "red"),
ylim = c(0, 100),
background.title = "white",
fontsize = 12,
col.title = "black"
)
dev.off() 

#############
## NOTCH1
png(filename = "./figs/NOTCH1.png", width = 2500, height = 2000, res = 300)
plotTracks(list(itrack, gtrack, grtrack, sgC_5mC_BigWig, sgPKCI_5mC_BigWig
                #, MRs_sgPKCI
), 
#type = "histogram",
showSampleNames = TRUE, 
separator = 1, from = 139383045, to = 139438660,
#sizes = c(25, 25, 100, 100),
chromosome = "chr9",
geneSymbols = T,
#collapseTranscripts = T,
transcriptAnnotation = "gene",
#family = "gaussian",
#evaluation = 50,
#col.histogram	= c("blue", "red"),
#fill.histogram	= c("blue", "red"),
ylim = c(0, 150),
background.title = "white",
fontsize = 12,
col.title = "black"
)
dev.off() 

##############
## CD44
png(filename = "./figs/CD44.png", width = 2500, height = 2000, res = 300)
plotTracks(list(itrack, gtrack, grtrack, sgC_5mC_BigWig, sgPKCI_5mC_BigWig
                #, MRs_sgPKCI
), 
#type = "histogram",
showSampleNames = TRUE, 
separator = 1, from = 35145980, to = 35222353,
#sizes = c(25, 25, 100, 100),
chromosome = "chr11",
geneSymbols = T,
#collapseTranscripts = T,
transcriptAnnotation = "gene",
#family = "gaussian",
#evaluation = 50,
#col.histogram	= c("blue", "red"),
#fill.histogram	= c("blue", "red"),
ylim = c(0, 100),
background.title = "white",
fontsize = 12,
col.title = "black"
)
dev.off() 

#################
## MTOR
png(filename = "./figs/MTOR.png", width = 2500, height = 2000, res = 300)
plotTracks(list(itrack, gtrack, grtrack, sgC_5mC_BigWig, sgPKCI_5mC_BigWig
                #, MRs_sgPKCI
), 
#type = "histogram",
showSampleNames = TRUE, 
separator = 1, from = 11258279, to = 11328780,
#sizes = c(25, 25, 100, 100),
chromosome = "chr1",
geneSymbols = T,
grid = T,
#collapseTranscripts = T,
transcriptAnnotation = "gene",
#family = "gaussian",
#evaluation = 50,
#col.histogram	= c("blue", "red"),
#fill.histogram	= c("blue", "red"),
ylim = c(0, 80),
background.title = "white",
fontsize = 12,
col.title = "black"
)
dev.off() 


#################
## ADCY5
png(filename = "./figs/ADCY5.png", width = 2500, height = 2000, res = 300)
plotTracks(list(itrack, gtrack, grtrack, sgC_5mC_BigWig, sgPKCI_5mC_BigWig
                #, MRs_sgPKCI
), 
#type = "histogram",
showSampleNames = TRUE, 
separator = 1, from = 122945725, to = 123167496,
#sizes = c(25, 25, 100, 100),
chromosome = "chr3",
geneSymbols = T,
grid = T,
#collapseTranscripts = T,
transcriptAnnotation = "gene",
#family = "gaussian",
#evaluation = 50,
#col.histogram	= c("blue", "red"),
#fill.histogram	= c("blue", "red"),
ylim = c(0, 120),
background.title = "white",
fontsize = 12,
col.title = "black"
)
dev.off() 

#################
## ROBO2
png(filename = "./figs/ADCY5.png", width = 2500, height = 2000, res = 300)
plotTracks(list(itrack, gtrack, grtrack, sgC_5mC_BigWig, sgPKCI_5mC_BigWig
                #, MRs_sgPKCI
), 
#type = "histogram",
showSampleNames = TRUE, 
separator = 1, from = 122945725, to = 123167496,
#sizes = c(25, 25, 100, 100),
chromosome = "chr3",
geneSymbols = T,
grid = T,
#collapseTranscripts = T,
transcriptAnnotation = "gene",
#family = "gaussian",
#evaluation = 50,
#col.histogram	= c("blue", "red"),
#fill.histogram	= c("blue", "red"),
ylim = c(0, 120),
background.title = "white",
fontsize = 12,
col.title = "black"
)
dev.off() 





