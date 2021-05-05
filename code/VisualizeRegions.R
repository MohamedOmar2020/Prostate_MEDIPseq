
rm(list = ls())

library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicFeatures)

#########
## Change the plotting parameters
getOption("Gviz.scheme")
scheme <- getScheme()
scheme$GeneRegionTrack$fill <- "salmon"
scheme$GeneRegionTrack$col <- NULL
scheme$GeneRegionTrack$transcriptAnnotation <- "symbol"
addScheme(scheme, "myScheme")
options(Gviz.scheme = "myScheme")
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
grtrack <- GeneRegionTrack(txdb, genome = "hg19", name = "Gene Model")

# Build sequence track
strack <- SequenceTrack(Hsapiens)

##########################################
## Load data and build data tracks

## First, the bigwig files, converted from bedgraph files generated during peak calling by macs3
sgC_5mC_BigWig <- import.bw("./MACS/bigwig/sgC_treat_pileup.bw")
sgC_input_BigWig <- import.bw("./MACS/bigwig/sgC_control_lambda.bw")
sgPKCI_5mC_BigWig <- import.bw("./MACS/bigwig/sgPKCI_treat_pileup.bw")
sgPKCI_input_BigWig <- import.bw("./MACS/bigwig/sgPKCI_control_lambda.bw")


sgC_5mC_BigWig <- DataTrack(range = sgC_5mC_BigWig, genome = "hg19", name = "sgC5mC")
sgC_input_BigWig <- DataTrack(range = sgC_input_BigWig, genome = "hg19", name = "sgCinput")
sgPKCI_5mC_BigWig <- DataTrack(range = sgPKCI_5mC_BigWig, genome = "hg19", name = "sgPKCI5mC")
sgPKCI_input_BigWig <- DataTrack(range = sgPKCI_input_BigWig, genome = "hg19", name = "sgPKCIinput")


## Then the bed files resulting from the differntail methylation analysis using macs3 bdgdiff
MRs_Common <- "./MACS/Diff_cutoff_3/diff_sgC_vs_sgPRKCI_c3.0_common.bed"
MRs_sgC <- "./MACS/Diff_cutoff_3/diff_sgC_vs_sgPRKCI_c3.0_cond1.bed"
MRs_sgPKCI <- "./MACS/Diff_cutoff_3/diff_sgC_vs_sgPRKCI_c3.0_cond2.bed"

MRs_Common <- DataTrack(range = MRs_Common, genome = "hg19", name = "MRs_Common")
MRs_sgC <- DataTrack(range = MRs_sgC, genome = "hg19", name = "MRs_sgC")
MRs_sgPKCI <- DataTrack(range = MRs_sgPKCI, genome = "hg19", name = "MRs_sgPKCI")


##########################################################
## Plot
plotTracks(list(grtrack, gtrack, sgC_5mC_BigWig, sgPKCI_5mC_BigWig), type = "histogram",
           showSampleNames = TRUE, cex.sampleNames = 0.6,
           separator = 1, from = 49369522, to = 49378986,
           #sizes = c(25, 25, 100, 100),
           chromosome = "chr12",
           geneSymbols = T,
           transcriptAnnotation = "gene",
           #family = "gaussian",
           #evaluation = 50,
           col.histogram	= c("blue", "red"),
           fill.histogram	= c("blue", "red"),
           ylim = c(0, 100)
)








