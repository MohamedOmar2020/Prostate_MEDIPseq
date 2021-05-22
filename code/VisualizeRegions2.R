rm(list = ls())

library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicFeatures)
library(annotatr)
library(bsseq)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

#########
## Change the plotting parameters
#getOption("Gviz.scheme")
#scheme <- getScheme()
#scheme$GeneRegionTrack$fill <- "salmon"
#scheme$GeneRegionTrack$col <- NULL
#scheme$GeneRegionTrack$transcriptAnnotation <- "symbol"
# scheme$GdObject$background.title <- "white"
# scheme$GdObject$cex.axis <- 12
# addScheme(scheme, "myScheme")
# options(Gviz.scheme = "myScheme")

##############
## Load the annotation list
annotations <- load("./objs/annListGR.rda")

############################
## Build annotation tracks
#txdb <- makeTxDbFromGFF("./data/ForIndex_hg19/hg19.refGene.gtf", format="gtf", organism = "Homo sapiens")

# Select annotations for intersection with regions
#annots <- c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
#            'hg19_genes_intronexonboundaries')


# Build the cpg track
cpgTrack <- build_annotations(genome = 'hg19', annotations = "hg19_cpgs")

cpgTrack <- AnnotationTrack(cpgTrack, name = "CpG",
                            fill = "black", col="black",
                            collapse=TRUE,
                            stacking= "squish")

# Build a genome axis track
# This indicate the genomic coordinates we are currently looking at in order to provide some reference.
gTrack <- GenomeAxisTrack(add53 = TRUE, add35 = TRUE, littleTicks = TRUE)


### Further annotation tracks

### Prepare annotations
tmp <- annList[c("enhancers fantom", "genes promoters", "CpG islands")]
names(tmp) <- c("Enhancers", "Promoters", "CpG island")
annCol <- c("lavender", "azure", "lightgreen")

trackAnnList <- mapply(x=tmp, nms=names(tmp), myCol=annCol,
                       FUN=function(x, nms, myCol, ...) {
                         AnnotationTrack(x, trackType = "AnnotationTrack",
                                         fill = myCol, genome = "hg19", name = nms,
                                         id=mcols(x)$name,
                                         featureAnnotation = "id",
                                         cex=0.8, fontcolor.feature = "dodgerblue4", fontface=2)
                       })

##########################################
## Load data and build data tracks

## First, the bigwig files, converted from bedgraph files generated during peak calling by macs3
sgC_5mC_BigWig <- import.bw("./MACS/bigwig/sgC_treat_pileup.bw")
#sgC_input_BigWig <- import.bw("./MACS/bigwig/sgC_control_lambda.bw")
sgPKCI_5mC_BigWig <- import.bw("./MACS/bigwig/sgPKCI_treat_pileup.bw")
#sgPKCI_input_BigWig <- import.bw("./MACS/bigwig/sgPKCI_control_lambda.bw")


sgC_5mC_Track <- DataTrack(range = sgC_5mC_BigWig, genome = "hg19", name = "sgC5mC", fill.histogram = "blue", col.histogram = "blue", cex.sampleNames = 12, type = "histogram")
#sgC_input_BigWig <- DataTrack(range = sgC_input_BigWig, genome = "hg19", name = "sgCinput")
sgPKCI_5mC_Track <- DataTrack(range = sgPKCI_5mC_BigWig, genome = "hg19", name = "sgPKCI5mC", fill.histogram = "red", col.histogram = "red", cex.sampleNames	= 12, type = "histogram")
#sgPKCI_input_BigWig <- DataTrack(range = sgPKCI_input_BigWig, genome = "hg19", name = "sgPKCIinput")


## Then the bed files resulting from the differntail methylation analysis using macs3 bdgdiff
MRs_Common_path <- "./MACS/Diff_cutoff_3/diff_sgC_vs_sgPRKCI_c3.0_common.bed"
MRs_sgC_path <- "./MACS/Diff_cutoff_3/diff_sgC_vs_sgPRKCI_c3.0_cond1.bed"
MRs_sgPKCI_path <- "./MACS/Diff_cutoff_3/diff_sgC_vs_sgPRKCI_c3.0_cond2.bed"

MRs_Common <- GRanges(import.bed(MRs_Common_path))
MRs_sgC <- GRanges(import.bed(MRs_sgC_path))
MRs_sgPKCI <- GRanges(import.bed(MRs_sgPKCI_path))

# Modify the names
MRs_sgPKCI$name <- gsub("diff_sgC_vs_sgPRKCI_cond2", "DMRs_sgPKCI_vs_sgC", MRs_sgPKCI$name)

#############
## Filter the regions based on the score
MRs_sgPKCI_ord <- MRs_sgPKCI[order(MRs_sgPKCI$score, decreasing = T), ]

# Keep the top 500
MRs_sgPKCI_Fil <- MRs_sgPKCI_ord[1:500, ] 

#MRs_sgPKCI_Fil <- subset(MRs_sgPKCI, score >= 50) # From 9535 to 3522

### Drop unused levels in seqnames
keep <- seqlevels(MRs_sgPKCI_Fil) %in% unique(seqnames(MRs_sgPKCI_Fil))
seqlevels(MRs_sgPKCI_Fil) <- seqlevels(MRs_sgPKCI_Fil)[keep]

# sgC_5mC_BigWig_fil <- sgC_5mC_BigWig[seqnames(sgC_5mC_BigWig) %in% seqnames(MRs_sgPKCI_Fil), ]
# seqlevels(sgC_5mC_BigWig_fil, pruning.mode="coarse") <- seqlevels(sgC_5mC_BigWig_fil)[keep]
# sgC_5mC_Track_fil <- DataTrack(range = sgC_5mC_BigWig_fil, genome = "hg19", name = "sgC5mC", fill.histogram = "blue", col.histogram = "blue", cex.sampleNames = 12, type = "histogram")
# 
# sgPKCI_5mC_BigWig_fil <- sgPKCI_5mC_BigWig[seqnames(sgPKCI_5mC_BigWig) %in% seqnames(MRs_sgPKCI_Fil), ]
# seqlevels(sgPKCI_5mC_BigWig_fil, pruning.mode="coarse") <- seqlevels(sgPKCI_5mC_BigWig_fil)[keep]
# sgPKCI_5mC_Track_fil <- DataTrack(range = sgPKCI_5mC_BigWig_fil, genome = "hg19", name = "sgCPKCI5mC", fill.histogram = "blue", col.histogram = "blue", cex.sampleNames = 12, type = "histogram")


###########################################################################
### Chormosomes with regions of interest
#allChr_MRs_Common <- unique(MRs_Common@seqnames@values)
#allChr_MRs_sgC <- unique(MRs_sgC@seqnames@values)
#allChr_MRs_sgPKCI <- unique(MRs_sgPKCI@seqnames@values)

allChr_MRs_sgPKCI_fil <- unique(MRs_sgPKCI_Fil@seqnames@values)
#allChr_MRs_sgPKCI_fil <- droplevels(allChr_MRs_sgPKCI_fil)

### Distance to offSet upstream and downstream
maxDist <- 1000

## Load data and build data tracks

## First, the bigwig files, converted from bedgraph files generated during peak calling by macs3
# sgC_5mC_BigWig <- import.bw("./MACS/bigwig/sgC_treat_pileup.bw")
# sgPKCI_5mC_BigWig <- import.bw("./MACS/bigwig/sgPKCI_treat_pileup.bw")
# 
# 
# sgC_5mC_BigWig_fil <- sgC_5mC_BigWig
# sgC_5mC_BigWig_fil <- sgC_5mC_BigWig_fil[seqnames(sgC_5mC_BigWig_fil)%in% seqnames(MRs_sgPKCI_Fil)]
# seqlevels(sgC_5mC_BigWig_fil) <- seqlevels(sgC_5mC_BigWig_fil)[seqlevels(sgC_5mC_BigWig_fil) %in% unique(seqnames(sgC_5mC_BigWig_fil)) ]
# 
# sgPKCI_5mC_BigWig_fil <- sgPKCI_5mC_BigWig
# sgPKCI_5mC_BigWig_fil <- sgPKCI_5mC_BigWig_fil[seqnames(sgPKCI_5mC_BigWig_fil)%in% seqnames(MRs_sgPKCI_Fil)]
# seqlevels(sgPKCI_5mC_BigWig_fil) <- seqlevels(sgPKCI_5mC_BigWig_fil)[seqlevels(sgPKCI_5mC_BigWig_fil) %in% unique(seqnames(sgPKCI_5mC_BigWig_fil)) ]
# 
# 
# sgC_5mC_Track <- DataTrack(range = sgC_5mC_BigWig_fil, genome = "hg19", name = "sgC5mC", fill.histogram = "blue", col.histogram = "blue", cex.sampleNames = 12, type = "histogram")
# #sgC_input_BigWig <- DataTrack(range = sgC_input_BigWig, genome = "hg19", name = "sgCinput")
# sgPKCI_5mC_Track <- DataTrack(range = sgPKCI_5mC_BigWig_fil, genome = "hg19", name = "sgPKCI5mC", fill.histogram = "red", col.histogram = "red", cex.sampleNames	= 12, type = "histogram")
# #sgPKCI_input_BigWig <- DataTrack(range = sgPKCI_input_BigWig, genome = "hg19", name = "sgPKCIinput")

###########################################################################
### Open device
pdf("figs/sgPKCI_Top500DMRs2.pdf", width=11, height=8.5)

###########################################################################
### Start chromosome loop
for ( j in seq_along(allChr_MRs_sgPKCI_fil) ) {
  
  ## Set variable for chromosome of interest
  myChr <- as.character(allChr_MRs_sgPKCI_fil[j])
  ## Create an Ideogram
  ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = myChr)
  ## Get the Gene models for the chromosome
  genesTrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                name="Transcripts", 
                                chromosome=myChr,
                                fill="lightpink1", col.line="darkgreen", lwd=1.5)
  ## Add Gene symbols
  myIds <- gene(genesTrack)
  myIds <- select(org.Hs.eg.db, keys=myIds, keytype="ENTREZID", columns = "SYMBOL")
  myIds <- apply(myIds, 1, function(x) paste(x[2], x[1], sep="\n") )
  symbol(genesTrack) <- myIds
  
  ## Set the chromosomes for the data tracks
  ### chromosome(trackCov) <- myChr
  chromosome(sgC_5mC_Track) <- myChr
  chromosome(sgPKCI_5mC_Track) <- myChr
  
  ## Set Chr for annotations
  trackAnnList <- lapply(trackAnnList, function(x, myChr) {
    chromosome(x) <- myChr
    x
  }, myChr=myChr)
  
  ###########################################################################
  ## Get the regions
  sgPKCI_DMRsByChr <- MRs_sgPKCI_Fil[seqnames(MRs_sgPKCI_Fil) == myChr]
  
  ### Drop unused levels in seqnames
  #keep <- seqlevels(sgPKCI_DMRsByChr) %in% unique(seqnames(sgPKCI_DMRsByChr))
  #seqlevels(sgPKCI_DMRsByChr) <- seqlevels(sgPKCI_DMRsByChr)[keep]
  
  ## Start gene loop
  for (i in seq_along(sgPKCI_DMRsByChr)) {
     ## Define the genomic window to be plotted
     wStart <- start(sgPKCI_DMRsByChr)[i] - maxDist
     wEnd <-  end(sgPKCI_DMRsByChr)[i] + maxDist
     
    
     ## Highlight track
     ht <- HighlightTrack(trackList = c(
       trackAnnList, 
       list(genesTrack, sgC_5mC_Track, sgPKCI_5mC_Track)),
                         start = start(sgPKCI_DMRsByChr)[i], 
                         end = end(sgPKCI_DMRsByChr)[i],
                         chromosome = myChr,
                         inBackground=TRUE, 
                         genome="hg19")
     
     ## Plot
     plotTracks(list(ideoTrack, gTrack, ht),
               add=FALSE, from=wStart, to=wEnd, chromosome=myChr, grid=TRUE,
               transcriptAnnotation = "symbol",  #shape = "fixedArrow",
               stacking = "squish", ## jitter.x=TRUE, amount=10, #factor=10,
               background.panel = rgb(0.99, 0.99, 0.01, 0.1),
               background.title = "darkblue",
               main=paste(MRs_sgPKCI_Fil$name[i], "", "Genomic Region")
                )
    
  }
}

### close device
dev.off()






