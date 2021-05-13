###########################################################################
### Luigi Marchionni
### Bisulfite Sequencing analysis to find DMRs
### Metastasis vs Primary baldder cancer, matched samples from Indiana
### Collaboration with Noah Hahn

###########################################################################
### Set wd and clean ws
setwd("~/Research/Bladder/labHahn/DMRs/")
rm(list=ls())

### Libraries
require(bsseq)
require(Gviz)
require(doParallel)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(org.Hs.eg.db)


###########################################################################
### Load separate sets
load("objs/bsObjsForPlots.rda")
load("objs/dmrsDatLarge.rda")
load("~/Research/Annotation/fromAnnotatr/objs/annListGR.rda")


###########################################################################
### Finding DMRs

### Find number of min methylation locations
nMin <- 10

### Subset
dmrsFil <- subset(dmrs, n >= nMin & abs(meanDiff) >= 0.33)
dim(dmrsFil)

### Remove the NAs Symbol
dmrsFil[is.na(dmrsFil[,"name"]) , "name"] <- "Closest gene not annotated in NCBI Gene"

### Regions
myDMRs <- GRanges(dmrsFil[, c("chr", "start", "end", "name")])

### Sort regions
myDMRs <- sort(myDMRs)


###########################################################################
### Prepare objects

### Genome axis
gTrack <- GenomeAxisTrack(add53 = TRUE, add35 = TRUE, littleTicks = TRUE)


###########################################################################
### Prepare data tracks

### Create groups variable
DiseaseStatus <- bsCovDmrs$Group

### Drop unused levels in seqnames
keep <- seqlevels(bsCovDmrs) %in% unique(seqnames(bsCovDmrs))
seqlevels(bsCovDmrs) <- seqlevels(bsCovDmrs)[keep]

### ### Create GenomicRanges for coverage
### covGR <- GRanges(rowRanges(bsCovDmrs))
### mcols(covGR) <- getCoverage(bsCovDmrs)
###
### ### Create tracks
### trackCov <- DataTrack(covGR, name = "RRBS coverage", groups=DiseaseStatus,
###                        genome="hg19", type = c("a", "p", "confint"))

### Create GenomicRanges for methylation
methGR <- GRanges(rowRanges(bsCovDmrs))
mcols(methGR) <- getMeth(bsCovDmrs)

### Create tracks
trackMeth <- DataTrack(methGR, name = "RRBS Methylation", groups=DiseaseStatus,
                       genome="hg19", type = c("a", "p", "confint"),
                       pch=16, cex=0.5)

### CpG locations
trackCpG <- AnnotationTrack(reduce(rowRanges(bsCovDmrs)),
                            fill = "black", col="black",
                            collapse=TRUE,
                            genome = "hg19", name = "CpG", stacking= "squish")

### ### CpG locations
### CpG <- reduce(rowRanges(bsCovDmrs)) ; mcols(CpG) <- 1
### trackCpG <- DataTrack(CpG, fill = "black", genome = "hg19", name = "CpG",
###                       cex=0.1, type="histo", ylim=c(0,0.1))


###########################################################################
### Annotation tracks

### Prepare annotations
tmp <- annList[c("enhancers fantom", "genes promoters", "CpG islands")]
names(tmp) <- c("Enhancers", "Promoters", "CpG island")
annCol <- c("lavender","azure", "lightgreen")

trackAnnList <- mapply(x=tmp, nms=names(tmp), myCol=annCol,
                       FUN=function(x, nms, myCol, ...) {
    AnnotationTrack(x, trackType = "AnnotationTrack",
                    fill = myCol, genome = "hg19", name = nms,
                    id=mcols(x)$name,
                    featureAnnotation = "id",
                    cex=0.8, fontcolor.feature = "dodgerblue4", fontface=2)
})


###########################################################################
### Chormosomes with regions of interest
allChr <- unique(dmrsFil$chr)

### Distance to offSet upstream and downstream
maxDist <- 1000


###########################################################################
### Open device
pdf("figs/METSvsPRIMARY_RRBS_DMRs.pdf", width=11, height=8.5)

###########################################################################
### Start chromosome loop
for ( j in seq_along(allChr) ) {

    ## Set variable for chromosome of interest
    myChr <- allChr[j]
    ## Create an Ideogram
    ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = myChr)
    ## Get the Gene models for the chromosome
    genesTrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                  name="Transcripts", chromosome=myChr,
                                  fill="lightpink1", col.line="darkgreen", lwd=1.5)
    ## Add Gene symbols
    myIds <- gene(genesTrack)
    myIds <- select(org.Hs.eg.db, keys=myIds, keytype="ENTREZID",columns = "SYMBOL")
    myIds <- apply(myIds, 1, function(x) paste(x[2], x[1], sep="\n") )
    symbol(genesTrack) <- myIds

    ## Set the chromosomes for the data tracks
    ### chromosome(trackCov) <- myChr
    chromosome(trackMeth) <- myChr
    chromosome(trackCpG) <- myChr

    ## Set Chr for annotations
    trackAnnList <- lapply(trackAnnList, function(x, myChr) {
        chromosome(x) <- myChr
        x
    }, myChr=myChr)

###########################################################################
    ## Get the regions
    myDMRsByChr <- myDMRs[seqnames(myDMRs) == myChr]

    ## Start gene loop
    for (i in seq_along(myDMRsByChr)) {
        ## Define the genomic window to be plotted
        wStart <- start(myDMRsByChr)[i] - maxDist
        wEnd <-  end(myDMRsByChr)[i] + maxDist

        ## Highlight track
        ht <- HighlightTrack(trackList = c(trackAnnList, list(genesTrack, trackMeth, trackCpG)),
                             start = start(myDMRsByChr)[i], end = end(myDMRsByChr)[i],
                             chromosome = myChr, inBackground=TRUE, genome="hg19")

        ## Plot
        plotTracks(list(ideoTrack, gTrack, ht),
                   add=FALSE, from=wStart, to=wEnd, chromosome=myChr, grid=TRUE,
                   transcriptAnnotation = "symbol",  #shape = "fixedArrow",
                   stacking = "squish", ## jitter.x=TRUE, amount=10, #factor=10,
                   background.panel = rgb(0.99, 0.99, 0.01, 0.1),
                   background.title = "darkblue",
                   main=paste(myDMRsByChr$name[i], "", "Genomic Region - RRBS in bladder cancer")
                   )

    }
}

### close device
dev.off()


### ###########################################################################
### ### Plotting using bsseq
###
### ### Open device
### pdf("figs/topDMRsLarge.pdf", width=10, height=10)
### ### Plot
### system.time(invisible(sapply(1:nrow(dmrsFil[1:100,]), function(x, bs, dmr, ...) {
###         plotRegion(bs, dmrsFil[x,], extend=1000, addRegions = dmr,
###                    main=dmrsFil[x,"name"], ...)
### }, bs=bsCovDmrs, dmr=dmrsFil, annoTrack=annList)))
### ### Close device
### dev.off()


###########################################################################
### Session and quit

### Session
sessionInfo()


### Quite
q("no")
