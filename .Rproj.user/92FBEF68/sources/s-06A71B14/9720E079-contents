
rm(list = ls())

library(qsea)
library(BSgenome.Hsapiens.UCSC.hg19)
library(MEDIPSData)
library(GenomeInfoDb)
#library(biomaRt)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#library(BSgenome.Hsapiens.UCSC.hg38)

#library(BiocParallel)
#register(MulticoreParam(workers=4))

## Create sample table
Samples <- data.frame(sample_name = c("sgC", "sgPKCI"),
                      file_name = c("./data/bam2/sgC5mC/sgC5mC_Aligned.sortedByCoord.out.bam",
                                    "./data/bam2/sgPKCI5mC/sgPKCI5mC_Aligned.sortedByCoord.out.bam"),
                      group = c("control", "PKCIknocked"),
                      input_files = c("./data/bam2/sgCinput/sgCinput_Aligned.sortedByCoord.out.bam",
                                      "./data/bam2/sgPKCIinput/sgPKCIinput_Aligned.sortedByCoord.out.bam"))

## Create qsea set
qseaSet <- createQseaSet(sampleTable = Samples, 
                      BSgenome = "BSgenome.Hsapiens.UCSC.hg19", 
                      window_size = 400)
qseaSet


## compute the MeDIP coverage for each window
qseaSet <- addCoverage(qseaSet, uniquePos = TRUE, paired = FALSE, fragment_length = 80)

## Compute CNVs
qseaSet <- addCNV(qseaSet, file_name = "input_files", window_size=10000, 
               paired=FALSE, parallel=FALSE, MeDIP=FALSE, fragment_length = 80)


## Scaling Library Factor
qseaSet <- addLibraryFactors(qseaSet)

###Estimating model parameters for transformation to absolute methylation values
qseaSet <- addPatternDensity(qseaSet, "CG", name="CpG")

# From the regions without CpGs we can estimate the coverage offset from background reads.
qseaSet <- addOffset(qseaSet, enrichmentPattern = "CpG")


wd <- which(getRegions(qseaSet)$CpG_density>1 &
           getRegions(qseaSet)$CpG_density<15)

signal <- (15-getRegions(qseaSet)$CpG_density[wd])*.55/15+.25
qseaSet <- addEnrichmentParameters(qseaSet, enrichmentPattern="CpG", 
                                      windowIdx=wd, signal=signal)

getOffset(qseaSet, scale="fraction")

plotEPmatrix(qseaSet)

############################################

plotCNV(qseaSet)


## Differential Methylation Analysis
design <- model.matrix(~group, getSampleTable(qseaSet) )
qseaGLM <- fitNBglm(qseaSet, design, norm_method="beta")
qseaGLM <- addContrast(qseaSet, qseaGLM, coef = 2, name = "PKCIkOvsNorm")

##################
##Annotating, Exploring and Exporting Results

pca_cgi<- getPCA(qseaSet, norm_method="beta")

col <- rep(c("red", "green"), 3)

plotPCA(pca_cgi, bg = col, main="PCA plot based on CpG Island Promoters")


library(GenomicRanges)

## Keep only significant windows
sig <- isSignificant(qseaGLM, fdr_th=.01)

#Annotation <- GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg19))

## Make annotation object:

# Construct txdb from the genome gtf file
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb <- makeTxDbFromGFF("./data/ForIndex_hg19/hg19.refGene.gtf", format="gtf", organism = "Homo sapiens")

transcript <- transcripts(txdb)
exon <-  exons(txdb)
gene <-  genes(txdb)
Annotation <- list(gene = gene, transcript = transcript, exon = exon)

#################################
## Results all (FDR = 0.01)
result <- makeTable(qseaSet, 
                 glm=qseaGLM, 
                 groupMeans=getSampleGroups(qseaSet), 
                 keep=sig, 
                 annotation=Annotation, 
                 norm_method="beta")

knitr::kable(head(result))


###################################
## Results by gain vs loss (FDR = 0.1)
sigList <- list(gain=isSignificant(qseaGLM, fdr_th=.1,direction="gain"),
             loss=isSignificant(qseaGLM, fdr_th=.1,direction="loss"))

result_gain <- makeTable(qseaSet, 
                    glm=qseaGLM, 
                    groupMeans=getSampleGroups(qseaSet), 
                    keep=sigList$gain, 
                    annotation=Annotation, 
                    norm_method="beta")

result_loss <- makeTable(qseaSet, 
                         glm=qseaGLM, 
                         groupMeans=getSampleGroups(qseaSet), 
                         keep=sigList$loss, 
                         annotation=Annotation, 
                         norm_method="beta")

###################################
## ROIs stats > all
roi_stats <- regionStats(qseaSet, subsets=sigList, ROIs=Annotation)



ts_rel <- roi_stats[,-1]/roi_stats[,1]
x <- barplot(t(ts_rel)*100,ylab="fraction of ROIs[%]",
        names.arg=rownames(roi_stats), beside=TRUE, legend=TRUE, 
          las=2, args.legend=list(x="topleft"), 
          main="Feature enrichment PRKCI_KO vs WT")

#text(x=x[2,],y=-.15,labels=rownames(ts_rel), xpd=TRUE, srt=30, cex=1, adj=c(1,0))



########################################
## Look at genes in DMRs


AR <- result_gain[grep("CDKN", result_gain$gene), ]















