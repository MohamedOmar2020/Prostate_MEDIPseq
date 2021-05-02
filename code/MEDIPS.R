
rm(list = ls())

library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
library(AnnotationDbi)

## The bam files
bam_sgC5mC <- "./data/bam/sgC5mC/sgC5mC_Aligned.sortedByCoord.out.bam"

bam_sgCinput <- "./data/bam/sgCinput/sgCinput_Aligned.sortedByCoord.out.bam"

bam_sgPKCI5mC <- "./data/bam/sgPKCI5mC/sgPKCI5mC_Aligned.sortedByCoord.out.bam"

bam_sgPKCIinput <- "./data/bam/sgPKCIinput/sgPKCIinput_Aligned.sortedByCoord.out.bam"

## The human genome
BSgenome="BSgenome.Hsapiens.UCSC.hg38"


#######################################
## Create MEDIPS sets

sgC5mC <- MEDIPS.createSet(file = bam_sgC5mC,
                               BSgenome = BSgenome, extend = 300, shift = 0, uniq = 1e-3,
                               window_size = 300)

sgCinput <- MEDIPS.createSet(file = bam_sgCinput,
                                  BSgenome = BSgenome, extend = 300, shift = 0, uniq = 1e-3,
                                  window_size = 300)


sgPKCI5mC <- MEDIPS.createSet(file = bam_sgPKCI5mC,
                                   BSgenome = BSgenome, extend = 300, shift = 0, uniq = 1e-3,
                                   window_size = 300)

sgPKCIinput <- MEDIPS.createSet(file = bam_sgPKCIinput,
                                   BSgenome = BSgenome, extend = 300, shift = 0, uniq = 1e-3,
                                   window_size = 300)

## Generate coupling set
CS <- MEDIPS.couplingVector(pattern = "CG", refObj = sgC5mC)

###########################################################################
## Calculate differential methylation

mr.edgeR <- MEDIPS.meth(MSet1 = sgPKCI5mC, MSet2 = sgC5mC,
                       CSet = CS, ISet1 = sgPKCIinput, ISet2 = sgCinput, p.adj = "bonferroni",
                       diff.method = "edgeR", MeDIP = T, CNV = F, minRowSum = 10)

###############
## selecting significant windows
mr.edgeR.sig <- MEDIPS.selectSig(results = mr.edgeR, p.value = 0.1,
                              adj = T, ratio = NULL, bg.counts = NULL, CNV = F)


##############
## Merging neighboring significant windows 

# extract significant windows in PRKCI versus control 
mr.edgeR.s.gain <- mr.edgeR.sig[which(mr.edgeR.sig[, grep("logFC",
                                                     colnames(mr.edgeR.sig))] > 0), ]

# merge all adjacent significant regions into one region which will be finally regarded as one event of differential coverage:

mr.edgeR.s.gain.merged <- MEDIPS.mergeFrames(frames = mr.edgeR.s.gain,
                                       distance = 1)


#################################################################
#################################################################
## Extracting data at regions of interest: 

# extract only the columns that contain count values
columns <- names(mr.edgeR)[grep("counts", names(mr.edgeR))]

# select all genomic windows in the original result table (here mr.edgeR) which are included in the boundaries defined by the genomic coordinates of the given set of ROIs
ROIs <- MEDIPS.selectROIs(results = mr.edgeR, rois = mr.edgeR.s.gain.merged,
                           columns = columns, summarize = NULL)

# As an alternative, it is also possible to calculate mean values over the extracted windows for each ROI, a behaviour that is controlled by the parameter summarize:
rois.summarized <- MEDIPS.selectROIs(results = mr.edgeR, rois = mr.edgeR.s.gain.merged,
                               columns = columns, summarize = "avg")


##################################################################
##################################################################
## Quality controls

# Saturation analysis ???? 
si.bsg <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
si.bam <- seqinfo(BamFile(bam_sgC5mC))

sr_sgC5mC <- MEDIPS.saturation(file = bam_sgC5mC, BSgenome = BSgenome,
                       uniq = 1, extend = 300, shift = 0, window_size = 300,
                       nit = 10, nrit = 1, empty_bins = TRUE,
                       rank = FALSE)
sr_sgC5mC

#########
## Correlation between samples
cor.matrix = MEDIPS.correlation(MSets = c(sgC5mC, sgPKCI5mC,
                                          sgCinput, sgPKCIinput), plot = T, method = "pearson")

############
## Sequence Pattern Coverage

cr = MEDIPS.seqCoverage(file = bam_sgC5mC, pattern = "CG",
                        BSgenome = BSgenome, extend = 300,
                        shift = 0, uniq = 0.001)



############
# CpG enrichment
er = MEDIPS.CpGenrich(file = bam_sgC5mC, BSgenome = BSgenome,
                      extend = 300, shift = 0,
                      uniq = 0.1)


###################################################################
####################################################################
## Annotation of significant windows
anno.mart.gene = MEDIPS.getAnnotation(dataset = c("hsapiens_gene_ensembl"),
                                      annotation = c("GENE"))
anno.mart.gene <- anno.mart.gene$Gene
mr.edgeR.s.gain.merged = MEDIPS.setAnnotation(regions = mr.edgeR.s.gain.merged, annotation = anno.mart.gene)


Up_Symbols <- mapIds(org.Hs.eg.db,
                              keys=mr.edgeR.s.gain.merged[ ,'1_id'],
                              column="SYMBOL",
                              keytype="ENSEMBL",
                              multiVals="first")

mr.edgeR.s.gain.merged$symbol <- Up_Symbols































