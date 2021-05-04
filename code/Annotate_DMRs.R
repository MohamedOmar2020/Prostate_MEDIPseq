

library(annotatr)
library(GenomicRanges)
library(regioneR)

########################################
## Load DMRs

# Unique to sgPKCI
DMRs_sgPKCI_vs_sgC <- read_regions(con = "./MACS/Diff/diff_sgC_vs_sgPRKCI_c1.0_cond2.bed", genome = 'hg19', format = 'bed')
DMRs_sgPKCI_vs_sgC

DMRs_sgPKCI_vs_sgC <- GRanges(DMRs_sgPKCI_vs_sgC)

# Unique to sgC
DMRs_sgC_vs_sgPKCI <- read_regions(con = "./MACS/Diff/diff_sgC_vs_sgPRKCI_c1.0_cond1.bed", genome = 'hg19', format = 'bed')
DMRs_sgC_vs_sgPKCI

DMRs_sgC_vs_sgPKCI <- GRanges(DMRs_sgC_vs_sgPKCI)

########################################
## Annotating Regions

# see available annotations
builtin_annotations()

# Select annotations for intersection with regions
# Note inclusion of custom annotation, and use of shortcuts
annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
           'hg19_genes_intronexonboundaries')

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg19', annotations = annots)

# Intersect the regions we read in with the annotations
DMRs_sgPKCI_vs_sgC_annotated = annotate_regions(
  regions = DMRs_sgPKCI_vs_sgC,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

# A GRanges object is returned
print(dm_annotated)





















