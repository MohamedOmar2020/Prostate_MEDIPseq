
rm(list = ls())

library(annotatr)
library(GenomicRanges)
library(regioneR)

########################################
## Load DMRs

# Unique to sgPKCI
DMRs_sgPKCI_vs_sgC <- read_regions(con = "./MACS/Diff_cutoff_3/diff_sgC_vs_sgPRKCI_c3.0_cond2.bed", genome = 'hg19', format = 'bed')
DMRs_sgPKCI_vs_sgC

DMRs_sgPKCI_vs_sgC <- GRanges(DMRs_sgPKCI_vs_sgC)

# Unique to sgC
DMRs_sgC_vs_sgPKCI <- read_regions(con = "./MACS/Diff_cutoff_3/diff_sgC_vs_sgPRKCI_c3.0_cond1.bed", genome = 'hg19', format = 'bed')
DMRs_sgC_vs_sgPKCI

DMRs_sgC_vs_sgPKCI <- GRanges(DMRs_sgC_vs_sgPKCI)

########################################
## Annotating Regions

# see available annotations
builtin_annotations()

# Select annotations for intersection with regions
# Note inclusion of custom annotation, and use of shortcuts
annots <- c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
           'hg19_genes_intronexonboundaries')

# Build the annotations (a single GRanges object)
annotations <- build_annotations(genome = 'hg19', annotations = annots)

# Intersect the regions we read in with the annotations
DMRs_sgPKCI_vs_sgC_annotated <- annotate_regions(
  regions = DMRs_sgPKCI_vs_sgC,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

# Coerce to a data.frame
DMRs_sgPKCI_vs_sgC_annotated <- data.frame(DMRs_sgPKCI_vs_sgC_annotated)


DMRs_sgC_vs_sgPKCI_annotated <- annotate_regions(
  regions = DMRs_sgC_vs_sgPKCI,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

# Coerce to a data.frame
DMRs_sgC_vs_sgPKCI_annotated <- data.frame(DMRs_sgC_vs_sgPKCI_annotated)

############
# Subset based on a gene symbol, in this case NOTCH1
notch1_subset <- subset(DMRs_sgPKCI_vs_sgC_annotated, annot.symbol == 'NOTCH1')
print(head(notch1_subset))

##################################################
## Randomizing Regions

# Randomize the input regions
sgPKCIvs_sgC_random_regions <- randomize_regions(
  regions = DMRs_sgPKCI_vs_sgC,
  allow.overlaps = TRUE,
  per.chromosome = TRUE)

# Annotate the random regions using the same annotations as above
# These will be used in later functions
sgPKCIvs_sgC_random_regions_annotated = annotate_regions(
  regions = sgPKCIvs_sgC_random_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = TRUE)

## the same for sgC
sgCvs_sgPKCI_random_regions <- randomize_regions(
  regions = DMRs_sgC_vs_sgPKCI,
  allow.overlaps = TRUE,
  per.chromosome = TRUE)

sgCvs_sgPKCI_random_regions_annotated = annotate_regions(
  regions = sgCvs_sgPKCI_random_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = TRUE)

###########################
## Summarizing Over Annotations

# Find the number of regions per annotation type
# and the number of random regions per annotation type
sgPKCIvs_sgC_annsum <- summarize_annotations(
  annotated_regions = DMRs_sgPKCI_vs_sgC_annotated,
  annotated_random = sgPKCIvs_sgC_random_regions_annotated,
  quiet = TRUE)

print(sgPKCIvs_sgC_annsum)


sgCvs_sgPKCI_annsum <- summarize_annotations(
  annotated_regions = DMRs_sgC_vs_sgPKCI_annotated,
  annotated_random = sgCvs_sgPKCI_random_regions_annotated,
  quiet = TRUE)

print(sgCvs_sgPKCI_annsum)

##################################################
##################################################
## Plotting

# Plotting Regions per Annotation
# View the number of regions per annotation. This function
# is useful when there is no classification or data
# associated with the regions.
annots_order = c(
  'hg19_genes_1to5kb',
  'hg19_genes_promoters',
  'hg19_genes_5UTRs',
  'hg19_genes_exons',
  'hg19_genes_intronexonboundaries',
  'hg19_genes_introns',
  'hg19_genes_3UTRs',
  'hg19_genes_intergenic')

sgPKCIvs_sgC_Plot = plot_annotation(
  annotated_regions = DMRs_sgPKCI_vs_sgC_annotated,
  annotation_order = annots_order,
  plot_title = '# unique DMRs in sgPKCI',
  x_label = 'knownGene Annotations',
  y_label = 'Count')

print(sgPKCIvs_sgC_Plot)

sgCvs_sgPKCI_Plot = plot_annotation(
  annotated_regions = DMRs_sgC_vs_sgPKCI_annotated,
  annotation_order = annots_order,
  plot_title = '# unique DMRs in sgC',
  x_label = 'knownGene Annotations',
  y_label = 'Count')

print(sgCvs_sgPKCI_Plot)

#############################
## Plotting Regions Occurring in Pairs of Annotations

# View a heatmap of regions occurring in pairs of annotations
annots_order = c(
  'hg19_genes_promoters',
  'hg19_genes_5UTRs',
  'hg19_genes_exons',
  'hg19_genes_introns',
  'hg19_genes_3UTRs',
  'hg19_genes_intergenic')

sgPKCIvs_sgC_coannotations <- plot_coannotations(
  annotated_regions = DMRs_sgPKCI_vs_sgC_annotated,
  annotation_order = annots_order,
  axes_label = 'Annotations',
  plot_title = 'Regions in Pairs of Annotations (sgPKCI)')

print(sgPKCIvs_sgC_coannotations)


sgCvs_sgPKCI_coannotations <- plot_coannotations(
  annotated_regions = DMRs_sgC_vs_sgPKCI_annotated,
  annotation_order = annots_order,
  axes_label = 'Annotations',
  plot_title = 'Regions in Pairs of Annotations (sgC)')

print(sgCvs_sgPKCI_coannotations)

########################################
## Plotting Numerical Data Over Regions

sgPKCIvs_sgC_score <- plot_numerical(
  annotated_regions = DMRs_sgPKCI_vs_sgC_annotated,
  x = 'score',
  facet = 'annot.type',
  facet_order = c('hg19_genes_1to5kb','hg19_genes_promoters',
                  'hg19_genes_5UTRs','hg19_genes_3UTRs',
                  'hg19_genes_intergenic', 'hg19_cpg_islands'),
  bin_width = 5,
  plot_title = 'sgPKCI Region Methylation In Genes',
  x_label = 'sgPKCI')

print(sgPKCIvs_sgC_score)


sgCvs_sgPKCI_score <- plot_numerical(
  annotated_regions = DMRs_sgC_vs_sgPKCI_annotated,
  x = 'score',
  facet = 'annot.type',
  facet_order = c('hg19_genes_1to5kb','hg19_genes_promoters',
                  'hg19_genes_5UTRs','hg19_genes_3UTRs',
                  'hg19_genes_intergenic', 'hg19_cpg_islands'),
  bin_width = 5,
  plot_title = 'sgC Region Methylation In Genes',
  x_label = 'sgC')

print(sgCvs_sgPKCI_score)









