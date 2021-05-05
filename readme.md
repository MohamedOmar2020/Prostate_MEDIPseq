This repository contains the scripts used to analyze the MEDIP-seq data from PKCI KO (sgPKCI) versus WT (sgC) samples from prostate c42b cell line.

- Alignment to hg19: STAR v2.7.8a

- Fastq quality control: FASTqc

- BAM quality control: qualimap v.2.2.2-dev

- Analysis tools: MEDIPS and qsea (R) - macs3 (command line)

##################
## For macs3:

- Bedgraph to bigwig: ucsc-bedgraphtobigwig v4
>>Note: Chromosome sizes file can be downloaded from 
http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/latest/ 

- Convert the annotation gtf to bed file: convert2bed v: 2.4.39

- Bed annotation: bedmap v2.4.39 

#################
## Further analysis and visualization: 

annotatr and Gviz (R)
