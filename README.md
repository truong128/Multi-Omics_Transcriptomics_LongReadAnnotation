# Multi-Omics_Transcriptomics_LongReadSequencing

This pipeline was developed by leveraging advanced Multiomics sequencing technologies, including RNA-seq, Oxford Nanopore Technologies (ONT) long-read sequencing, and Capillary Analysis of Gene Expression (CAGE-seq) sequencing to perform "Comprehensive resolution and classification of the Epstein Barr virus transcriptome". README files are included in each stage/folder of peak_caller, convert_CAGE, and ONT_Read_Validation for instructions on how to run them.

# Citation

Truong D Nguyen

Erik K Flemington

Link to manuscript: https://www.researchsquare.com/article/rs-5079871/v1 (Under Revision at Nature Communications)

# Requirements
Perl 5

# Dependencies
File::Basename

# Installation
```
git clone https://github.com/truong128/Multi-Omics_Transcriptomics_LongReadAnnotation.git
```

# File formats

| File    | Description     |
|:---------------|:---------------:|
|Wiggle  | Standard wiggle format |
| BED12 | Standard BED12 format  |
| BED6 | Standard BED6 format  |
| Genome fasta | Standard genome fasta (can be wrapped or unwrapped)|
