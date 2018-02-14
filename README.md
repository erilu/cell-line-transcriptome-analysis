# R-Cell-Line-HPA-Transcript-Analysis

**Use R to analyze RNA-seq data from 64 cell lines from The Human Protein Atlas**

The RNA-seq data can be downloaded from the [Human Protein Atlas website](https://www.proteinatlas.org/about/download). Look under item **4 RNA Gene Data**. A copy of the cell line RNA data is also included in this repository. The data contains RNA-seq read counts from 64 different human cell lines. 

A description of the lines can be found here: (https://www.proteinatlas.org/humancell/cell+line). This is a very useful resource for examining gene expression differences between cell lines from different origin tissues. These include many commonly used cell lines--researchers using these cell lines for experimental purposes will likely want to know what genes they express and to what extent they express them at.

This is an R script that will:

1. Clean and organize the data. Currently, the data is provided as a tab-delimited file with dimensions [1255232 x 5]. There are 64 cell lines, and roughly 19,600 genes per line. For each gene, the expression values for each cell line are provided in consecutive rows. This is not a great way to organize the data, since some researchers may want to quickly pull out the gene expression profiles for only a few of their cell lines of interest. I will transform the data such that the expression values for the 19,000+ genes will be categorized by cell line, by column. The dimensions of the cleaned up data are [19613 x 66].

2. Explore the differences in transcription between cell lines of interest. Pull out commonly expressed genes as well as differentially expressed genes.
