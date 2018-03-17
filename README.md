# Cell Line Transcriptome Analysis Using R

**Use R to clean, analyze, and visualize RNA-seq data from 64 cell lines from The Human Protein Atlas**

---
## The data

The RNA-seq data can be downloaded from the [Human Protein Atlas website](https://www.proteinatlas.org/about/download). Look under item **4 RNA Gene Data**. A copy of the cell line RNA data is also included in this repository,  [data_rna_celline_copy](https://github.com/erilu/R-Cell-Line-Transcriptome-Analysis/blob/master/data_rna_celline_copy.tsv). The data contains RNA-seq read counts from 64 different human cell lines.

A description of the lines can be found here: (https://www.proteinatlas.org/humancell/cell+line). This is a very useful resource for examining gene expression differences between cell lines from different origin tissues. These include many commonly used cell lines--researchers using these cell lines for experimental purposes will likely want to know what genes they express and to what extent they express them at. This project will show you how to identify differentially expressed genes between different groups of cell lines, which is a useful technique that can be used to identify enriched pathways and possible drug targets.

---
## Organizing the data

Currently, the data is provided as a tab-delimited file with dimensions [1255232 x 5]. There are 64 cell lines, and roughly 19,600 genes per line. For each gene, the expression values for each cell line are provided in consecutive rows:

```
# data is organized in the "long" format
#   Gene            Gene.name    Sample Value
# 1 ENSG00000000003    TSPAN6     A-431  27.8
# 2 ENSG00000000003    TSPAN6      A549  37.6
# 3 ENSG00000000003    TSPAN6      AF22 108.1
# 4 ENSG00000000003    TSPAN6    AN3-CA  51.8
# 5 ENSG00000000003    TSPAN6  ASC diff  32.3
# 6 ENSG00000000003    TSPAN6 ASC TERT1  17.7
...
```

Although the way the data is organized might be useful for certain data analysis packages, some researchers may want to quickly pull out the gene expression profiles for only a few of their cell lines of interest. Furthermore, the popular RNA-seq data analysis tool [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) from [Bioconductor](https://www.bioconductor.org/help/workflows/rnaseqGene/) will require the data in a different format. I will transform the data such that the expression values for the 19,000+ genes will be categorized by cell line, by column. The dimensions of the cleaned up data are [19613 x 66]. Using the dcast() function from reshape2 will do the trick (see file to understand variable names).

```
# use the dcast function to put cell types as columns, using unique ensembl id as rows
byCellrawData = dcast (rawEnsemblFrame, Gene~Sample, value.var = c("Value"))
```
This results in the reorganized data file, [data_clean_rna_cellline.txt](https://github.com/erilu/R-Cell-Line-Transcriptome-Analysis/blob/master/data_clean_rna_cellline.txt).

---

## Differential Gene Analysis

We can now read in the cleaned up data file and explore the differences in transcription between cell lines of interest. We will cluster the cell lines and pull out differentially expressed genes between subtypes of cells. To do this, we will use the [Bioconductor](https://www.bioconductor.org/help/workflows/rnaseqGene/) package [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). This package has the capability to read in raw data files from each individual biological replicate, as well as the ability to analyze the read counts from an already compiled data file (which is what the Human Protein Atlas provides).

The file [cellatlas_analysis.R](https://github.com/erilu/R-Cell-Line-Transcriptome-Analysis/blob/master/cellatlas_analysis.R) will take you through how to perform the analysis (reading the data in, making a DESeq object, and annotating, exporting, and plotting the results). The exported list of differentially expressed genes between hematopoietic cells and non-hematopoietic cels is called [results_hemato_vs_non_DEGs.csv](https://github.com/erilu/R-Cell-Line-Transcriptome-Analysis/blob/master/results_hemato_vs_non_DEGs.csv). Below are some sample plots that can be used to visualize the data.

### PCA clustering of cell lines

![PCA clustering of cell lines](https://github.com/erilu/R-Cell-Line-Transcriptome-Analysis/blob/master/results_PCA_cluster_hemato_vs_non.png)

Above is a principal components analysis (PCA) plot on a subset of the cell lines in the dataset. We can observe that the hematopoietic cell lines cluster away from the non-hematopoietic cell lines in the PCA plot. This suggests that they have different gene expression profiles. The next plot will show us some of the genes that contribute to the clustering we see here.

### Volcano plot to visualize differentially expressed genes with p-value cutoff

![Volcano plot cell line DEGs](https://github.com/erilu/R-Cell-Line-Transcriptome-Analysis/blob/master/results_volcano_plot_DEGs.png)

In accord with the clustering analysis, there are a lot of genes that are differentially expressed in hematopoietic vs non-hematopoietic cells.

### Heatmap to display top differentially expressed genes

![Heatmap of top DEGs](https://github.com/erilu/R-Cell-Line-Transcriptome-Analysis/blob/master/results_heatmap_top50_DEGs_ggplot2.png)

The differentially expressed genes can also be visualized using a heatmap.

---

If you are interested in learning how to perform a full RNA-seq pipeline analysis, you can look at my other [repo](https://github.com/erilu/Complete-RNA-seq-Pipeline-Transcriptome-Analysis) where I align raw .fastq sequencing files to a mouse reference genome, then use Bioconductor to find differentially expressed genes in activated vs. un-activated dendritic cells.
