
#cellatlas_analysis.R
#This code performs some exploratory data analysis on the RNAseq data provided by the Human Protein Atlas
#We will use the cleaned-up data file: clean_rna_cellline.txt and Bioconductor package DESeq2
#At the end we will also perform some very rudimentary analysis for fun.

######################################################################
# Bioconductor DESeq2 analysis
######################################################################

#initialize packages
library(DESeq2)
library(ggplot2)
library(RColorBrewer)

#Reading in data file and creating DEseq2 object for further analysis

celldata = read.table ( "clean_rna_cellline.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
head(celldata)

#Use the DESeq2 command, DESeqDataSetFromMatrix, to create a DEseq2 object from the raw data
#DeSeq object requires a specific format: colData (different cell type condition), and countData (values);

coldata = data.frame(colnames(celldata[,-c(1,2)]))
rownames(coldata) = coldata[,1]

#the countdata must be in dataframe format with rownames as unique genes, and colnames as the cell types from coldata

countdata = celldata[,-c(1,2)]
rownames(countdata) = celldata[,2]

#for whatever reason, the values must also be integers (no decimal places allowed)
#I will use the ceiling because there is quite a discrepancy between cells that have 0
#expression of a gene and those that have 0.1-0.9 expression units.
countdata = ceiling(countdata)

#create the DESeq2 object
ddsMat = DESeqDataSetFromMatrix (countData = countdata,
                                 colData = coldata,
                                 design = ~ 1)








