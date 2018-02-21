
#cellatlas_analysis.R
#This code performs some exploratory data analysis on the RNAseq data provided by the Human Protein Atlas
#We will use the cleaned-up data file: clean_rna_cellline.txt and Bioconductor package DESeq2
#At the end we will also perform some very rudimentary analysis for fun.

######################################################################
# Bioconductor DESeq2 analysis
# Info can be found at: https://www.bioconductor.org/help/workflows/rnaseqGene/
# resources: https://www.biostars.org/p/152033/
######################################################################

#initialize packages
library("DESeq2")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")

setwd("~/cellatlas_erick")

#Reading in data file and creating DEseq2 object for further analysis

celldata = read.table ( "clean_rna_cellline.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
head(celldata)

#Use the DESeq2 command, DESeqDataSetFromMatrix, to create a DEseq2 object from the raw data
#DeSeq object requires a specific format: colData (different cell type condition), and countData (values);

#Making the coldata parameter:
#categorize each cell line as hematopoietic vs non-hematopoietically derived:

#Find indicies of all the lines which are hematopoietically derived:
hemato = c("Daudi", "U.698", "MOLT.4", "REH","HEL","K.562","HMC.1","HL.60","U.937","NB.4","THP.1" )
hematoIndex = which(cellLines %in% hemato)

#make a data.frame categorizing each cell line with their respective cell type
celltype = rep (c("non"), 64)
celltype [hematoIndex]= "hemato"

coldata = data.frame(factor(celltype))
rownames(coldata) = colnames(celldata[,-c(1,2)])
colnames(coldata) = "celltype"

#the countdata must be in dataframe format with rownames as unique genes, and colnames as the cell types from coldata
countdata = celldata[,-c(1,2)]
rownames(countdata) = celldata[,2]

#for whatever reason, the values must also be integers (no decimal places allowed)

#Two solutions:
# (1) Can round up using ceiling() because there is quite a discrepancy between cells that have 0
#expression of a gene and those that have 0.1-0.9 expression units.
countdata = ceiling(countdata)

# (2) another way to do this would be to multiply all values by 10, preserving the decimal place
countdata10 = countdata*10

#create the DESeq2 object, indicating that the cell lines are categorized by cell type
#modified this to include categorization by cell name as well as type
dds = DESeqDataSetFromMatrix (countData = countdata10,
                              colData = coldata,
                              design = ~ celltype)

#filter out rows (genes) with no reads at all:
nrow(dds)
# [1] 19613
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
# [1] 18738, around 1K genes filtered out

# can perform regularized-log transformation before plotting data on PCA plot:
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
save(rld, file = "rld_dds.Robj")

#alternatively you could use a variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
save(vsd, file = "vsd_dds.Robj")



sampleDists <- dist(t(assay(rld)))
sampleDists



#https://support.bioconductor.org/p/90791/
plotPCA(rld, intgroup = c("celltype"))


######################################################################
# Performing differential gene analysis between cell lines / groups
######################################################################

#This requires performing the analysis on the raw counts (the 'dds' object)
#I have seen others performing estimateSizeFactors() and estimateDispersions() before running DESeq
#I will test whether the same results are obtained with or without pre-performing these functions.

dds

# class: DESeqDataSet 
# dim: 18738 64 
# metadata(1): version
# assays(1): counts
# rownames(18738): ENSG00000000003 ENSG00000000005 ... ENSG00000284552 ENSG00000284554
# rowData names(0):
#   colnames(64): A.431 A549 ... U.937 WM.115
# colData names(1): celltype

#This estimates size factors to account for sequencing depth
dds <- estimateSizeFactors(dds)

# class: DESeqDataSet 
# dim: 18738 64 
# metadata(1): version
# assays(1): counts
# rownames(18738): ENSG00000000003 ENSG00000000005 ... ENSG00000284552 ENSG00000284554
# rowData names(0):
#   colnames(64): A.431 A549 ... U.937 WM.115
# colData names(2): celltype sizeFactor

#obtain dispersion estimates for negative binomial distributed data
dds <- estimateDispersions(dds)

# class: DESeqDataSet 
# dim: 18738 64 
# metadata(1): version
# assays(2): counts mu
# rownames(18738): ENSG00000000003 ENSG00000000005 ... ENSG00000284552 ENSG00000284554
# rowData names(9): baseMean baseVar ... dispOutlier dispMAP
# colnames(64): A.431 A549 ... U.937 WM.115
# colData names(2): celltype sizeFactor







