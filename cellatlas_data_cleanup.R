
#cellatlas_data_cleanup.R
#The raw data is provided as a .tsa file ("rna_cellline_copy.tsv"), which lists the TPM counts for each gene for each cell line
#This code will clean up and reorganize the data into a friendlier format for data analysis

#initialize packages
library(reshape2)

#set the working directory to where the data and scripts are stored.
setwd("~/bioinformatics/cellatlas/")
getwd()

#read in the .tsv file, tab delimited
rawData = read.table("rna_celline_copy.tsv", sep = "\t", header = T, stringsAsFactors = FALSE)
head(rawData)

#remove the last column, "unit," because all values are 'TPM'
rawData = rawData[,-5]
head(rawData)

# data is organized in 64-cell line chunks (for each gene, values for 64 cell lines)
#   Gene            Gene.name    Sample Value
# 1 ENSG00000000003    TSPAN6     A-431  27.8
# 2 ENSG00000000003    TSPAN6      A549  37.6
# 3 ENSG00000000003    TSPAN6      AF22 108.1
# 4 ENSG00000000003    TSPAN6    AN3-CA  51.8
# 5 ENSG00000000003    TSPAN6  ASC diff  32.3
# 6 ENSG00000000003    TSPAN6 ASC TERT1  17.7

#get list of cell line names
cellNames = unique(rawData$Sample)

#get a list of the gene names and the ensembl IDs
geneNames = unique(rawData$Gene.name)
ensemblNames = unique(rawData$Gene)

length(geneNames)
#[1] 19600
length(ensemblNames)
#[1] 19613

#problem here is that dim geneNames is 19600, but dim ensemblNames is 19613. there are 
#ensembl genes listed that have duplicate gene names.

#Ideally I want to have columns for each of the cell lines, and each row is a gene
#right now the genes are grouped and the "Values" column contains a looping of genes
#this StackOverflow link shows how to fix this: https://stackoverflow.com/questions/35213415/transpose-in-r-grouping-by-row-and-column
#we will use the reshape2 package, function dcast

#get rid of the gene.names column, which has non-unique names
rawEnsemblFrame = data.frame(rawData[,-2])

#use the dcast function to put cell types as columns, using unique ensembl id as rows
byCellrawData = dcast (rawEnsemblFrame, Gene~Sample, value.var = c("Value"))

#I would still like to have the gene name, so I will now
#Map the gene names to the ensembl ids
ensemblMap = array()
geneMap = array()

#testing code before running loop, making sure pairing of gene names is done correctly
rawData[which (rawData$Gene.name == 'TSPAN6'),1]
rawData[which (rawData$Gene == 'ENSG00000000003'),2]

#loop through the ensembl ids in the data and pair them with the corresponding gene name
for ( eName in byCellrawData$Gene ) {
    geneName = rawData[which (rawData$Gene == eName)[1],2]
    ensemblMap = append(ensemblMap, eName)
    geneMap = append(geneMap,geneName)
}


#combine the mapped ensembl names and gene names
nameMap = cbind(ensemblMap[-1], geneMap[-1])
colnames(nameMap) = c ("EnsemblMap", "Gene.name")

#combine the gene names with the data frame
mapped_byCell_raw = cbind (nameMap, byCellrawData)

#double check that the gene names were matched correctly
matches = array()
for ( i in c(1:length(mapped_byCell_raw$Gene)) ) {
  if ( mapped_byCell_raw$EnsemblMap[i] == mapped_byCell_raw$Gene[i] ) {
    matches = append(matches, TRUE)
  }
  else
    matches = append(matches,FALSE)
}

#check that all the gene.names match with the ensembl ids
length(which(matches == TRUE))
#[1] 19613
length(which(matches == FALSE))
#[1] 0

#write the reorganized data to a new file, deleting the redunant Ensembl names column
write.table (mapped_byCell_raw[,-1], "clean_rna_cellline.txt", sep = '\t')


