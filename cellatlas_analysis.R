
#cellatlas_analysis.R
#This code performs some exploratory data analysis on the RNAseq data provided by the Human Protein Atlas
#We will use the cleaned-up data file: clean_rna_cellline.txt and Bioconductor package DESeq2
#Will visualize the data using PCA plot and Volcano plot, and also identify differentially expressed genes (csv file).
#We will also repeat our analysis, narrowed down on enzyme expression

######################################################################
# Bioconductor DESeq2 analysis
# Info can be found at: https://www.bioconductor.org/help/workflows/rnaseqGene/
######################################################################

#initialize packages
library("DESeq2")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("ggrepel") 

setwd("~/cellatlas")

#Reading in data file and creating DEseq2 object for further analysis
celldata = read.table ( "data_clean_rna_cellline.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
head(celldata)

#Use the DESeq2 command, DESeqDataSetFromMatrix, to create a DEseq2 object from the raw data
#DeSeq object requires a specific format: colData (different cell type condition), and countData (values);

#Making the coldata parameter:
#first categorize each cell line as hematopoietic vs non-hematopoietically derived:

#Find indicies of all the lines which are hematopoietically derived:
hemato = c( "HEL","NB.4","HAP1","HL.60","HMC.1","K.562","THP.1","U.937",
            "REH","Daudi","HDLM.2","Karpas.707","MOLT.4","RPMI.8226","U.266.70","U.266.84","U.698" )
hematoIndex = which(cellLines %in% hemato)

#Here are some classifications for the other cell lines. During this analysis I will only focus on
#the hematopoietically derived vs non-hematopoietically derived lines.

# liver = c( "Hep.G2", "HHSteC")
# brain = c( "AF22", "SH.SY5Y","U.251.MG","U.138.MG","U.87.MG")
# urinary = c( "HEK.293", "RPTEC.TERT1","RT4","NTERA-2","PC.3")
# skin = c( "A.431", "HaCaT","SK.MEL.30","WM.115")
# cervical = c( "HeLa", "SiHa")
# lung = c( "A549", "HBEC3.KT","SCLC.21H")

#make a data.frame categorizing each cell line with their respective cell type
celltype = rep (c("misc"), 64)
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
# vsd <- vst(dds, blind = FALSE)
# head(assay(vsd), 3)
# save(vsd, file = "vsd_dds.Robj")


#Use the rlog transformed data to cluster the cell lines based on similarity
plotDists = function (rld.obj) {
  sampleDists <- dist(t(assay(rld.obj)))
  
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste( rld.obj$celltype )
  
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
}

plotDists(rld)

#https://support.bioconductor.org/p/90791/
#plotPCA() shows the PCA plot of the log transformed data
#This code adds the names of the samples ontop of the PCA plot dots
name.plotPCA = function (rld.obj) {
  p <- plotPCA(rld.obj,  intgroup = c("celltype"))
  p <- p + geom_text_repel(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
  p = p+ ggtitle("PCA Plot - Clustering analysis of
                 \nHematopoietic vs Non-Hematopoietic cell lines")
  print(p)
}

name.plotPCA(rld)

#can observe that the cell lines: Hep.G2, HDLM.2, HMC.1 and SH.SY5Y have transcriptional profiles
#that vastly differ from the majority of the cells

rld_removedoutliers = rld[,-c(21,25,28, 48)]
name.plotPCA(rld_removedoutliers)

#now we observe that ASC.diff, MOLT.4, and HEL are also quite separated from the pack
#also remove the BJ-derivative cell lines, since these all originate from a common line (the BJ line)
#and we would be over-representing these cells in the dataset in the context of this analysis

rld_removedoutliers_2nd = rld[,-c(5,9,10,11,18,21,23,25,28,38,48)]
name.plotPCA(rld_removedoutliers_2nd)

#re-running the sample dists on the outlier removed matrix:
plotDists(rld_removedoutliers_2nd)

######################################################################
# Performing differential gene analysis between cell lines / groups
######################################################################

#This requires performing the analysis on the raw counts (the 'dds' object) and the DESeq function
#DEseq performs estimateSizeFactors() and estimateDispersions(), so pre-running these commands are not necessary
#finally, will perform the analysis on only the cell lines of interest (+/- HepG2);

#use the non-outlier cells for the analysis
ddsclean = dds[,-c(5,9,10,11,18,21,23,25,28,38,48)]
ddsclean = DESeq (ddsclean)
res_clean = results(ddsclean)

# Sample output:
# log2 fold change (MLE): celltype misc vs hemato 
# Wald test p-value: celltype misc vs hemato 
# DataFrame with 18738 rows and 6 columns
# baseMean log2FoldChange     lfcSE        stat       pvalue         padj
# <numeric>      <numeric> <numeric>   <numeric>    <numeric>    <numeric>
#   ENSG00000000003    291.99934    6.626358377 0.4757870 13.92715323 4.332839e-44 2.854648e-41
# ENSG00000000005      0.00000    0.000000000 0.0000000  0.00000000 1.000000e+00           NA
# ENSG00000000419    859.52378    0.005919813 0.1760574  0.03362433 9.731767e-01 9.832052e-01
# ENSG00000000457     56.76461   -0.475573276 0.1599503 -2.97325568 2.946589e-03 1.079717e-02
# ENSG00000000460    155.16068   -0.541859654 0.2926190 -1.85175800 6.406058e-02 1.336638e-01
# ...                      ...            ...       ...         ...          ...          ...
# ENSG00000284512   3.40815722    -0.57261642 0.4252183 -1.34664092    0.1780959    0.2944312
# ENSG00000284526 128.84405899    -0.24041343 0.6612565 -0.36357060    0.7161787    0.8075150
# ENSG00000284546   0.05234491     0.09012843 3.5115811  0.02566605    0.9795237           NA
# ENSG00000284552   0.04810316    -0.60188290 3.5115811 -0.17139940    0.8639097           NA
# ENSG00000284554   3.88451639    -0.86557310 1.1975708 -0.72277405    0.4698187    0.6008764

#run the analysis on subset of the samples:
dds_select = dds[,c(14,22, 24, 25, 53)]
dds_select = DESeq(dds_select)
res_select = results (dds_select)
#if you want to include the raw counts side by side with the results
#(so you can verify the values per sample for interesting genes):
res_select_withcount = cbind(res_select, assay(dds[,c(14, 53,22, 24, 25)]))


dds_select_noHep = dds[,c(14,22, 24, 53)]
dds_select_noHep = DESeq(dds_select_noHep)
res_select_noHep = results (dds_select_noHep)
res_select_noHep_withcount = cbind(res_select_noHep, assay(dds[,c(14, 53,22, 24)]))


#run the analysis all all the samples (including outliers)
dds_full = DESeq(dds)
res_full = results(dds_full)

# 4 result tables:
# 1. res_full : all the cell lines, comparing hematopoietic vs non-hematopoietic
# 2. res_clean : removed outlier cell lines shown in the PCA plot
# 3. res_select : kept only a small subset of the cell lines
# 4. res_select_noHep : the small subset with the Hep.G2 cell line removed

#now we can export the files 
#we want to annotate the files with gene symbol, annotated enzymes, and GO terms
#https://www.bioconductor.org/packages/3.7/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf

library("AnnotationDbi")
library("org.Hs.eg.db")

#This function adds columns representing gene symbol, enzyme ID, and GO term
my.mapids = function (res) {
  res$symbol <- mapIds(org.Hs.eg.db,
                       keys=row.names(res),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
  res$enzyme <- mapIds(org.Hs.eg.db,
                       keys=row.names(res),
                       column="ENZYME",
                       keytype="ENSEMBL",
                       multiVals="first")
  res$go <- mapIds(org.Hs.eg.db,
                   keys=row.names(res),
                   column="GO",
                   keytype="ENSEMBL",
                   multiVals="first")
  resOrdered <- res[order(res$pvalue),]
  return(resOrdered)
}

#Run the function my.mapids()
res_clean = my.mapids(res_clean)
res_full = my.mapids(res_full)
res_select = my.mapids(res_select)
res_select_noHep = my.mapids(res_select_noHep)

#save the RObjects to files incase we want to load them for analysis later on.
save(res_clean, file = "res_clean.Robj")
save(res_full, file = "res_full.Robj")
save(res_select_withcount, file = "res_select.Robj")
save(res_select_noHep_withcount, file = "res_select_noHep.Robj")

#export gene lists, top 2000 differentially expressed genes
resOrderedDF <- as.data.frame(res_clean)[1:2000, ]
write.csv(resOrderedDF, file = "results_hemato_vs_non_DEG.csv")

resOrderedDF <- as.data.frame(res_full)[1:2000, ]
write.csv(resOrderedDF, file = "res_full_results.csv")

resOrderedDF <- as.data.frame(res_select_withcount)[1:2000, ]
write.csv(resOrderedDF, file = "res_select_results.csv")

resOrderedDF <- as.data.frame(res_select_noHep_withcount)[1:2000, ]
write.csv(resOrderedDF, file = "res_select_noHep_results.csv")


# Volcano plot to visualize top 20 differentially expressed genes in the cleaned up dataset

plot.volcano = function (res) {
  input <- mutate(data.frame(res), sig=ifelse(data.frame(res)$padj<0.0001, "padj < 0.0001", "Not Sig"))
  #dim(input)
  input = input[!is.na(input$sig),]
  #dim(input)
  volc = ggplot(input, aes(log2FoldChange, -log10(pvalue))) +
    geom_point(aes(col=sig)) + 
    scale_color_manual(values=c("black", "red")) + 
    ggtitle("Volcano Plot - Differentially Expressed Genes 
            \nEnriched in hematopoietic (left) vs non-hematopoietic (right) cell lines") +
    geom_text_repel(data=head(input, 20), aes(label=symbol))
  print(volc)
}

plot.volcano(res_clean)


######################################################################
# Repeat the analysis, this time performing Enzyme-only comparison
######################################################################

map.enzyme = function (res) {
  res$enzyme <- mapIds(org.Hs.eg.db,
                       keys=row.names(res),
                       column="ENZYME",
                       keytype="ENSEMBL",
                       multiVals="first")
  return(res)
}

#tag all the ensembl genes that have an enzyme entry
enzyme_countdata10 = map.enzyme(countdata10)

#remove all the ensembl genes that are missing an enzyme entry, leaving only annotated enzymes in the matrix
enzyme_countdata10 = na.omit(enzyme_countdata10)

dim(countdata10)
#[1] 19613    64
dim(enzyme_countdata10)
#[1] 2222   65

#it looks like there were 2222 annotated enzymes. now create the DESeq object:
dds_e = DESeqDataSetFromMatrix (countData = enzyme_countdata10[,-65],
                                colData = coldata,
                                design = ~ celltype)

#filter out rows (enzymes) with no reads at all:
nrow(dds_e)
# [1] 2222
dds_e <- dds_e[ rowSums(counts(dds_e)) > 1, ]
nrow(dds_e)
# [1] 2200, around 22 enzymes filtered out

# perform regularized-log transformation before plotting data on PCA plot:
rld_e <- rlog(dds_e, blind = FALSE)
head(assay(rld_e), 3)
save(rld_e, file = "rld_e_dds.Robj")

plotDists(rld_e)
name.plotPCA(rld_e)
#observe less variability when only looking at enzymes - try removing ASC.diff, A549, RT4, Hep.G2

rld_e_removedoutliers = rld_e[,-c(2,5,25, 46)]
name.plotPCA(rld_e_removedoutliers)

#now we observe that THP.1 is highy separated from pack

rld_e_removedoutliers_2nd = rld_e[,-c(2,5, 18, 23, 25, 46, 53)]
name.plotPCA(rld_e_removedoutliers_2nd)

#re-running the sample dists on the outlier removed matrix:
plotDists(rld_e_removedoutliers_2nd)


#############################################################################
# Perform DESeq2 differential gene analysis on the enzyme DEseq object
#############################################################################

ddsclean = dds_e[,-c(2,5,9,10,11, 18, 23, 25, 46, 53)]
ddsclean = DESeq (ddsclean)
res_clean = results(ddsclean)

dds_select = dds_e[,c(14,22, 24, 25, 53)]
dds_select = DESeq(dds_select)
res_select = results (dds_select)

head(assay(dds_e[,c(14, 53,22, 24, 25)]))
res_select_withcount = cbind(res_select, assay(dds_e[,c(14, 53,22, 24, 25)]))

dds_select_noHep = dds_e[,c(14,22, 24, 53)]
dds_select_noHep = DESeq(dds_select_noHep)
res_select_noHep = results (dds_select_noHep)
res_select_noHep_withcount = cbind(res_select_noHep, assay(dds_e[,c(14, 53,22, 24)]))

dds_full = DESeq(dds_e)
res_full = results(dds_full)

res_clean = my.mapids(res_clean)
res_full = my.mapids(res_full)
res_select_withcount = my.mapids(res_select_withcount)
res_select_noHep_withcount = my.mapids(res_select_noHep_withcount)

save(res_clean, file = "res_clean_enzonly.Robj")
save(res_full, file = "res_full_enzonly.Robj")
save(res_select_withcount, file = "res_select_enzonly_withcount.Robj")
save(res_select_noHep_withcount, file = "res_select_noHep_enzonly_withcount.Robj")

#export sorted enzyme lists:
write.csv(res_clean, file = "enzymeonly_res_clean_results.csv")
write.csv(res_full, file = "enzymeonly_res_full_results.csv")
write.csv(res_select_withcount, file = "enzymeonly_res_select_results_withcount.csv")
write.csv(res_select_noHep_withcount, file = "enzymeonly_res_select_noHep_results_withcount.csv")

plot.volcano(res_clean)
plot.volcano(res_select_noHep_withcount)

#taking a closer look at the data - some genes are greatly enriched in only one of the non-hematopoetic cell types
#this causes the algorithm to pull it out as a differentially expressed gene. for example:

res_select_withcount[which (res_select_withcount$Hep.G2 > 800),]

#                  symbol    Daudi     THP.1     HEK.293  HeLa    Hep.G2
#ENSG00000115718     PROC   0         0         2        48      1055   

#may be useful to then filter the results such that there is a minimum cutoff required for each sample
#here is how to filter out genes that have less than 10 transcripts for a certain cell line:

filtered = res_select_withcount[which (res_select_withcount$Hep.G2 > 10),]


sessionInfo()
#######
# Session Info
######
# R version 3.4.2 (2017-09-28)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 16.04.3 LTS
# 
# Matrix products: default
# BLAS: /usr/lib/libblas/libblas.so.3.6.0
# LAPACK: /usr/lib/lapack/liblapack.so.3.6.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pheatmap_1.0.8             dplyr_0.7.4                RColorBrewer_1.1-2         ggplot2_2.2.1              DESeq2_1.18.1             
# [6] SummarizedExperiment_1.8.1 DelayedArray_0.4.1         matrixStats_0.53.1         Biobase_2.38.0             GenomicRanges_1.30.2      
# [11] GenomeInfoDb_1.14.0        IRanges_2.12.0             S4Vectors_0.16.0           BiocGenerics_0.24.0        BiocInstaller_1.28.0      
# 
# loaded via a namespace (and not attached):
#   [1] bit64_0.9-7            splines_3.4.2          Formula_1.2-2          assertthat_0.2.0       latticeExtra_0.6-28    blob_1.1.0            
# [7] GenomeInfoDbData_1.0.0 pillar_1.1.0           RSQLite_2.0            backports_1.1.2        lattice_0.20-35        glue_1.2.0            
# [13] digest_0.6.15          XVector_0.18.0         checkmate_1.8.5        colorspace_1.3-2       htmltools_0.3.6        Matrix_1.2-12         
# [19] plyr_1.8.4             XML_3.98-1.10          pkgconfig_2.0.1        genefilter_1.60.0      zlibbioc_1.24.0        xtable_1.8-2          
# [25] scales_0.5.0           BiocParallel_1.12.0    htmlTable_1.11.2       tibble_1.4.2           annotate_1.56.1        nnet_7.3-12           
# [31] lazyeval_0.2.1         survival_2.41-3        magrittr_1.5           memoise_1.1.0          foreign_0.8-69         tools_3.4.2           
# [37] data.table_1.10.4-3    stringr_1.3.0          munsell_0.4.3          locfit_1.5-9.1         cluster_2.0.6          AnnotationDbi_1.40.0  
# [43] bindrcpp_0.2           compiler_3.4.2         rlang_0.2.0            grid_3.4.2             RCurl_1.95-4.10        rstudioapi_0.7        
# [49] htmlwidgets_1.0        labeling_0.3           bitops_1.0-6           base64enc_0.1-3        gtable_0.2.0           DBI_0.7               
# [55] reshape2_1.4.3         R6_2.2.2               gridExtra_2.3          knitr_1.20             bit_1.1-12             bindr_0.1             
# [61] Hmisc_4.1-1            stringi_1.1.6          Rcpp_0.12.15           geneplotter_1.56.0     rpart_4.1-12           acepack_1.4.1






