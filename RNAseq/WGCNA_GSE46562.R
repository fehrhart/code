# This script is to do WGCNA analysis with GSE46562 data
#clean workspace
rm(list=ls())

#install required packages and load them
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("WGCNA")

library(WGCNA)
options(stringsAsFactors = FALSE)

#working directory
setwd("C:/Users/friederike.ehrhart/Documents/Streamline/RNAseq stay 25.9. - 1.12.23/GSE46562")

#input data and overview
data = read.table(file = 'GSE46562_FPKM_expression_table.txt')
dim(data)
names(data)

#Transpose 
data1 = as.data.frame(t(data[, -c(1)]))
names(data1) = data$V1
rownames(data1) = names(data)[-c(1)]

#Check for genes with missing values - f the last statement returns TRUE, all genes have passed the cuts. If not use code below.
gsg = goodSamplesGenes(data1, verbose = 3);
gsg$allOK

#If FALSE do this:
#if (!gsg$allOK)
#{
  # Optionally, print the gene and sample names that were removed:
 # if (sum(!gsg$goodGenes)>0)
  #  printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  #if (sum(!gsg$goodSamples)>0)
   # printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  #datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
#}

#Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
sampleTree = hclust(dist(data1), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
