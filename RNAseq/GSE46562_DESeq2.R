# This script is to analyse GSE46562 data with DESeq2
#clean workspace
rm(list=ls())

#install required packages and load them
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

#working directory
setwd("C:/Users/friederike.ehrhart/Documents/Streamline/RNAseq stay 25.9. - 1.12.23/GSE46562")

library(DESeq2)
library(ggplot2)

#input data and overview
data = read.table(file = 'GSE46562_FPKM_expression_table.txt', row.names = 1, header = TRUE)
summary(data)
boxplot(log(data))

#create metadata file
metadata <-cbind(sampleID=colnames(data),Condition=c("control", "control", "patient", "control", "patient", "patient", "patient", "patient", "patient", "control", "control", "patient", "patient", "patient", "patient", "control", "control", "patient", "patient"))

#create "dds" item from data and metadata - needs round function as DESeq2 works with non-normalised counts which are integers
dds <- DESeqDataSetFromMatrix(countData = round(data),
                              colData = metadata,
                             design = ~ Condition, tidy = FALSE)
#run function
dds <- DESeq(dds)
#output results
res <- results(dds)
head(results(dds, tidy=TRUE))

#sort by log2FC
res <- res[order(res$log2FoldChange),]
head(res)

#export to CSV
write.csv(res, file = "results.csv")

#Volcano plot
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,6)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
