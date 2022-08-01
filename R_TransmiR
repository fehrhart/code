#clean workspace
rm(list=ls())

setwd("C:/Users/friederike.ehrhart/Documents/actual_linksets/miRNA linkset update 14.7.22")

#import data
input = read.table(file = 'mmu.txt', sep = '\t', header = FALSE)
#select those rows with evidence "level 2" OR "literature"
input3 <- subset(input, V7=="level 2"|V7=="literature")
# remove not needed columns
input4 <- subset(input3, select = -c(V3,V4,V9))
#change column names
colnames(input4) <-c('Gene','miRNA','Interaction','Interaction_evidence','Evidence_score','Tisse')
#check result
head(input4)

#export table as tab separated text file, without "" and row names
#write.table(input4, file="TransmiR_input.txt", sep = "\t", quote=FALSE, row.names = FALSE)

#download mapping file from biomart
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

#how to identify datasets and attributes 
#attributes = listAttributes(ensembl)
#datasets <- listDatasets(mart)
#mart <- useMart("ensembl")

#query a list of mmu genes with 2 attributes 
Ensembl_mmu_genes <- getBM(attributes = c("external_gene_name","entrezgene_id"),mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"))

#change column nambe to be able to join them
colnames(Ensembl_mmu_genes) <-c('Gene','ID')
head(Ensembl_mmu_genes)

#merge
total <- merge(input4, Ensembl_mmu_genes, by="Gene")

#export table as tab separated text file, without "" and row names
write.table(total, file="mmuinput2.txt", sep = "\t", quote=FALSE, row.names = FALSE)
