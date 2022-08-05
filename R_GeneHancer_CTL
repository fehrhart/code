#clean workspace
rm(list=ls())

setwd("C:/Users/friederike.ehrhart/Documents/actual_linksets/GeneHancer")

library(tidyr)
library(dplyr)
library(stringr)

#import data
TFBS = read.table(file = 'GeneHancer_TFBSs_v5.10.txt', sep = '\t', header = TRUE)
GeneTargets = read.table(file = 'GeneHancer_v5.10.gff', sep = '\t', header = FALSE)

# remove not needed columns from GeneTargets
GeneTargets1 <- subset(GeneTargets, select = -c(V1,V2,V4,V5,V6,V7,V8))

# remove unnecessary text from string
GeneTargets1$V9<-gsub("genehancer_id=","",as.character(GeneTargets1$V9))
GeneTargets1$V9<-gsub("connected_gene=","",as.character(GeneTargets1$V9))
GeneTargets1$V9<-gsub(";score=","=",as.character(GeneTargets1$V9))
# create separate column for GHids and then remove it from string in V9
GeneTargets1$V10=substr(GeneTargets1$V9,1,11)
GeneTargets1$V9=str_sub(GeneTargets1$V9, 13,)

# break string in V9 in separate rows
GeneTargets2 <- separate_rows(GeneTargets1,V9,sep = ";")

# break genes and scores in separate columns
GeneTargets3=separate(GeneTargets2,V9, c("GeneTarget", "Score"), sep = "=")
GeneTargets3$Score=as.numeric(GeneTargets3$Score)

#change column names
colnames(GeneTargets3) <-c('Regulator','GeneTarget','Score','GHid')

#map gene names to NCBI gene identifiers
#download mapping file from biomart
library(biomaRt)
mart <- useMart("ensembl")
#how to identify datasets and attributes
#datasets <- listDatasets(mart)
#attributes = listAttributes(ensembl)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#query a list of hsa genes with 2 attributes 
Ensembl_hsa_genes <- getBM(attributes = c("external_gene_name","entrezgene_id"),mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"))

#change column name to be able to join them
colnames(Ensembl_hsa_genes) <-c('GeneTarget','NCBIID')
head(Ensembl_hsa_genes)

#merge
GeneTargets3E <- merge(GeneTargets3, Ensembl_hsa_genes, by="GeneTarget")

# create small subsets for merge testing
#GeneTargetsSmall <- GeneTargets3E[1:5,]
#GeneTargetsSmall$GHid=c('GH01J149156','GH01J149156','GH01J149156','GH01J149156','GH01J149156')
#TFBSmall <- TFBS[1:5,]
#TFBSmall$GHid=c('GH01J149156','GH01J149156','GH01J149156','GH01J149156','GH01J149100')
#totalsmall <- merge(GeneTargetsSmall, TFBSmall, by="GHid")

#merge ALL
#total <- merge(GeneTargets3E, TFBS, by="GHid")

# identify score threshold
#ScoreCheck <- as.numeric(unlist(GeneTargets3E[ ,"Score"]))
#hist(ScoreCheck, breaks = 250, col = "black")

#select by score
GeneTargets4 <- subset(GeneTargets3E, Score > 5)

#select by score
GeneTargets5 <- subset(GeneTargets3E, Score > 200)

#merge with score filter
totalScore5 <- merge(GeneTargets4, TFBS, by="GHid")

#export table as tab separated text file, without "" and row names
write.table(totalScore5, file="GeneHancerScore5_input.txt", sep = "\t", quote=FALSE, row.names = FALSE)
