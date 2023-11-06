#Monte Carlo simulation test "small" - this one compares how many genes (exact number) overlap with the BDNF pathway (WP2380) when replaced by a random selection of all human genes available in WP
#clean workspace
rm(list=ls())
setwd("C:/...")

library(dplyr)

#1. import data
#import gene-pathway list derived from SPARQL query
countlist <- read.csv("GeneID_independent.csv", header= TRUE, fill = TRUE)  
as.data.frame.matrix(countlist)

#subset list of available genes 
Available_genes <- (countlist$GeneName)

#2. create original subsets
WP5412 <- data.frame(c("PTK2B",	"ATG101",	"TP53",	"TNF",	"TRAF2",	"GABARAPL2",	"RB1",	"RB1CC1",	"ULK1",	"MAP3K5",	"GABARAP",	"ATG16L1",	"ATG13",	"PTK2",	"WDR45B",	"GABARAPL1",	"RB1CC1",	"RB1"))
colnames(WP5412) <-("gene")

WP2380 <- data.frame(c("ACACB","ADAM17","AKT1","ALPL","APC","BAD","BCL2L11","BDNF","BMP2","CAMK1","CAMK2A","CAMK4","CASP3","CDC42","CDH2","CDK5","CDK5R1","CDKL5","CFL1","CHUK","CNR1","CREB1","CRTC1","CSNK2A1","CTNNB1","CYFIP1","DLG1","DOCK3","DOK5","DPYSL2","EEF2","EGR1","EGR2","EIF2S1","EIF2S2","EIF4E","EIF4EBP1","ELK1","FOS","FOXO3","FRS2","FRS3","FYN","GABRB3","GRB2","GRIA1","GRIA2","GRIA3","GRIN1","GRIN2B","GRIP1","GSK3B","HRAS","IGF2BP1","IKBKB","IKBKG","IRS1","IRS2","JAK2","JUN","KCNA3","KCNN2","KIDINS220","KSR1","LINGO1","MAP2K1","MAP2K2","MAP2K5","MAP3K1","MAP3K2","MAPK1","MAPK10","MAPK14","MAPK3","MAPK7","MAPK8","MAPK9","MAPT","MARCKS","MEF2A","MEF2C","MTOR","NCAM1","NCF1","NCF2","NCK1","NCK2","NFATC4","NFKB1","NFKBIA","NGF","NGFR","NSF","NTF3","NTRK1","NTRK2","NTRK3","PDPK1","PIK3R1","PIK3R2","PLCG1","PPP2CA","PRKAA1","PRKAA2","PRKCD","PTK2B","PTPN11","PTPRF","RAB3A","RAC1","RACK1","RAF1","RANBP9","RAP1A","RASGRF1","RELA","RHOG","RPS6","RPS6KA1","RPS6KA3","RPS6KA5","RPS6KB1","SH2B1","SH2B2","SHC1","SHC2","SHC3","SHC4","SIRPA","SORT1","SPP1","SQSTM1","SRC","STAT1","STAT3","STAT5A","STAT5B","SYN1","TIAM1","TRAF6","TSC2","VAV2","VAV3","YBX1"))
colnames(WP2380) <-("gene")

#Repetition
i <- 1
resultList = list()
repeat  {
  
  count_WP5412 <- c(dim(WP5412)[1])
  n_WP5412 <- data.frame(sample(Available_genes, size = count_WP5412, replace = FALSE, prob = NULL))
  colnames(n_WP5412) <-("gene")
  
  WP2380_WP5412 <- inner_join(WP2380,n_WP5412)
  result <- c(dim(WP2380_WP5412)[1])

  resultList[[length(resultList)+1]] = result
  i <- i+1
  if (i>100) {break}
  }

#repetition end - number of repetition is set to 100 for testing. 10 000 takes a while.

#Results are in resultList, as double list, this step does type conversion and extraction of results in form of mean, sd and histogram plot
resultList <- as.numeric(resultList)
hist(resultList)
mean(resultList)
sd(resultList)
