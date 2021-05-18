#clean workspace
rm(list=ls())

#install required packages and load them
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("Glimma")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("gplots")
BiocManager::install("RColorBrewer")

library(limma)
library(edgeR)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)

#set working directory - path to the folder where stuff is and should be stored
setwd("C:/Users/friederike.ehrhart/Documents/Teaching/Bachelor thesis/2021/2021 (1) Natalie Hinkova")

#Import table file
DF <- read.table("GSE119501_data.txt", header= TRUE, sep = "",fill = TRUE)  
as.data.frame.matrix(DF) 

#assign the first column with gene IDs as rownames
DF2 <- DF[,-1]
rownames(DF2) <- DF[,1]

#check for data and data structure
head (DF2) 
str(DF2)

#create design matrix by assigning experiment and control columns
design<-cbind(XPC=c(1,1,1,0,0,0),Con=c(0,0,0,1,1,1))
design

#voom transformation of data - make a "voom object" of the data and design matrix
par(mfrow=c(1,1))
v <- voom(DF2,design,plot = TRUE)
head(v)

#fit linear model and plot the names of the columns 
fit <- lmFit(v)
names(fit)

#define comparision XPC vs Con for statistics
cont.matrix <- makeContrasts(B.XPCVsCon=XPC - Con,levels=design)
cont.matrix

#calculate
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)

#print summary
summa.fit <- decideTests(fit.cont)
summary(summa.fit)

#plots for checking if everything worked well
par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=summa.fit[,"B.XPCVsCon"], values = c(-1, 1), hl.col=c("blue","red"))
volcanoplot(fit.cont,coef=1,highlight=100,main="B.XPCVsCon")

#export to txt
write.table(fit.cont, file="fit.tsv")