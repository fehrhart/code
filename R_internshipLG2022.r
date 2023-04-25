---
title: "Internship"
output: html_document
---

```{r}
#Select necessary packages
library(readxl)
library(dplyr)
library(Cairo)
```

```{r}
#Import all Homo sapiens genes
all_genes <- read_excel("Data files/All_genes.xlsx")
all_genes = all_genes %>% select(1,2,3) #Select filtered columns
colnames(all_genes) <- c("Ensembl ID","HGNC Symbol","EntrezID") #name columns
all_genes <- filter(all_genes, rowSums(is.na(all_genes)) != ncol(all_genes)) #remove rows containing only N/A value 
all_genes
```

```{r}
#Import the list of schizophrenia-associated genes obtained from the Psychiatric Genomics Consortium
Ensembl_export_list_schizophrenia <- read_excel("Data files/Ensembl export lists complete/Ensembl export list schizophrenia.xls")
schizophrenia_genes = Ensembl_export_list_schizophrenia %>% select(1,2,3) #Select filtered columns
colnames(schizophrenia_genes) <- c("Ensembl ID","HGNC Symbol","EntrezID") #name columns
schizophrenia_genes <- filter(schizophrenia_genes, rowSums(is.na(schizophrenia_genes)) != ncol(schizophrenia_genes)) #remove rows containing only N/A value 
schizophrenia_genes
```

```{r}
#Import the list of genes involved in brain development from GeneOntology
Ensembl_export_list_GeneOntology_brain_development <- read_excel("Data files/Ensembl export lists complete/Ensembl export list GeneOntology brain development.xls")
brain_development_genes = Ensembl_export_list_GeneOntology_brain_development %>% select(c(1,2,3,4)) #Select filtered columns
colnames(brain_development_genes) <- c("Ensembl ID","HGNC Symbol","EntrezID","UniProtKB gene name ID") #name columns
brain_development_genes <- filter(brain_development_genes, rowSums(is.na(brain_development_genes)) != ncol(brain_development_genes)) #remove rows containing only N/A value
brain_development_genes
```

```{r}
#Import the list of genes involved in cerebral cortex development from GeneOntology
Ensembl_export_list_GeneOntology_cerebral_cortex_development <- read_excel("Data files/Ensembl export lists complete/Ensembl export list GeneOntology cerebral cortex development.xls")
cc_development_genes = Ensembl_export_list_GeneOntology_cerebral_cortex_development %>% select(c(1,2,3,4)) #Select filtered columns
colnames(cc_development_genes) <- c("Ensembl ID","HGNC Symbol","EntrezID","UniProtKB gene name ID") #name columns
cc_development_genes <- filter(cc_development_genes, rowSums(is.na(cc_development_genes)) != ncol(cc_development_genes)) #remove rows containing only N/A value
cc_development_genes
```

```{r}
#Import the list of genes involved in cerebral cortex development from Grasby et al. (2021)

#Import genes associated with surface area
Ensembl_export_list_Grasby_SA <- read_excel("Data files/Ensembl export lists complete/Ensembl export list Grasby SA.xls")
Grasby_genes_SA = Ensembl_export_list_Grasby_SA %>% select(c(1,2,3)) #Select filtered columns
colnames(Grasby_genes_SA) <- c("Ensembl ID","HGNC Symbol","EntrezID") #Name columns
Grasby_genes_SA <- filter(Grasby_genes_SA, rowSums(is.na(Grasby_genes_SA)) != ncol(Grasby_genes_SA)) #remove rows containing only N/A value
Grasby_genes_SA

#Import genes associated with thickness
Grasby_genes_TH <- read_excel("Data files/Ensembl export lists complete/Ensembl export list Grasby TH.xls")
colnames(Grasby_genes_TH) <- c("Ensembl ID","HGNC Symbol","EntrezID") #Name columns
Grasby_genes_TH <- filter(Grasby_genes_TH, rowSums(is.na(Grasby_genes_TH)) != ncol(Grasby_genes_TH)) #remove rows containing only N/A value
Grasby_genes_TH
```

```{r}
#Extracting which genes are involved in both brain and cerebral cortex development (GeneOntology)
b_vs_cc = inner_join(brain_development_genes, cc_development_genes)
b_vs_cc
```

```{r}
#Extracting which genes are involved in brain development and have been linked to schizophrenia
b_vs_s = inner_join(brain_development_genes, schizophrenia_genes)
b_vs_s
```

```{r}
#Extracting which genes are involved in cerebral cortex development (GeneOntology) and have been linked to schizophrenia
cc_vs_s = inner_join(cc_development_genes,schizophrenia_genes)
cc_vs_s
```

```{r}
#Extracting which genes are involved in brain development and surface area (Grasby)
b_vs_sa = inner_join(brain_development_genes, Grasby_genes_SA)
b_vs_sa

install.packages("writexl")
library("writexl")
write_xlsx(b_vs_sa,"Data files/RStudio export files\\brain development vs surface area.xlsx")
```

```{r}
#Extracting which genes are involved in brain development and thickness (Grasby)
b_vs_th = inner_join(brain_development_genes, Grasby_genes_TH)
b_vs_th

write_xlsx(b_vs_th,"Data files/RStudio export files\\brain development vs thickness.xlsx")
```

```{r}
#Extracting which genes are involved in cerebral cortex development and surface area
cc_vs_sa = inner_join(cc_development_genes, Grasby_genes_SA)
cc_vs_sa

write_xlsx(cc_vs_sa,"Data files/RStudio export files\\cerebral cortex development vs surface area.xlsx")
```

```{r}
#Extracting which genes are involved in cerebral cortex development and thickness
cc_vs_th = inner_join(cc_development_genes, Grasby_genes_TH)
cc_vs_th

write_xlsx(cc_vs_th,"Data files/RStudio export files\\cerebral cortex development vs thickness.xlsx")
```

```{r}
#Extracting which genes are involved in surface area (Grasby) and have been linked to schizophrenia
s_vs_sa = inner_join(schizophrenia_genes, Grasby_genes_SA)
s_vs_sa
```

```{r}
#Extracting which genes are involved in surface area (Grasby) and have been linked to schizophrenia
s_vs_th = inner_join(schizophrenia_genes, Grasby_genes_TH)
s_vs_th
```

```{r}
#Extracting which genes are involved in surface area and thickness (Grasby)
sa_vs_th = inner_join(Grasby_genes_SA,Grasby_genes_TH)
sa_vs_th
```

```{r}
#Extracting which genes are involved in brain development, cerebral cortex development, and schizophrenia
b_vs_cc_vs_s = inner_join(b_vs_cc, schizophrenia_genes)
b_vs_cc_vs_s
```

```{r}
#Extracting which genes are involved in brain development, cerebral cortex development, and surface area
b_vs_cc_vs_sa = inner_join(b_vs_cc,Grasby_genes_SA)
b_vs_cc_vs_sa
```

```{r}
#Extracting which genes are involved in brain development, cerebral cortex development, and thickness
b_vs_cc_vs_th = inner_join(b_vs_cc,Grasby_genes_TH)
b_vs_cc_vs_th
```

```{r}
#Extracting which genes are involved in brain development, schizophrenia, and surface area
b_vs_s_vs_sa = inner_join(b_vs_s,Grasby_genes_SA)
b_vs_s_vs_sa
```

```{r}
#Extracting which genes are involved in brain development, schizophrenia, and thickness
b_vs_s_vs_th = inner_join(b_vs_s,Grasby_genes_TH)
b_vs_s_vs_th
```

```{r}
#Extracting which genes are involved in brain development, surface area, and surface area
b_vs_sa_vs_th = inner_join(b_vs_sa,Grasby_genes_TH)
b_vs_sa_vs_th
```

```{r}
#Extracting which genes are involved in cerebral cortex development, schizophrenia, and surface area
cc_vs_s_vs_sa = inner_join(cc_vs_s,Grasby_genes_SA)
cc_vs_s_vs_sa
```

```{r}
#Extracting which genes are involved in cerebral cortex development, schizophrenia, and thickness
cc_vs_s_vs_th = inner_join(cc_vs_s,Grasby_genes_TH)
cc_vs_s_vs_th
```

```{r}
#Extracting which genes are involved in cerebral cortex development, surface area, and thickness
cc_vs_sa_vs_th = inner_join(cc_vs_sa,Grasby_genes_TH)
cc_vs_sa_vs_th
```

```{r}
#Extracting which genes are involved in schizophrenia, surface area, and thickness
s_vs_sa_vs_th = inner_join(s_vs_sa,Grasby_genes_TH)
s_vs_sa_vs_th
```

```{r}
b_vs_cc_vs_s_vs_sa = inner_join(b_vs_cc, s_vs_sa)
b_vs_cc_vs_s_vs_sa
```

```{r}
b_vs_s_vs_sa_vs_th = inner_join(b_vs_s,sa_vs_th)
b_vs_s_vs_sa_vs_th
```

```{r}
#download RCy3 package from BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RCy3")

library(RCy3) #add RCy3 to local library
```

```{r}
#check to see if RStudio is connected to Cytoscape
cytoscapePing ()
cytoscapeVersionInfo ()
```

Code for nodes, edges, and network creation obtained from https://kozo2.github.io/cyautoworkshop/articles/cyautoworkshop.html (23-5-2022).

```{r}
#Add row to include schizophrenia node
schizophrenia_cytoscape = schizophrenia_genes[,1:2] %>% add_row(`Ensembl ID` = "Schizophrenia", `HGNC Symbol` = "Schizophrenia")

#create nodes
nodes_schizophrenia <- data.frame(id = c(schizophrenia_cytoscape[,2]),
    group = c("Network"),
    score = as.integer(1),
    BP = "schizophrenia",
    stringsAsFactors=FALSE)
colnames(nodes_schizophrenia) <- c("id","group","score")
nodes_schizophrenia

#create edges
edges_schizophrenia <- data.frame(source=c(schizophrenia_cytoscape[,2]),
    target = "Schizophrenia",
    interaction = "Associated with",
    weight = c(1),
    stringsAsFactors=FALSE)
colnames(edges_schizophrenia) <- c("source","target","interaction","weight")
edges_schizophrenia

#create Network
createNetworkFromDataFrames(nodes_schizophrenia, edges_schizophrenia, title="Schizophrenia", collection="Internship Network")
```

```{r}
#Add row to include Grasby SA node
Grasby_SA_cytoscape = Grasby_genes_SA[,1:2] %>% add_row(`Ensembl ID` = "Grasby_SA", `HGNC Symbol` = "Grasby_SA")

#create nodes
nodes_grasby_SA <- data.frame(id = c(Grasby_SA_cytoscape[,2]),
    group = c("Network"),
    score = as.integer(1),
    BP = "Surface Area",
    stringsAsFactors=FALSE)
colnames(nodes_grasby_SA) <- c("id","group","score")
nodes_grasby_SA

#create edges
edges_grasby_SA <- data.frame(source=c(Grasby_SA_cytoscape[,2]),
    target = "Grasby_SA",
    interaction = "Associated with",
    weight = c(1),
    stringsAsFactors=FALSE)
colnames(edges_grasby_SA) <- c("source","target","interaction","weight")
edges_grasby_SA

#create Network
createNetworkFromDataFrames(nodes_grasby_SA, edges_grasby_SA, title="Grasby_SA", collection="Internship Network")
```

```{r}
#Add row to include Grasby SA node
Grasby_TH_cytoscape = Grasby_genes_TH[,1:2] %>% add_row(`Ensembl ID` = "Grasby_TH", `HGNC Symbol` = "Grasby_TH")

#create nodes
nodes_grasby_TH <- data.frame(id = c(Grasby_TH_cytoscape[,2]),
    group = c("Network"),
    score = as.integer(1),
    BP = "Thickness",
    stringsAsFactors=FALSE)
colnames(nodes_grasby_TH) <- c("id","group","score")
nodes_grasby_TH

#create edges
edges_grasby_TH <- data.frame(source=c(Grasby_TH_cytoscape[,2]),
    target = "Grasby_TH",
    interaction = "Associated with",
    weight = c(1),
    stringsAsFactors=FALSE)
colnames(edges_grasby_TH) <- c("source","target","interaction","weight")
edges_grasby_TH

#create Network
createNetworkFromDataFrames(nodes_grasby_TH, edges_grasby_TH, title="Grasby_TH", collection="Internship Network")
```

```{r}
#Add row to include schizophrenia node
brain_cytoscape = brain_development_genes[,1:2] %>% add_row(`Ensembl ID` = "Brain development", `HGNC Symbol` = "Brain development")

#create nodes
nodes_brain <- data.frame(id = c(brain_cytoscape[,2]),
    group = c("Network"),
    score = as.integer(1),
    BP = "Brain development",
    stringsAsFactors=FALSE)
colnames(nodes_brain) <- c("id","group","score")
nodes_brain

#create edges
edges_brain <- data.frame(source=c(brain_cytoscape[,2]),
    target = "Brain development",
    interaction = "Associated with",
    weight = c(1),
    stringsAsFactors=FALSE)
colnames(edges_brain) <- c("source","target","interaction","weight")
edges_brain

#create Network
createNetworkFromDataFrames(nodes_brain, edges_brain, title="Brain development", collection="Internship Network")
```

```{r}
#Add row to include schizophrenia node
cerebral_cytoscape = cc_development_genes[,1:2] %>% add_row(`Ensembl ID` = "Cerebral development", `HGNC Symbol` = "Cerebral development")

#create nodes
nodes_cerebral <- data.frame(id = c(cerebral_cytoscape[,2]),
    group = c("Network"),
    score = as.integer(1),
    BP = "Cerebral cortex development",
    stringsAsFactors=FALSE)
colnames(nodes_cerebral) <- c("id","group","score")
nodes_cerebral

#create edges
edges_cerebral <- data.frame(source=c(cerebral_cytoscape[,2]),
    target = "Cerebral development",
    interaction = "Associated with",
    weight = c(1),
    stringsAsFactors=FALSE)
colnames(edges_cerebral) <- c("source","target","interaction","weight")
edges_cerebral

#create Network
createNetworkFromDataFrames(nodes_cerebral, edges_cerebral, title="Cerebral cortex development", collection="Internship Network")
```

Test run for network creation in Cytoscape
```{r}
nodes_test <- data.frame(id=c("node 0","node 1","node 2","node 3"),
    group=c("A","A","B","B"), # categorical strings
    score=as.integer(c(20,10,15,5)), # integers
    stringsAsFactors=FALSE)
nodes_test

edges_test <- data.frame(source=c("node 0","node 0","node 0","node 2"),
    target=c("node 1","node 2","node 3","node 3"),
    interaction=c("inhibits","interacts","activates","interacts"),  # optional
    weight=c(5.1,3.0,5.2,9.9), # numeric
    stringsAsFactors=FALSE)
edges_test

createNetworkFromDataFrames(nodes_test, edges_test, title="Test", collection="Test")
```

```{r}
#install rWikiPathways
if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly = T))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways", update = F)
}
library(rWikiPathways)
```

```{r}
#Install necessary packages according to BigCat Pathway Analysis vignette
load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
    print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
    cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
    status
}
```

According to the vignette, we need to prepare the data set (e.g., up- and down-regulated genes). Due to the nature of the used data sets, this is not possible.

```{r}
#import node table from Cytoscape of complete network merge to isolate all genes associated with Schizophrenia and at least one other data set.
library(tidyverse)

#Transform all_genes data set into vector
genes_ensembl <- all_genes %>% pull("Ensembl ID")
genes_ensembl

#Transform Ensembl ID to Entrez ID
all_genes_entrez  <- clusterProfiler::bitr(genes_ensembl,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
all_genes_entrez

#Data frame of all genes associated with Schizophrenia
merge_schizophrenia <- read.csv("Data files/Cytoscape export files/Schizophrenia genes associated with biological processes.csv")
merge_schizophrenia

#Enrichment analysis -> Gene Ontology
egobp <- clusterProfiler::enrichGO(
        gene     = merge_schizophrenia[[2]],
        universe = all_genes_entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

head(egobp,10)

as.data.frame(egobp)

#Make and save plots
barplot1 <- barplot(egobp, showCategory = 20)
ggsave(barplot1, path = "Results/Enrichment analysis/Background = human genome/Schizophrenia + miscellaneous", filename = "GO barplot - Schizophrenia + miscellaneous.png", dpi = 500, type = "cairo")

dotplot1 <- dotplot(egobp, showCategory = 20)
ggsave(dotplot1, path = "Results/Enrichment analysis/Background = human genome/Schizophrenia + miscellaneous", filename = "GO dotplot - Schizophrenia + miscellaneous.png", dpi = 500,  type = "cairo")

goplot1 <- goplot(egobp)
ggsave(goplot1, path = "Results/Enrichment analysis/Background = human genome/Schizophrenia + miscellaneous", filename = "GO goplot - Schizophrenia + miscellaneous.png", dpi = 500, type = "cairo")


## extract a data frame with results from object of type enrichResult
egobp.results.df <- egobp@result

## create a new column for term size from BgRatio
egobp.results.df$term.size <- gsub("/(\\d+)", "", egobp.results.df$BgRatio)

## filter for term size to keep only term.size => 3, gene count >= 5 and subset
egobp.results.df <- egobp.results.df[which(egobp.results.df[,'term.size'] >= 3 & egobp.results.df[,'Count'] >= 5),]
egobp.results.df <- egobp.results.df[c("ID", "Description", "pvalue", "qvalue", "geneID")]

## format gene list column
egobp.results.df$geneID <- gsub("/", ",", egobp.results.df$geneID)

## add column for phenotype
egobp.results.df <- cbind(egobp.results.df, phenotype=1)
egobp.results.df <- egobp.results.df[, c(1, 2, 3, 4, 6, 5)]

## change column headers
colnames(egobp.results.df) <- c("Name","Description", "pvalue","qvalue","phenotype", "genes")

egobp.results.filename <-file.path(getwd(),paste("clusterprofiler_cluster_enr_results.txt",sep="_"))
write.table(egobp.results.df,egobp.results.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

em_command = paste('enrichmentmap build analysisType="generic" ', 
                   'pvalue=',"0.05", 'qvalue=',"0.05",
                   'similaritycutoff=',"0.25",
                   'coeffecients=',"JACCARD",
                   'enrichmentsDataset1=',egobp.results.filename ,
                   sep=" ")

  #enrichment map command will return the suid of newly created network.
  em_network_suid <- commandsRun(em_command)
  
  renameNetwork("Cluster1_enrichmentmap", network=as.numeric(em_network_suid))
```

Test run for EnrichmentMap in Cytoscape
```{r}
lung.expr <- read.csv(system.file("extdata","data-lung-cancer.csv", package="rWikiPathways"),stringsAsFactors = FALSE)
nrow(lung.expr)
head(lung.expr)

up.genes <- lung.expr[lung.expr$log2FC > 1 & lung.expr$adj.P.Value < 0.05, 1] 
dn.genes <- lung.expr[lung.expr$log2FC < -1 & lung.expr$adj.P.Value < 0.05, 1]
bkgd.genes <- lung.expr[,1]
bkgd.genes

up.genes.entrez <- clusterProfiler::bitr(up.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
up.genes.entrez
bkgd.genes.entrez <- bitr(bkgd.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez

egobp_test <- clusterProfiler::enrichGO(
        gene     = up.genes.entrez[[2]],
        universe = bkgd.genes.entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

head(egobp_test,10)

barplot(egobp_test, showCategory = 20)
dotplot(egobp_test, showCategory = 20)
goplot(egobp_test)

## extract a data frame with results from object of type enrichResult
egobp.test.results.df <- egobp_test@result

## create a new column for term size from BgRatio
egobp.test.results.df$term.size <- gsub("/(\\d+)", "", egobp.test.results.df$BgRatio)

## filter for term size to keep only term.size => 3, gene count >= 5 and subset
egobp.test.results.df <- egobp.test.results.df[which(egobp.test.results.df[,'term.size'] >= 3 & egobp.test.results.df[,'Count'] >= 5),]
egobp.test.results.df <- egobp.test.results.df[c("ID", "Description", "pvalue", "qvalue", "geneID")]

## format gene list column
egobp.test.results.df$geneID <- gsub("/", ",", egobp.test.results.df$geneID)

## add column for phenotype
egobp.test.results.df <- cbind(egobp.test.results.df, phenotype=1)
egobp.test.results.df <- egobp.test.results.df[, c(1, 2, 3, 4, 6, 5)]

## change column headers
colnames(egobp.test.results.df) <- c("Name","Description", "pvalue","qvalue","phenotype", "genes")

egobp.test.results.filename <-file.path(getwd(),paste("clusterprofiler_cluster_enr_results_test.txt",sep="_"))
write.table(egobp.test.results.df,egobp.test.results.filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

em_command = paste('enrichmentmap build analysisType="generic" ', 
                   'pvalue=',"0.05", 'qvalue=',"0.05",
                   'similaritycutoff=',"0.25",
                   'coeffecients=',"JACCARD",
                   'enrichmentsDataset1=',egobp.test.results.filename ,
                   sep=" ")

  #enrichment map command will return the suid of newly created network.
  em_network_suid <- commandsRun(em_command)
  
  renameNetwork("Cluster1_enrichmentmap_test", network=as.numeric(em_network_suid))
  

ewp.up <- clusterProfiler::enrichWP(
    up.genes.entrez[[2]],
    universe = bkgd.genes.entrez[[2]],
    organism = "Homo sapiens",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.1, #p.adjust cutoff; relaxed for demo purposes
)

head(ewp.up)

ewp.up <- DOSE::setReadable(ewp.up, org.Hs.eg.db, keyType = "ENTREZID")
head(ewp.up)

barplot(ewp.up, showCategory = 20)
dotplot(ewp.up, showCategory = 20)


ewp.dn <- enrichWP(
    dn.genes.entrez[[2]],
    #universe = bkgd.genes[[2]],  #hint: comment out to get any results for demo
    organism = "Homo sapiens",
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.1, #p.adjust cutoff; relaxed for demo purposes
)

 ewp.dn <- setReadable(ewp.dn, org.Hs.eg.db, keyType = "ENTREZID")
 head(ewp.dn)
 dotplot(ewp.dn, showCategory = 20)
```

```{r}
ewp <- clusterProfiler::enrichWP(
  merge_schizophrenia[[2]],
  universe = all_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.1, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

head(ewp)

ewp <- DOSE::setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID")
head(ewp)

barplot2 <- barplot(ewp, showCategory = 20)
ggsave(barplot2, path = "Results/Enrichment analysis/Background = human genome/Schizophrenia + miscellaneous", filename = "WP barplot - Schizophrenia + miscellaneous.png", dpi = 500, type = "cairo")

dotplot2 <- dotplot(ewp, showCategory = 20)
ggsave(dotplot2, path = "Results/Enrichment analysis/Background = human genome/Schizophrenia + miscellaneous", filename = "WP dotplot - Schizophrenia + miscellaneous.png", dpi = 500, type = "cairo")
```

```{r Enrichment analysis of surface area, thickness, and schizophrenia associated genes w/ universe = surface area genes}
#Enrichment analysis w/ gene = TH, SA, and S
#Extract gene list
s_sa_th_genes <- read.csv("Data files/Cytoscape export files/SA, TH, & S genes.csv")
s_sa_th_genes

#Transform Grasby_genes_SA data set into vector
sa_genes <- Grasby_genes_SA %>% pull("Ensembl ID")
sa_genes

#Transform Ensembl ID to Entrez ID
sa_genes_entrez  <- clusterProfiler::bitr(sa_genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
sa_genes_entrez

#Enrichment analysis GeneOntology
egosa <- clusterProfiler::enrichGO(
        gene     = s_sa_th_genes[[9]],
        universe = sa_genes_entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.1, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

head(egosa, 10)

#Enrichment analysis WikiPathways
ewpsa <- clusterProfiler::enrichWP(
  s_sa_th_genes[[9]],
  universe = sa_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.1, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

head(ewpsa)
```

```{r Enrichment analysis of surface area genes that overlap with schizophrenia genes w/ universe = all_genes}
#Import gene list
s_sa_genes <- read.csv("Data files/Cytoscape export files/SA & S genes.csv")
s_sa_genes <- s_sa_genes %>% pull("Ensembl")
sa_genes

#Transform Ensembl ID to Entrez ID
s_sa_genes_entrez  <- clusterProfiler::bitr(s_sa_genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
s_sa_genes_entrez

#GeneOntology enrichment analysis
ego_s_sa <- clusterProfiler::enrichGO(
        gene     = s_sa_genes_entrez[[2]],
        universe = all_genes_entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

head(ego_s_sa, 10)


#WikiPathways enrichment analysis
ewp_s_sa <- clusterProfiler::enrichWP(
  s_sa_genes_entrez[[2]],
  universe = all_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

head(ewp_s_sa)

as.data.frame(ego_s_sa)
as.data.frame(ewp_s_sa)

#Making and saving plots
barplot3 <- barplot(ewp_s_sa, showCategory = 20)
ggsave(barplot3, path = "Results/Enrichment analysis/Background = human genome/Schizophrenia & SA", filename = "WP barplot - Schizophrenia & SA.png", dpi = 500, type = "cairo")

dotplot3 <- dotplot(ewp_s_sa, showCategory = 20)
ggsave(dotplot3, path = "Results/Enrichment analysis/Background = human genome/Schizophrenia & SA", filename = "WP dotplot - Schizophrenia & SA.png", dpi = 500, type = "cairo")
```

```{r Enrichment analysis of total surface area genes that overlap with schizophrenia genes w/ universe = surface area genes}
#Import gene list
s_tsa_genes <- read_excel("Overlap schizophrenia genes and brain areas.xlsx")
s_tsa_genes <- s_tsa_genes %>% select(6) %>% rename("SYMBOL" = "Genes associated with schizophrenia and total surface area")
s_tsa_genes <- na.omit(s_tsa_genes)
s_tsa_genes <- s_tsa_genes %>% pull("SYMBOL")
s_tsa_genes

#Transform Ensembl ID to Entrez ID
s_tsa_genes_entrez  <- clusterProfiler::bitr(s_tsa_genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
s_tsa_genes_entrez

#GeneOntology enrichment analysis
ego_s_tsa_usa <- clusterProfiler::enrichGO(
        gene     = s_tsa_genes_entrez[[2]],
        universe = sa_genes_entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

as.data.frame(ego_s_tsa_usa)

#WikiPathways enrichment analysis
ewp_s_tsa_usa <- clusterProfiler::enrichWP(
  s_tsa_genes_entrez[[2]],
  universe = sa_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

head(ewp_s_tsa_usa)
```

```{r Enrichment analysis of total surface area genes that overlap with schizophrenia genes w/ universe = all_genes}
#Import gene list
s_tsa_genes <- read_excel("Overlap schizophrenia genes and brain areas.xlsx")
s_tsa_genes <- s_tsa_genes %>% select(6) %>% rename("SYMBOL" = "Genes associated with schizophrenia and total surface area")
s_tsa_genes <- na.omit(s_tsa_genes)
s_tsa_genes <- s_tsa_genes %>% pull("SYMBOL")
s_tsa_genes

#Transform Ensembl ID to Entrez ID
s_tsa_genes_entrez  <- clusterProfiler::bitr(s_tsa_genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
s_tsa_genes_entrez

#GeneOntology enrichment analysis
ego_s_tsa <- clusterProfiler::enrichGO(
        gene     = s_tsa_genes_entrez[[2]],
        universe = all_genes_entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.2, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

as.data.frame(ego_s_tsa)

#WikiPathways enrichment analysis
ewp_s_tsa <- clusterProfiler::enrichWP(
  s_tsa_genes_entrez[[2]],
  universe = all_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

head(ewp_s_tsa)
```

```{r Enrichment analysis of fundamental gene sets -> brain development}
#Transform brain development gene data frame
brain_development_genes
brain_development_genes_ensembl <- brain_development_genes %>% pull("Ensembl ID")
brain_development_genes_entrez <- clusterProfiler::bitr(brain_development_genes_ensembl,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

#Enrichment analysis of Gerontology vs human genome
ego_bd <- clusterProfiler::enrichGO(
        gene     = brain_development_genes_entrez [[2]],
        universe = all_genes_entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

as.data.frame(ego_bd)

ewp_bd <- clusterProfiler::enrichWP(
  brain_development_genes_entrez[[2]],
  universe = all_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

as.data.frame(ewp_bd)

#Plots
barplot4 <- barplot(ego_bd, showCategory = 20)
ggsave(barplot4, path = "Results/Enrichment analysis/Background = human genome/Brain development", filename = "GO barplot - Brain development.png", dpi = 500, type = "cairo")

dotplot4 <- dotplot(ego_bd, showCategory = 20)
ggsave(dotplot4, path = "Results/Enrichment analysis/Background = human genome/Brain development", filename = "GO dotplot - Brain development.png", dpi = 500, type = "cairo")

goplot2 <- goplot(ego_bd, showCategory = 20)
ggsave(goplot2, path = "Results/Enrichment analysis/Background = human genome/Brain development", filename = "GO goplot - Brain development.png", dpi = 500, type = "cairo")

barplot5 <- barplot(ewp_bd, showCategory = 20)
ggsave(barplot5, path = "Results/Enrichment analysis/Background = human genome/Brain development", filename = "WP barplot - Brain development.png", dpi = 700, type = "cairo", width = 10, height = 10)

dotplot5 <- dotplot(ewp_bd, showCategory = 20)
ggsave(dotplot5, path = "Results/Enrichment analysis/Background = human genome/Brain development", filename = "WP dotplot - Brain development.png", dpi = 500, type = "cairo")
```

```{r Enrichment analysis of fundamental gene sets -> cerebral cortex development}
#Transform brain development gene data frame
cc_development_genes
cc_development_genes_ensembl <- cc_development_genes %>% pull("Ensembl ID")
cc_development_genes_entrez <- clusterProfiler::bitr(cc_development_genes_ensembl,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

#Enrichment analysis of Gerontology vs human genome
ego_ccd <- clusterProfiler::enrichGO(
        gene     = cc_development_genes_entrez [[2]],
        universe = all_genes_entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

as.data.frame(ego_ccd)

ewp_ccd <- clusterProfiler::enrichWP(
  cc_development_genes_entrez[[2]],
  universe = all_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

as.data.frame(ewp_ccd)

#Making and saving plots
barplot6 <- barplot(ego_ccd, showCategory = 20)
ggsave(barplot6, path = "Results/Enrichment analysis/Background = human genome/Cerebral cortex development", filename = "GO barplot - Cerebral cortex development genes.png", dpi = 500, type = "cairo")

dotplot6 <- dotplot(ego_ccd, showCategory = 20)
ggsave(dotplot6, path = "Results/Enrichment analysis/Background = human genome/Cerebral cortex development", filename = "GO dotplot - Cerebral cortex development genes.png", dpi = 500, type = "cairo")

goplot3 <- goplot(ego_ccd)
ggsave(goplot3, path = "Results/Enrichment analysis/Background = human genome/Cerebral cortex development", filename = "GO goplot - Cerebral cortex development genes.png", dpi = 500, type = "cairo")

barplot7 <- barplot(ewp_ccd, showCategory = 20)
ggsave(barplot7, path = "Results/Enrichment analysis/Background = human genome/Cerebral cortex development", filename = "WP barplot - Cerebral cortex development genes.png", dpi = 500, type = "cairo", width = 12, height = 12)

dotplot7 <- dotplot(ewp_ccd, showCategory = 20)
ggsave(dotplot7, path = "Results/Enrichment analysis/Background = human genome/Cerebral cortex development", filename = "WP dotplot - Cerebral cortex development genes.png", dpi = 500, type = "cairo")
```

```{r Enrichment analysis of fundamental gene sets -> schizophrenia associated genes}
#Transform brain development gene data frame
schizophrenia_genes
schizophrenia_genes_ensembl <- schizophrenia_genes %>% pull("Ensembl ID")
schizophrenia_genes_entrez <- clusterProfiler::bitr(schizophrenia_genes_ensembl,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

#Enrichment analysis of Gerontology vs human genome
ego_s <- clusterProfiler::enrichGO(
        gene     = schizophrenia_genes_entrez [[2]],
        universe = all_genes_entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

as.data.frame(ego_s)

ewp_s <- clusterProfiler::enrichWP(
  schizophrenia_genes_entrez[[2]],
  universe = all_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

as.data.frame(ewp_s)

#Plots
barplot(ewp_s, showCategory = 20)
dotplot(ewp_s, showCategory = 20)
```

```{r Enrichment analysis of fundamental gene sets -> cortical surface area (Grasby)}
#Transform brain development gene data frame
Grasby_genes_SA
Grasby_genes_SA_ensembl <- Grasby_genes_SA %>% pull("Ensembl ID")
Grasby_genes_SA_entrez <- clusterProfiler::bitr(Grasby_genes_SA_ensembl,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

#Enrichment analysis of Gerontology vs human genome
ego_sa <- clusterProfiler::enrichGO(
        gene     = Grasby_genes_SA_entrez [[2]],
        universe = all_genes_entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

as.data.frame(ego_sa)

ewp_sa <- clusterProfiler::enrichWP(
  Grasby_genes_SA_entrez[[2]],
  universe = all_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

as.data.frame(ewp_sa)
```

```{r Enrichment analysis of fundamental gene sets -> cortical surface area (Grasby)}
#Transform brain development gene data frame
Grasby_genes_TH
Grasby_genes_TH_ensembl <- Grasby_genes_TH %>% pull("Ensembl ID")
Grasby_genes_TH_entrez <- clusterProfiler::bitr(Grasby_genes_TH_ensembl,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

#Enrichment analysis of Gerontology vs human genome
ego_th <- clusterProfiler::enrichGO(
        gene     = Grasby_genes_TH_entrez [[2]],
        universe = all_genes_entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

as.data.frame(ego_th)

ewp_th <- clusterProfiler::enrichWP(
  Grasby_genes_TH_entrez[[2]],
  universe = all_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

as.data.frame(ewp_th)
```

```{r Enrichment analysis genes = cerebral cortex development, background = brain development}

#Enrichment analysis of Gerontology vs human genome
ego_cc_b <- clusterProfiler::enrichGO(
        gene     = cc_development_genes_entrez [[2]],
        universe = brain_development_genes_entrez [[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

as.data.frame(ego_cc_b)

ewp_cc_b <- clusterProfiler::enrichWP(
  cc_development_genes_entrez[[2]],
  universe = brain_development_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.2, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

as.data.frame(ewp_cc_b)

#Plots
barplot(ego_cc_b, showCategory = 20)
dotplot(ego_cc_b, showCategory = 20)
goplot(ego_cc_b)
```

```{r Enrichment analysis genes = surface area, background = schizophrenia}

#Enrichment analysis of Geneontology vs human genome
ego_sa_s <- clusterProfiler::enrichGO(
        gene     = Grasby_genes_SA_entrez [[2]],
        universe = schizophrenia_genes_entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

as.data.frame(ego_sa_s)

ewp_sa_s <- clusterProfiler::enrichWP(
  Grasby_genes_SA_entrez[[2]],
  universe = schizophrenia_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

as.data.frame(ewp_sa_s)
```

```{r Enrichment analysis genes = thickness, background = schizophrenia}

#Enrichment analysis of Geneontology vs human genome
ego_th_s <- clusterProfiler::enrichGO(
        gene     = Grasby_genes_TH_entrez [[2]],
        universe = schizophrenia_genes_entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

as.data.frame(ego_th_s)

ewp_th_s <- clusterProfiler::enrichWP(
  Grasby_genes_TH_entrez[[2]],
  universe = schizophrenia_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

as.data.frame(ewp_th_s)
```

```{r Enrichment analysis genes = surface area, background = brain development}

#Enrichment analysis of Geneontology vs human genome
ego_sa_b <- clusterProfiler::enrichGO(
        gene     = Grasby_genes_SA_entrez [[2]],
        universe = brain_development_genes_entrez [[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.2, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

as.data.frame(ego_sa_b)

ewp_sa_b <- clusterProfiler::enrichWP(
 Grasby_genes_SA_entrez[[2]],
  universe = brain_development_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.2, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

as.data.frame(ewp_sa_b)
```

```{r Enrichment analysis genes = thickness, background = brain development}

#Enrichment analysis of Gerontology vs human genome
ego_th_b <- clusterProfiler::enrichGO(
        gene     = Grasby_genes_TH_entrez [[2]],
        universe = brain_development_genes_entrez [[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.2, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

as.data.frame(ego_th_b)

ewp_th_b <- clusterProfiler::enrichWP(
  Grasby_genes_TH_entrez[[2]],
  universe = brain_development_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.2, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

as.data.frame(ewp_th_b)
```

```{r Enrichment analysis genes = cerebral cortex development, background = schizophrenia}

#Enrichment analysis of Gerontology vs human genome
ego_cc_s <- clusterProfiler::enrichGO(
        gene     = cc_development_genes_entrez [[2]],
        universe = schizophrenia_genes_entrez[[2]],
        OrgDb    = org.Hs.eg.db,
        ont      = "BP",
        pAdjustMethod = "fdr",
        pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
        readable = TRUE)

as.data.frame(ego_cc_s)

ewp_cc_s <- clusterProfiler::enrichWP(
  cc_development_genes_entrez[[2]],
  universe = schizophrenia_genes_entrez[[2]],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.5, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
)

as.data.frame(ewp_cc_s)

#Plots
barplot(ego_cc_s, showCategory = 20)
dotplot(ego_cc_s, showCategory = 20)
goplot(ego_cc_s)
```
