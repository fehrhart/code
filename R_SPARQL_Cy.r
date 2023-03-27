#clean workspace
rm(list=ls())

#install required packages and load them
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# SPARQL package info: https://cran.r-project.org/web/packages/SPARQL/SPARQL.pdf
# (this version of) SPARQL package requires XML and RCurl
BiocManager::install("XML")
BiocManager::install("RCurl")
BiocManager::install("SPARQL")
BiocManager::install("RCy3")

library(XML)
library(RCurl)
library(SPARQL)
library(RCy3)


#define the SPARQL endpoint 
endpoint <- "https://sparql.wikipathways.org/sparql"
options <- NULL

#define the query itself
query <- ("

select distinct  (str(?title) as ?PathwayName) (str(?wpid) as ?PathwayID) (fn:substring(?genename,36) as ?GeneName) (fn:substring(?ncbiGeneId,33) as ?GeneID) where {
  ?gene a wp:DataNode ;
    dcterms:identifier ?id ;
    dcterms:isPartOf ?pathwayRes ;
    wp:bdbEntrezGene ?ncbiGeneId ;
    wp:bdbHgncSymbol ?genename .
  ?pathwayRes a wp:Pathway ;
    dcterms:identifier ?wpid ;
    dc:title ?title .
}
")

#execute query (will take a while to get 70 000 datapoints ;)
res <- SPARQL(endpoint,query)$results

#remove the / from the dataframe while keeping the structure
res[] <- lapply(res, gsub, pattern='/', replacement='')

#export results to file
write.table(res, file="sparql_res.tsv")

#arrange data for cytoscape import: 1. one column with all nodes
NodePathways <- res[,1:2]
NodeGenes <- res[,3:4]
names(NodeGenes) <- NULL
names(NodePathways) <- NULL
colnames(NodeGenes)<-c("Name", "ID")
colnames(NodePathways)<-c("Name", "ID")

Name_ID <- rbind(NodePathways, NodeGenes)

#make sure Cytoscape is up
cytoscapePing()

# Create network 
nodes <- data.frame(id=c(Name_ID[,2]),
                    group=c(Name_ID[,1]),
                    stringsAsFactors=FALSE)
edges <- data.frame(source=c(res[,4]),
                    target=c(res[,2]),
                    stringsAsFactors=FALSE)
createNetworkFromDataFrames(nodes, edges)



#Example from documentation
#nodes <- data.frame(id=c("node 0","node 1","node 2","node 3"),
                    #group=c("A","A","B","B"), #categorical strings
                   # score=as.integer(c(20,10,15,5)), #integers
                    #stringsAsFactors=FALSE)
#edges <- data.frame(source=c("node 0","node 0","node 0","node 2"),
                   # target=c("node 1","node 2","node 3","node 3"),
                   # interaction=c("inhibits","interacts",
                                  "activates","interacts"),  #optional
                  #  weight=c(5.1,3.0,5.2,9.9), #numerics
                  #  stringsAsFactors=FALSE)
#createNetworkFromDataFrames(nodes, edges)


