# This script follows the steps on this tutorial https://cytoscape.org/cytoscape-tutorials/protocols/wikipathways-app/#/ 
#clean workspace
rm(list=ls())

#install required packages and load them
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RCy3")
BiocManager::install("rWikiPathways")

library(RCy3)
library(rWikiPathways) 


#Check if Cytoscape is connected and active
cytoscapePing()
cytoscapeVersionInfo()

#search for a specific keyword - here "Dravet syndrome"
pathways <- findPathwaysByText('Dravet')
#Filter for human pathways
human.pathways <- subset(pathways, species == "Homo sapiens")
