# This script is to do WGCNA analysis with the tutorial data from https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
#clean workspace
rm(list=ls())

#install required packages and load them
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("WGCNA")

library(WGCNA)
options(stringsAsFactors = FALSE)

#working directory
setwd("C:/Users/friederike.ehrhart/Documents/Streamline/RNAseq stay 25.9. - 1.12.23/WGCNA example")

#Read in the female liver data set
femData = read.csv("LiverFemale3600.csv")
# Read in the male liver data set
maleData = read.csv("LiverMale3600.csv")

# Take a quick look at what is in the data sets (caution, longish output):
dim(femData)
names(femData)
dim(maleData)
names(maleData)

# We work with two sets:
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Female liver", "Male liver")
shortLabels = c("Female", "Male")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(t(femData[-c(1:8)])));
names(multiExpr[[1]]$data) = femData$substanceBXH;
rownames(multiExpr[[1]]$data) = names(femData)[-c(1:8)];
multiExpr[[2]] = list(data = as.data.frame(t(maleData[-c(1:8)])));
names(multiExpr[[2]]$data) = maleData$substanceBXH;
rownames(multiExpr[[2]]$data) = names(maleData)[-c(1:8)];
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3)
gsg$allOK

#We now cluster the samples on their Euclidean distance, separately in each set.
sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

#pdf(file = "Plots/SampleClustering.pdf", width = 5, height = 5);
sizeGrWindow(5,5)
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7)
#Makes graphics disappear
#dev.off()
#There is one outlier F2_221 - now how to remove that. First choose the "base" cut height for the female data set
baseHeight = 16
# Adjust the cut height for the male data set for the number of samples
cutHeights = c(16, 16*exprSize$nSamples[2]/exprSize$nSamples[1]);
# Re-plot the dendrograms including the cut lines
#pdf(file = "Plots/SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
{
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7)
  abline(h=cutHeights[set], col = "red")
}
#dev.off()

#Now comes the actual outlier removal:
  for (set in 1:nSets)
  {
    # Find clusters cut by the line
    labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
    # Keep the largest one (labeled by the number 1)
    keep = (labels==1)
    multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
  }
collectGarbage();
# Check the size of the leftover data
exprSize = checkSets(multiExpr)
exprSize

#Add clinical traits for analysis
traitData = read.csv("ClinicalTraits.csv");
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)]
allTraits = allTraits[, c(2, 11:36) ]
# See how big the traits are and what are the trait and sample names
dim(allTraits)
names(allTraits)
allTraits$Mice
# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data)
  traitRows = match(setSamples, allTraits$Mice)
  Traits[[set]] = list(data = allTraits[traitRows, -1])
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1]
}
collectGarbage();
# Define data set dimensions
nGenes = exprSize$nGenes
nSamples = exprSize$nSamples

#Save that
save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize,
     file = "Consensus-dataInput.RData")

# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData")
# The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets

# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
#pdf(file = "Plots/scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}

#Plot with Scale free topology model fit - median, mean and max connectivity is shown - from this a power of 6 is chosen for further analysis threshold
net = blockwiseConsensusModules(
  multiExpr, power = 6, minModuleSize = 30, deepSplit = 2,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)
names(net)

#set colour selection
consMEs = net$multiMEs
moduleLabels = net$colors
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]]

#create dendrogram plot
sizeGrWindow(8,6);
#pdf(file = "Plots/ConsensusDendrogram-auto.pdf", wi = 8, he = 6)
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
#save it
save(consMEs, moduleLabels, moduleColors, consTree, file = "Consensus-NetworkConstruction-auto.RData")

#Relating consensus modules to external microarray sample information and exporting network analysis results

# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Also load results of network analysis
lnames = load(file = "Consensus-NetworkConstruction-auto.RData")
lnames
exprSize = checkSets(multiExpr)
nSets = exprSize$nSets

# Set up variables to contain the module-trait correlations
moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in 1:nSets)
{
  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
}

# Convert numerical lables to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");
# Open a suitably sized window (the user should change the window size if necessary)
sizeGrWindow(10,7)
#pdf(file = "Plots/ModuleTraitRelationships-female.pdf", wi = 10, he = 7);

# Plot the module-trait relationship table for set number 1
set = 1
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
                   signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))

# Plot the module-trait relationship table for set number 2
set = 2
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
                   signif(moduleTraitPvalue[[set]], 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
#pdf(file = "Plots/ModuleTraitRelationships-male.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2))
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))

#Combine the plots and highlight the differences
# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
# Find consensus negative correlations
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative])
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative])
# Find consensus positive correlations
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive])
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive])

#create heatmap
textMatrix = paste(signif(consensusCor, 2), "\n(",
                   signif(consensusPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
#pdf(file = "Plots/ModuleTraitRelationships-consensus.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = consensusCor,
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus module--trait relationships across\n",
                            paste(setLabels, collapse = " and ")))

#Exporting results
file = gzfile(description = "GeneAnnotation.csv");
annot = read.csv(file = file);
# Match probes in the data set to the probe IDs in the annotation file
probes = names(multiExpr[[1]]$data)
probes2annot = match(probes, annot$substanceBXH)

consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
GS = list()
kME = list()
for (set in 1:nSets)
{
  GS[[set]] = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data)
  kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data)
}

GS.metaZ = (GS[[1]]$Z + GS[[2]]$Z)/sqrt(2)
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z)/sqrt(2)
GS.metaP = 2*pnorm(abs(GS.metaZ), lower.tail = FALSE)
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE)

GSmat = rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[1]]$p, GS[[2]]$p, GS.metaZ, GS.metaP)
nTraits = checkSets(Traits)$nGenes
traitNames = colnames(Traits[[1]]$data)
dim(GSmat) = c(nGenes, 6*nTraits)
rownames(GSmat) = probes
colnames(GSmat) = spaste(
  c("GS.set1.", "GS.set2.", "p.GS.set1.", "p.GS.set2.", "Z.GS.meta.", "p.GS.meta"),
  rep(traitNames, rep(6, nTraits)))
# Same code for kME:
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP)
MEnames = colnames(consMEs.unord[[1]]$data)
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = probes
colnames(kMEmat) = spaste(
  c("kME.set1.", "kME.set2.", "p.kME.set1.", "p.kME.set2.", "Z.kME.meta.", "p.kME.meta"),
  rep(MEnames, rep(6, nMEs)))

info = data.frame(Probe = probes, GeneSymbol = annot$gene_symbol[probes2annot],
                  EntrezID = annot$LocusLinkID[probes2annot],
                  ModuleLabel = moduleLabels,
                  ModuleColor = labels2colors(moduleLabels),
                  GSmat,
                  kMEmat)
write.csv(info, file = "consensusAnalysis-CombinedNetworkResults.csv",
          row.names = FALSE, quote = FALSE)
