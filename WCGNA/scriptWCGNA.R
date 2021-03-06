##########################
#   WCGNA-processing      #
##########################
#
#by Francisco Jose de Novais
#This script to perform modules of co-expression analysis.
#
#INSTALATION
#source("http://bioconductor.org/biocLite.R") 
#biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
#library(BiocInstaller)
#biocValid()
#install.packages("WGCNA")
#load
#
#
library(WGCNA)
options(stringsAsFactors = FALSE) 
#INPUTTRANSCRIPTOMICDATA <- read.csv("intensitiesData.csv"); APPLY THESE STEP IF NECESSARY
datExpr0 <- as.data.frame(t(intensitiesData[, -1]));
names(datExpr0) <- intensities$X1 
rownames(datExpr0) <- names(intensities)[-1];
#
#CHECK OUTLIERS OR GENES WITH A LOT OF MISSING VALUES ##CHECAR OUTLIERS OU GENES COM MUITOS MISSING VALUES
gsg <- goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK 
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
sampleTree <- hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)#graphic dimmension 
par(cex = 1);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
#
#EXCLUDE OUTLIER
# Plot a line to show the cut
abline(h = 80, col = "red");
# Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = 70, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples <- (clust==1)
datExpr <- datExpr0[keepSamples, ]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
#
#INPUTPHENOTYPEDATA = read.csv("ClinicalTraits.csv") #phenotype data
# remove columns that hold information we do not need.
allTraits <- data.frame(Traits);
allTraits <- allTraits[, c(2,4)]
dim(allTraits)
names(allTraits)
#
# Form a data frame analogous to expression data that will hold the phenotype traits.
femaleSamples <- rownames(datExpr);
traitRows <- match(femaleSamples, allTraits$ANIMAL_ID);
TRAIT <- allTraits[traitRows,-1]
datTraits <- data.frame(TRAIT)
rownames(datTraits) <- rownames(datExpr)
collectGarbage();
#
# Re-cluster samples
sampleTree2 <- hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(datTraits, signed = T);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
#
#
#################################
# NETWORK ANALYSIS STEP-BY-STEP #
#################################
# Step-by-step network construction and module detection
#Choosing the soft-thresholding power: analysis of network topology
# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 3)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Co-expression similarity and adjacency
softPower = USEHERETHENUMBER; #APPLY THE NUMBER FOR THE SOFTPOWER
adjacency <- adjacency(datExpr, power = softPower);
#Topological Overlap Matrix (TOM)
# Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency);
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Feature clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 20;
# Module identification using dynamic tree cut:
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#
varExpModules <- propVarExplained(datExpr,dynamicColors, MEs, corFnc = "cor", corOptions = "use = 'p'")##CHECKING THE VARIANCE EXPLAINED BY FIRST COMPONENT OF THE MODULE
#
MEDissThres = 0.25 #0.75 of correlation
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
##Put together similar modules
sizeGrWindow(7, 6)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(mergedColors),
                    c("ModuleColor"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
propVarExplained(datExpr,mergedColors, mergedMEs, corFnc = "cor", corOptions = "use = 'p'")#check the proportion of variance explained by first component in each module
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
#Relating modules to external traits and identifying important genes
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#
sizeGrWindow(7,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 2), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 9, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = "TRAIT",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.55,
               zlim = c(-1,1),
               main = paste("Module-trait relationship"))
#
#
#Gene relationship to trait and important modules: Gene Significance and Module merbership
# Define variable weight containing the weight column of datTrait
TRAIT <- as.data.frame(datTraits$TRAIT);
names(TRAIT) = "TRAIT"
# names (colors) of the modules
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) <- paste("MM", modNames, sep="");
names(MMPvalue) <- paste("p.MM", modNames, sep="");
geneTraitSignificance <- as.data.frame(cor(datExpr,TRAIT, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(TRAIT), sep="");
names(GSPvalue) = paste("p.GS.", names(TRAIT), sep="");
#
#Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules
modNames
module = "nameofmodule" #STUDIED MODULE
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership"),
                   ylab = "Gene Significance",
                   main = paste("STUDIED MODULE\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")
#
propVarExplained(datExprM, mergedColors, MEs, corFnc = "cor", corOptions = "use = 'p'")
#
# Create the starting data frame
geneInfo0 <- data.frame(featureID = names(datExpr),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for TRAIT
modOrder = order(-abs(cor(MEs, TRAIT, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.RFI));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "PATHWAY/FILE.csv")
#
#
#SELECTING MM > 0.6 intramodule and < 0.6 in others modules
#DECIDE IF YOU WILL APPLY THIS PART OR NOT!!!!
#MMturquoise would be replaced using the significance module
genes.total <- data.frame(geneModuleMembership$MMturquoise >0.5 & moduleColors=="turquoise" & (apply(!geneModuleMembership[,c(1:18,20:22)] <0.5,1,sum)))
#
colnames(genes.total) <- c('MMturquoise')
#
test.total<-cbind(genes.total, t(datExpr))
#Select only the genes which are in one of the modules:
genes.total <- test.total[test.total$MMturquoise== T,]
#
#select genes which are in one particular module:
genes.turquoise <-test.total[test.total$MMturquoise==T,]
#
#
#The first 3 columns are from the MM of the modules, the next N columns are my samples
genes.turquoise <- genes.turquoise[,2:14]
#
#
###Calculate module eigengene per module, to make a barplot of those module eigengenes per module###
#you select the column of that module color with e.g. [2,] Ã  brown is the second column
MEs.totalselect <- t(MEs)
MEs.total.blue <- MEs.totalselect[5,]
MEs.total.darktur <- MEs.totalselect[7,]
MEs.total.turq <- MEs.totalselect[18,]
#
#DECIDE IF YOU WILL USE THIS PART OR NOT
#Create a barplot of one module, with samples on the x-axis (one bar per sample), colored per group (in my case three groups of 12 samples).
color = c(rep ("black",8), rep ("grey60", 8))
barplot(MEs.total.grey, col=color,space=c(rep(0.1,7), 0.5, rep(0.1,6)), cex.axis=1.5)
title(main = list("Grey60 Module", font = 2), cex.main=2, ylab="intensity (module eigengene)", cex.lab=1.5)
#
color = c(rep ("darkgreen",8), rep ("green", 8))
barplot(MEs.total.darktur, col=color,space=c(rep(0.1,7), 0.5, rep(0.1,6)), cex.axis=1.5)
title(main = list("DarkTurq Module", font = 2), cex.main=2, ylab="intensity (module eigengene)", cex.lab=1.5)
#
color = c(rep ("brown",8), rep ("orange", 8))
barplot(MEs.total.turq, col=color,space=c(rep(0.1,7), 0.5, rep(0.1,6)), cex.axis=1.5)
title(main = list("Turquoise Module", font = 2), cex.main=2, ylab="intensity (module eigengene)", cex.lab=1.5)
#
write.table (genes.blue, "/Users//bluemodule.txt", sep="\t")
write.table (genes.grey60, "/Users/.../wgcna/genesGrey.txt", sep="\t")
write.table (genes.turquoise, "../genesTurquoise.txt", sep="\t")
#
#
#############################################
#Network visualization using WGCNA functions#
#############################################
#
#Visualizing the gene network
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
#dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = NUMBER);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
#plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
#diag(plotTOM) = NA;
# Call the plot function
#sizeGrWindow(9,9)
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
#
#
###############
#    TRAIT     #
###############
# Isolate YOUR TRAIT from the traits
#TRAIT = as.data.frame(datTraits$TRAIT);
#names(TRAIT) = "TRAIT"
# Add the RFI to existing module eigengenes
MET = orderMEs(cbind(MEs, RFI))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab= 0.8, xLabelsAngle= 90)
#
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#                   Exporting to Cytoscape                 #
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# Select modules
modules <- c("turquoise");#example
# Select module probes
probes <- names(datExpr)
inModule <- is.finite(match(moduleColors, modules));
modProbes <- probes[inModule];
#modGenes = annot$geneName[match(modProbes,)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];

dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);

#calcular intramodular connectivity
IntraModularConn <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = FALSE)
IntraModularConn$moduleColors <- moduleColors
write.table (IntraModularConn, "PATHWAY/FILE.txt", sep="\t")
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
#Exporting a gene network to external visualization software#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
#Developed by Pamela Almeida Alexandre
#
#names <- read.table ("names.txt", sep="\t")
#names (names) = c ("geneID", "geneNames")
#annot <- data.frame (geneInfo$featureID)
#names (annot) <- c ("geneID")
#annot <- merge (geneID, names, by="geneID", all=T)
#write.table (annot, "annot.txt", sep="\t")
#
# Read in the annotation file
#annot = read.table ("annot.txt", sep="\t", header=T);
#

#That's all Folks!!
q()




