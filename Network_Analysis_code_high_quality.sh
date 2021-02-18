cd /Volumes/Hyden/blh226/WGCNA 
R
###########################installation###################################
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R") 
biocLite(c("GO.db", "preprocessCore", "impute")) 
orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg", "At", "Pt"); 
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6)); 
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep=""); 
biocLite(c("GO.db", "KEGG.db", "topGO", packageNames, "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList"))
install.packages("BiocManager") 
library(BiocManager)
library(WGCNA)
options(stringsAsFactors = FALSE)
###########read in files####################
#FPKM<-read.table('FPKM_combined_high_quality.txt', header=FALSE)
#Samples<-FPKM[1,]
FPKM<-read.table('FPKM_combined_high_quality.txt')
dim(FPKM)
names(FPKM)
datExpr0<-as.data.frame(t(FPKM))

gsg<-goodSamplesGenes(datExpr0, verbose=3)
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
##################filter######################
traitData<-read.table('/Volumes/Hyden/blh226/WGCNA/sex_high_quality.txt', row.names=1)
allTraits<-as.data.frame(t(traitData))
names(allTraits)<-c('plant', 'sex')
Samples<-rownames(datExpr0)
traitRows = match(Samples, allTraits$plant)
datTraits = allTraits[traitRows,]
datExpr<-datExpr0

########################cluster########################
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(as.numeric(as.factor(datTraits[,2])), signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf('plotDendroAndColors_HQ.pdf', width=24, height=12)
plotDendroAndColors(sampleTree2, traitColors,
                  groupLabels = names(datTraits),
                  main = "Sample dendrogram and trait heatmap")

# Plot a line to show the cut
abline(h = 30000, col = "red")
dev.off()
# Determine cluster under the line
clust = cutreeStatic(sampleTree2, cutHeight = 30000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

powers = c(c(1:10), seq(from = 12, to=20, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
pdf('Mean Connectivity and Soft Threshold HQ.pdf', height=12, width=24)
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
dev.off()
softPower = 15;
adjacency = adjacency(datExpr, power = softPower);
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
save(dissTOM, file="networkanalysisHQ_dissTOM.RData")

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
pdf('Gene clustering TOM HQ.pdf', width=24, height=12)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE, hang = 0.04)
dev.off()
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
              deepSplit = 2, pamRespectsDendro = FALSE,
              minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf('Dendro and Colors HQ.pdf', height=12, width=24)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                  dendroLabels = FALSE, hang = 0.03,
                  addGuide = TRUE, guideHang = 0.05,
                  main = "Gene dendrogram and module colors")
dev.off()
##############merge co-expressed modules##################
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('Clustering of Module Eigengenes HQ.pdf', height=12, width=24)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

pdf('Dynamic Tree Cut HQ.pdf', height=12, width=24)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(200));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

colorOrder2 = c("grey", standardColors(200))
moduleLabels2 = match(dynamicColors, colorOrder)-1
# Save module colors and labels for use in subsequent parts

save(MEs, moduleLabels, moduleLabels2, moduleColors, dynamicColors, geneTree, file = "networkanalysisHQ_merged.RData")

#################Relating Modules to external information###################

#lnames = load(file = "networkanalysisfinal_merged.RData")
#lnames

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs1 = moduleEigengenes(datExpr, dynamicColors)$eigengenes

#filter datTraits to correspond with removed sample 
datTraits<-datTraits[-80,]
datTraits<-datTraits[-15,]
datTraits2<-datTraits[,2]
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, as.numeric(as.factor(datTraits2))-2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


pdf('Module-trait Relationships Heatmap merged HQ extra narrow.pdf', width=4, height=12)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
             xLabels = 'sex',
             yLabels = names(MEs),
             ySymbols = names(MEs),
             colorLabels = FALSE,
             colors = blueWhiteRed(200),
             textMatrix = textMatrix,
             setStdMargins = FALSE,
             cex.text = 0,
             zlim = c(-1,1),
             main = paste("Module-trait relationships"))

dev.off()


MEs1 = orderMEs(MEs1)
moduleTraitCor = cor(MEs1, as.numeric(as.factor(datTraits2))-2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


pdf('Module-trait Relationships Heatmap unmerged HQ extra narrow.pdf', width=4, height=15)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
             xLabels = 'sex',
             yLabels = names(MEs1),
             ySymbols = names(MEs1),
             colorLabels = FALSE,
             colors = blueWhiteRed(200),
             textMatrix = textMatrix,
             setStdMargins = FALSE,
             cex.text = 0,
             zlim = c(-1,1),
             main = paste("Module-trait relationships"))
dev.off()




# Define variable weight containing the weight column of datTrait
sex<-as.numeric(as.factor(datTraits2))-2
names(sex)<-datTraits[,1]
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, sex, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
#names(geneTraitSignificance) = paste("GS.", names(sex), sep="");
#names(GSPvalue) = paste("p.GS.", names(sex), sep="");

modNames1 = substring(names(MEs1), 3)
geneModuleMembership1 = as.data.frame(cor(datExpr, MEs1, use = "p"));
MMPvalue1 = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership1), nSamples));
names(geneModuleMembership1) = paste("MM", modNames1, sep="");
names(MMPvalue1) = paste("p.MM", modNames1, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, sex, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
#names(geneTraitSignificance) = paste("GS.", names(sex), sep="");
#names(GSPvalue) = paste("p.GS.", names(sex), sep="");


for(i in levels(as.factor(moduleColors))){
  module = i
  column = match(module, modNames1);
  moduleGenes = dynamicColors==module;
  pdf(paste(module, "membership vs gene significance.pdf"))
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership1[moduleGenes, column]),                
                  abs(geneTraitSignificance[moduleGenes, 1]),
                  xlab = paste("Module Membership in", module, "module"),
                  ylab = "Gene significance for sex",
                  main = paste("Module membership vs. gene significance\n"),
                  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
}
for(i in levels(as.factor(dynamicColors))){
  module = i
  column = match(module, modNames1);
  moduleGenes = dynamicColors==module;
  pdf(paste(module, "membership vs gene significance_unmerged.pdf"))
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership1[moduleGenes, column]),                
                  abs(geneTraitSignificance[moduleGenes, 1]),
                  xlab = paste("Module Membership in", module, "module"),
                  ylab = "Gene significance for sex",
                  main = paste("Module membership vs. gene significance\n"),
                  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
}



# Read in the probe annotation
annot = read.delim(file = "Spurpurea_519_v5.1.defline.txt", sep='\t', header=FALSE);
newnames<-paste0(annot[,1], '.v5.1')
annotnew<-cbind(newnames, annot)

FPKM2<-read.table('FPKM_and_miRNA_combined_counts.txt')

probes = row.names(FPKM2)[as.numeric(names(datExpr))]
probes2annot<-match(probes, annotnew[,1])
# Create the starting data frame
geneInfo0 = data.frame(GeneID = probes,
										annotation = annot[probes2annot,3],
                    moduleColormerged = moduleColors,
                    moduleColorunmberged = dynamicColors,
                    geneTraitSignificance,
                    GSPvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, sex, use = "p")));
names(geneInfo0)<-c('GeneID', 'defline', 'moduleColormerged', 'moduleColorunmerged', 'GeneTraitSignificance', 'GSPvalue')

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
geneOrder = order(geneInfo0$moduleColormerged, -abs(geneInfo0$GeneTraitSignificance));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "WGCNA_HQ_results.csv")

#Make figures
RNA_Seq<-read.csv('All_RNASeq_analysis_combined.csv', header=TRUE, stringsAsFactors=FALSE)
library(pheatmap)
DE<-read.csv('DESEQ2 M vs F high quality with counts.csv', header=TRUE, stringsAsFactors=FALSE, row.names=1)
RNA_Seq_counts<-match(RNA_Seq$transcriptID, DE$Row.names)
All<-cbind(RNA_Seq, DE[RNA_Seq_counts, ])
write.csv(All, 'All_data_and_counts.csv')
Diff<-All[All[,10]<0.05,]
Diff<-Diff[!is.na(Diff$padj),]
Extreme1<-Diff[Diff$log2FoldChange.M.F.HQ>1.0,]
Extreme2<-Diff[Diff$log2FoldChange.M.F.HQ<-1.0,]
Extreme<-rbind(Extreme1,Extreme2)
purple<-Extreme[Extreme$moduleColormergedHQ=='purple',]
cyan<-Extreme[Extreme$moduleColormergedHQ=='cyan',]
quantile.range <- quantile(as.matrix(Extreme[,27:183]), probs = seq(0, 1, 0.01))
palette.breaks <- seq(quantile.range["0%"], quantile.range["90%"], 0.1)
color.palette  <- colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(length(palette.breaks) - 1)
pheatmap(cyan[,27:183], border_color=NA, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, col=color.palette, breaks=palette.breaks)


