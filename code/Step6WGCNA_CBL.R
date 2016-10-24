## Code related to "Genome-wide changes in lncRNA, alternative splicing, and cortical patterning in autism"
## This script will perform network analysis in CB samples
## Note that, for ease of sharing the code, the data are adjusted without bootstrapping.
## Additionally, WGCNA is done without bootstrapping - for the study we constructed 100 topological overlap matrices by resampling the full sample set and took the median edge strength, providing a robust output of modules.
## At the end of the script, we compare with what is published and match up modules.

rm(list=ls())
options(stringsAsFactors=FALSE)
library(WGCNA)

## Load FPKM expression matrix
load("../data/modified/AdjustedExpressionData.Rdata") ##datExpr.CBL.reg and datMeta.CBL

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(t(datExpr.CBL.reg), powerVector = powers, verbose = 5,networkType="signed",corFnc="bicor")
par(mfrow = c(1,2)); # Plot the results:
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

## Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
softPower=9 ## Choose a power of 9 as it's the first on where the R^2 is > 0.8

## Run an automated network analysis
net3 <- blockwiseModules(t(datExpr.CBL.reg),power=softPower,deepSplit=4,minModuleSize=100,
                         mergeCutHeight=0.2, detectCutHeight=0.9999,
                         corType="bicor",networkType="signed",pamStage=FALSE,pamRespectsDendro=TRUE,
                         verbose=3,saveTOMs=FALSE,maxBlockSize=30000)
table(net3$colors)
length(table(net3$colors))

## Match module Colors with rWGCNA module definitions
moduleColors_CBL.rWGCNA <- read.csv('../data/provided/SuppTable2B.csv') # load rWGCNA CBL module Definitions
gnS <- intersect(rownames(datExpr.CBL.reg),moduleColors_CBL.rWGCNA$ENSEMBL.ID)
moduleColors_CBL.rWGCNA <- moduleColors_CBL.rWGCNA[match(gnS,moduleColors_CBL.rWGCNA$ENSEMBL.ID),]
moduleColors_CBL.rWGCNA_cols <- as.character(moduleColors_CBL.rWGCNA$WGCNA.module.color)
moduleColors <- net3$colors[match(gnS,rownames(datExpr.CBL.reg))]
datExpr.CBL.reg <- datExpr.CBL.reg[match(gnS,rownames(datExpr.CBL.reg)),]

matchLabel1 <- matchLabels(moduleColors, moduleColors_CBL.rWGCNA_cols) #matchLabels(source,reference)
moduleColors.old <- moduleColors
merge <- mergeCloseModules(t(datExpr.CBL.reg), matchLabel1, cutHeight = 0.2, verbose = 3)

# The merged module colors
mergedColors <- merge$colors;
moduleColors <- mergedColors
table(moduleColors)

# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs;
MEs <- mergedMEs;
MEs <- orderMEs(MEs)
KMEs <- signedKME(t(datExpr.CBL.reg), MEs,outputColumnName = "kME",corFnc = "bicor")

## Annotate the output
library(biomaRt)
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","band","gene_biotype")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl") # add host="sep2013.archive.ensembl.org" for Gencode v18
geneAnno1 <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values=rownames(datExpr.CBL.reg),mart=mart)
geneAnno2 <- geneAnno1[match(rownames(datExpr.CBL.reg),geneAnno1[,1]),]
geneInfo <- as.data.frame(cbind(geneAnno2,moduleColors, KMEs))

## Comparison with rWGCNA Modules
geneInfo.rWGCNA <- moduleColors_CBL.rWGCNA
background.rWGCNA <- as.character(unique(geneInfo.rWGCNA$ENSEMBL.ID))
geneInfo.ASD <- geneInfo
background.ASD <- as.character(unique(geneInfo.ASD$ensembl_gene_id))
uniquemodcolors.rWGCNA <- as.character(unique(geneInfo.rWGCNA$WGCNA.Module.Label))
uniquemodcolors.rWGCNA <- uniquemodcolors.rWGCNA[order(uniquemodcolors.rWGCNA)]
uniquemodcolors.ASD <- as.character(unique(geneInfo.ASD$moduleColors))
uniquemodcolors.ASD <- uniquemodcolors.ASD[order(uniquemodcolors.ASD)]

ORmat <- Pmat <- matrix(NA,nrow=length(uniquemodcolors.ASD),ncol=1)
source("./ORA.R")
for (i in 1:length(uniquemodcolors.rWGCNA)){
  thismod <- uniquemodcolors.rWGCNA[i]
  thisGene <- geneInfo.rWGCNA[geneInfo.rWGCNA$WGCNA.Module.Label==thismod,"ENSEMBL.ID"]
  testpath <- as.character(unique(thisGene)) ## Module Genes
  oraMat1 <- matrix(NA,ncol=9,nrow=length(uniquemodcolors.ASD))
  
  for(j in 1:length(uniquemodcolors.ASD)){
    thismod1 <- uniquemodcolors.ASD[j]
    thisGene1 <- geneInfo.ASD[geneInfo.ASD$moduleColors==thismod1,"ensembl_gene_id"]
    refpath <- as.character(unique(thisGene1)) ## Module Genes
    
    testbackground <- background.rWGCNA
    refbackground <- background.ASD
    
    oraout <- ORA(testpath,refpath,testbackground,refbackground)
    
    oraMat1[j,] <- oraout
  }
  ORmat <- cbind(ORmat,as.numeric(oraMat1[,1]))
  Pmat <- cbind(Pmat,as.numeric(oraMat1[,2]))
}

ORmat.ASD <- ORmat[,-1]
Pmat.ASD <- Pmat[,-1]
FDRmat.ASD <- matrix(p.adjust(Pmat.ASD,method="BH"),nrow=nrow(Pmat.ASD),ncol=ncol(Pmat.ASD))
colnames(ORmat.ASD) <- colnames(Pmat.ASD) <- colnames(FDRmat.ASD) <- paste('rWGCNA.CBL.M',uniquemodcolors.rWGCNA,sep='')
rownames(ORmat.ASD) <- rownames(Pmat.ASD) <- rownames(FDRmat.ASD) <- paste('CBL_ASD.',uniquemodcolors.ASD,sep='')
dispMat <- -log10(Pmat.ASD)*sign(log2(ORmat.ASD)) ## You can change this to be just log2(Bmat) if you want the color to reflect the odds ratios

## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g.
txtMat <- ORmat.ASD
txtMat[FDRmat.ASD>=0.05] <- ""
txtMat[FDRmat.ASD <0.05&FDRmat.ASD >0.01] <- ""
txtMat[FDRmat.ASD <0.01&FDRmat.ASD >0.005] <- ""
txtMat[FDRmat.ASD <0.005] <- ""
txtMat1 <- signif(ORmat.ASD,2)
txtMat1[txtMat1<5] <- ""
textMatrix1 <- paste( txtMat1, '\n', txtMat , sep = '');
textMatrix1 <- matrix(textMatrix1,ncol=ncol(Pmat.ASD),nrow=nrow(Pmat.ASD))

## Plot how the modules match up between what is computed here, and what is provided in the supplement
par(mfrow=c(1,1))
labeledHeatmap(Matrix=ORmat.ASD,
               yLabels= rownames(ORmat.ASD),
               yColorLabels=F,
               xLabels= colnames(ORmat.ASD),
               colors=blueWhiteRed(100),
               textMatrix = textMatrix1,
               cex.lab.x=1.0,
               zlim=c(-100,100),
               main="Signed -log10(p-value) Heatmap")
## The high overlaps (OR > 5) represent modules that are the same or, in some cases, split or merge between the analyses. We recommend using the rWGCNA results for downstream analyses (e.g. comparing with future studies) as the network analysis output computed here is a less robust version.