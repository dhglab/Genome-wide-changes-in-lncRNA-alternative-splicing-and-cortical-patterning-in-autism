## Code related to "Genome-wide changes in lncRNA, alternative splicing, and cortical patterning in autism"
## This script takes count data and filters for expressed genes using the criteria described in the Supplementary Information.
## It then uses conditional quantile normalization to perform gene length and GC normalization across samples and adjusts the data.
## Note that, for ease of sharing the code, the data are normalized without additional samples in the set (resulting in differences) at the GC normalization step and bootstrapping/permutation steps are not performed. Results are compared with those reported in the manuscript to demonstrate there are only minor differences and the final results are robust to these modifications on the data.

## Set current R session to be in the appropriate working directory

##### Set libraries and defaults and load provided files
rm(list=ls())
options(stringsAsFactors=FALSE)
load(file="../data/provided/ASD_RNAseq_ExpressionData.Rdata")

## Check the contents
summary(datMeta)

##### Split data into CTX and CBL and identify expressed genes
## Mark cortical vs cerebeallar samples
keepCTX <- datMeta[,"RegionID"] != "vermis"
keepCBL <- datMeta[,"RegionID"] == "vermis"

## Filtering Cufflinks gene expression - Cufflinks provides a lower bound for expression. If the lower bound of expression is > 0 in 80% of samples, we'll keep the genes. We filter in cortex and cerebellum separately, as they might have different genes expressed at high levels.
datExpr <- datExpr.Cufflinks[,keepCTX]
passvec <- apply(datLB.gene[,keepCTX]>0,1,sum)>(0.8*sum(keepCTX)) ## OK in 80% of samples
datExpr.Cufflinks.CTX <- datExpr[passvec,]

datExpr <- datExpr.Cufflinks[,keepCBL]
passvec <- apply(datLB.gene[,keepCBL]>0,1,sum)>(0.8*sum(keepCBL)) ## OK in 80% of samples
datExpr.Cufflinks.CBL <- datExpr[passvec,]

length(intersect(rownames(datExpr.Cufflinks.CTX),rownames(datExpr.Cufflinks.CBL))) ## 19161 overlap

length(setdiff(rownames(datExpr.Cufflinks.CTX),rownames(datExpr.Cufflinks.CBL))) ## 19402 overlap 1654 are "cortex specific"
ctxmean <- apply(datExpr.Cufflinks.CTX,1,mean)

## For the HTSeq counts data, we use a simple threshold of observing at least 10 fragments in 80% of the samples. This is done within cortex and cerebellum for the whole gene (including introns) level data and the exon union model data. We also compare the two data types at the beginning
geneIntersect <- intersect(rownames(datExpr.HTSC.unionexon),substr(rownames(datExpr.HTSC.wholegene),1,15))
geneInt1 <- datExpr.HTSC.unionexon[match(geneIntersect,rownames(datExpr.HTSC.unionexon)),]
geneInt2 <- datExpr.HTSC.wholegene[match(geneIntersect,substr(rownames(datExpr.HTSC.wholegene),1,15)),]

datExpr <- datExpr.HTSC.wholegene[,keepCTX]
passvec <- apply(datExpr,1,quantile,0.8) > 10
datExpr.HTSC.wholegene.CTX <- datExpr[passvec,]

datExpr <- datExpr.HTSC.wholegene[,keepCBL]
passvec <- apply(datExpr,1,quantile,0.8) > 10
datExpr.HTSC.wholegene.CBL <- datExpr[passvec,]

length(intersect(rownames(datExpr.HTSC.wholegene.CTX),rownames(datExpr.HTSC.wholegene.CBL))) ## 21206 intersect

datExpr <- datExpr.HTSC.unionexon[,keepCTX]
passvec <- apply(datExpr,1,quantile,0.8) > 10
datExpr.HTSC.unionexon.CTX <- datExpr[passvec,]

datExpr <- datExpr.HTSC.unionexon[,keepCBL]
passvec <- apply(datExpr,1,quantile,0.8) > 10
datExpr.HTSC.unionexon.CBL <- datExpr[passvec,]

length(intersect(rownames(datExpr.HTSC.unionexon.CTX),rownames(datExpr.HTSC.unionexon.CBL))) ##  22995 intersect

##### Apply conitional quantile normalization to obtained normalized FPKM values
library(cqn)
load("../data/provided/GC18unionAnno.Rdata") ## Modified version of Gencode v18 .gtf file to include union exon gene lengths
gc18unionAnno <- gc18unionAnno[gc18unionAnno[,1]>200,] ## Keep only the genes with length > 200 - the others are hard to assess and would be difficult to map to anyway                                                                                                                                                  

keepvec <- intersect(substr(rownames(gc18unionAnno),1,15),rownames(datExpr.HTSC.unionexon.CTX))
geneAnno <- gc18unionAnno[match(keepvec,substr(rownames(gc18unionAnno),1,15)),]
rownames(geneAnno) <- substr(rownames(geneAnno),1,15)

## For cortex, HTSeq union exon
datExpr.HTSC.CTX <- datExpr.HTSC.unionexon.CTX[match(rownames(geneAnno),rownames(datExpr.HTSC.unionexon.CTX)),]

cqn.dat <- cqn(datExpr.HTSC.CTX,lengths = as.numeric(geneAnno[,1]), x = as.numeric(geneAnno[,2]),lengthMethod=c("smooth"),sqn=FALSE) ## Run cqn with specified depths and no quantile normalization                                                                                                                      
cqn.dat <- cqn.dat$y + cqn.dat$offset ## Get the log2(Normalized FPKM) values                                                                                

columnNames <- colnames(datExpr.HTSC.CTX)
rowNames <- rownames(datExpr.HTSC.CTX)
datExpr.HTSC.CTX <- cqn.dat
rownames(datExpr.HTSC.CTX) <- rowNames
colnames(datExpr.HTSC.CTX) <- columnNames

datMeta.CTX <- datMeta[match(colnames(datExpr.HTSC.CTX),rownames(datMeta)),] ## Already filtered for include/exclude                                         

## For cortex, Cufflinks union exon
keepvec <- intersect(substr(rownames(gc18unionAnno),1,15),rownames(datExpr.Cufflinks.CTX))
geneAnno <- gc18unionAnno[match(keepvec,substr(rownames(gc18unionAnno),1,15)),]
rownames(geneAnno) <- substr(rownames(geneAnno),1,15)

datExpr.Cufflinks.CTX <- datExpr.Cufflinks.CTX[match(rownames(geneAnno),rownames(datExpr.Cufflinks.CTX)),]

cqn.dat <- cqn(datExpr.Cufflinks.CTX,lengths = 1000, x = as.numeric(geneAnno[,2]),lengthMethod=c("fixed"),sqn=FALSE) ## Run cqn with specified depths and no quantile normalization                                                                                                                                      
cqn.dat <- cqn.dat$y + cqn.dat$offset ## Get the log2(Normalized FPKM) values                                                                                

columnNames <- colnames(datExpr.Cufflinks.CTX)
rowNames <- rownames(datExpr.Cufflinks.CTX)
datExpr.Cufflinks.CTX <- cqn.dat
rownames(datExpr.Cufflinks.CTX) <- rowNames
colnames(datExpr.Cufflinks.CTX) <- columnNames

datMeta.CTX <- datMeta[match(colnames(datExpr.Cufflinks.CTX),rownames(datMeta)),] ## Already filtered for include/exclude                                    

## Intersect to get normalized FPKM matrix with expression values that pass filtering criteria
htscint.CTX <- intersect(substr(rownames(datExpr.HTSC.wholegene.CTX),1,15),
                         intersect(rownames(datExpr.HTSC.CTX),
                                   rownames(datExpr.Cufflinks.CTX)))
length(htscint.CTX) ## 16310 genes

## For cerebellum, HTSeq union exon
keepvec <- intersect(substr(rownames(gc18unionAnno),1,15),rownames(datExpr.HTSC.unionexon.CBL))
geneAnno <- gc18unionAnno[match(keepvec,substr(rownames(gc18unionAnno),1,15)),]
rownames(geneAnno) <- substr(rownames(geneAnno),1,15)

datExpr.HTSC.CBL <- datExpr.HTSC.unionexon.CBL[match(rownames(geneAnno),rownames(datExpr.HTSC.unionexon.CBL)),]

cqn.dat <- cqn(datExpr.HTSC.CBL,lengths = as.numeric(geneAnno[,1]), x = as.numeric(geneAnno[,2]),lengthMethod=c("smooth"),sqn=FALSE) ## Run cqn with specified depths and no quantile normalization                                                                                                                      
cqn.dat <- cqn.dat$y + cqn.dat$offset ## Get the log2(Normalized FPKM) values                                                                                

columnNames <- colnames(datExpr.HTSC.CBL)
rowNames <- rownames(datExpr.HTSC.CBL)
datExpr.HTSC.CBL <- cqn.dat
rownames(datExpr.HTSC.CBL) <- rowNames
colnames(datExpr.HTSC.CBL) <- columnNames

datMeta.CBL <- datMeta[match(colnames(datExpr.HTSC.CBL),rownames(datMeta)),] ## Already filtered for include/exclude                                         

## For cerebellum, Cufflinks union exon
keepvec <- intersect(substr(rownames(gc18unionAnno),1,15),rownames(datExpr.Cufflinks.CBL))
geneAnno <- gc18unionAnno[match(keepvec,substr(rownames(gc18unionAnno),1,15)),]
rownames(geneAnno) <- substr(rownames(geneAnno),1,15)

datExpr.Cufflinks.CBL <- datExpr.Cufflinks.CBL[match(rownames(geneAnno),rownames(datExpr.Cufflinks.CBL)),]

cqn.dat <- cqn(datExpr.Cufflinks.CBL,lengths = 1000, x = as.numeric(geneAnno[,2]),lengthMethod=c("fixed"),sqn=FALSE) ## Run cqn with specified depths and no quantile normalization                                                                                                                                      
cqn.dat <- cqn.dat$y + cqn.dat$offset ## Get the log2(Normalized FPKM) values                                                                                

columnNames <- colnames(datExpr.Cufflinks.CBL)
rowNames <- rownames(datExpr.Cufflinks.CBL)
datExpr.Cufflinks.CBL <- cqn.dat
rownames(datExpr.Cufflinks.CBL) <- rowNames
colnames(datExpr.Cufflinks.CBL) <- columnNames

datMeta.CBL <- datMeta[match(colnames(datExpr.Cufflinks.CBL),rownames(datMeta)),] ## Already filtered for include/exclude                                    

## Intersect to get normalized FPKM matrix with expression values that pass filtering criteria
htscint.CBL <- intersect(substr(rownames(datExpr.HTSC.wholegene.CBL),1,15),
                      intersect(rownames(datExpr.HTSC.CBL),
                                 rownames(datExpr.Cufflinks.CBL)))
length(htscint.CBL) ## 15824

## Save cortex and cerebellum data
datExpr.CTX <- datExpr.HTSC.CTX[match(htscint.CTX,rownames(datExpr.HTSC.CTX)),]
datExpr.CBL <- datExpr.HTSC.CBL[match(htscint.CBL,rownames(datExpr.HTSC.CBL)),]

save(datExpr.CTX,datMeta.CTX,datExpr.CBL,datMeta.CBL,file="../data/modified/NormalizedExpressionData.Rdata")

##### Identify covariates to compute adjusted FPKM values
## Note that this step was performed with a bootstrapped estimate of covariates computed on the subset of matched samples for the published analysis.

## Cortical data
## Extract the covariate data                                                                                                                                    
condition <- 2-as.numeric(as.factor(datMeta.CTX[,"ASD.CTL"]))
age <- as.numeric(datMeta.CTX[,"Age"])
sex <- as.numeric(as.factor(datMeta.CTX[,"Sex"]))-1
region <- as.numeric(as.factor(datMeta.CTX[,"RegionID"]))-1

batch <- as.numeric(as.factor(datMeta.CTX[,"SeqBatch"]))-1
RIN <- as.numeric(datMeta.CTX[,"RIN"])
bank <- as.numeric(as.factor(datMeta.CTX[,"BrainBank"]))-1

datSeq <- data.matrix(datMeta.CTX[,c(25:43)])
datSeq <- datSeq[,c(6:19,c(5))]
datSeqNorm <- t(scale(datSeq,scale=F))
PC.datSeq <- prcomp(datSeqNorm);
varexp <- (PC.datSeq$sdev)^2 / sum(PC.datSeq$sdev^2)
print(varexp[1:5])
topPC.datSeq <- PC.datSeq$rotation[,1:2];
colnames(topPC.datSeq) <- c("SeqSV1","SeqSV2") ## Recompute since those in datMeta were with additional samples included
seqStatPC1 <- as.numeric(topPC.datSeq[,1])
seqStatPC2 <- as.numeric(topPC.datSeq[,2])

regvars <- as.data.frame(cbind(condition,age,sex,region,batch,RIN,bank,seqStatPC1,seqStatPC2))

## Run the regression and make the adjusted FPKM matrix                                                                                                                                        
datExpr.reg <- matrix(NA,nrow=nrow(datExpr.CTX),ncol=ncol(datExpr.CTX))
rownames(datExpr.reg) <- rownames(datExpr.CTX)
colnames(datExpr.reg) <- colnames(datExpr.CTX)
coefmat <- matrix(NA,nrow=nrow(datExpr.CTX),ncol=ncol(regvars)+1)

for (i in 1:nrow(datExpr.CTX)) {
  lmmod1 <- lm(as.numeric(datExpr.CTX[i,])~condition+age+sex+region+batch+RIN+bank+seqStatPC1+seqStatPC2,data=regvars)
  coef <- coef(lmmod1)
  coefmat[i,] <- coef
  datExpr.reg[i,] <- datExpr.CTX[i,] - coef["RIN"]*regvars[,"RIN"] - coef["bank"]*regvars[,"bank"] - coef["seqStatPC1"]*regvars[,"seqStatPC1"] - coef["seqStatPC2"]*regvars[,"seqStatPC2"]
}
datExpr.CTX.reg <- datExpr.reg ## This datExpr.CTX.reg is now a technical variable corrected matrix.
rownames(datExpr.CTX.reg) <- rownames(datExpr.CTX)
colnames(datExpr.CTX.reg) <- rownames(datMeta.CTX)

## Cerebellar data
## Extract the covariate data                                                                                                                                    
condition <- 2-as.numeric(as.factor(datMeta.CBL[,"ASD.CTL"]))
age <- as.numeric(datMeta.CBL[,"Age"])
sex <- as.numeric(as.factor(datMeta.CBL[,"Sex"]))-1

batch <- as.numeric(as.factor(datMeta.CBL[,"SeqBatch"]))-1
RIN <- as.numeric(datMeta.CBL[,"RIN"])
bank <- as.numeric(as.factor(datMeta.CBL[,"BrainBank"]))-1

datSeq <- data.matrix(datMeta.CBL[,c(25:43)])
datSeq <- datSeq[,c(6:19,c(5))]
datSeqNorm <- t(scale(datSeq,scale=F))
PC.datSeq <- prcomp(datSeqNorm);
varexp <- (PC.datSeq$sdev)^2 / sum(PC.datSeq$sdev^2)
print(varexp[1:5])
topPC.datSeq <- PC.datSeq$rotation[,1:2];
colnames(topPC.datSeq) <- c("SeqSV1","SeqSV2") ## Recompute since those in datMeta were with additional samples included
seqStatPC1 <- as.numeric(topPC.datSeq[,1])
seqStatPC2 <- as.numeric(topPC.datSeq[,2])

regvars <- as.data.frame(cbind(condition,age,sex,batch,RIN,bank,seqStatPC1,seqStatPC2))

## Run the regression and make the adjusted FPKM matrix                                                                                                                                        
datExpr.reg <- matrix(NA,nrow=nrow(datExpr.CBL),ncol=ncol(datExpr.CBL))
rownames(datExpr.reg) <- rownames(datExpr.CBL)
colnames(datExpr.reg) <- colnames(datExpr.CBL)
coefmat <- matrix(NA,nrow=nrow(datExpr.CBL),ncol=ncol(regvars)+1)

for (i in 1:nrow(datExpr.CBL)) {
  lmmod1 <- lm(as.numeric(datExpr.CBL[i,])~condition+age+sex+batch+RIN+bank+seqStatPC1+seqStatPC2,data=regvars)
  coef <- coef(lmmod1)
  datExpr.reg[i,] <- datExpr.CBL[i,] - coef["RIN"]*regvars[,"RIN"] - coef["bank"]*regvars[,"bank"] - coef["seqStatPC1"]*regvars[,"seqStatPC1"] - coef["seqStatPC2"]*regvars[,"seqStatPC2"]
}
datExpr.CBL.reg <- datExpr.reg ## This datExpr.CBL.reg is now a technical variable corrected matrix.
rownames(datExpr.CBL.reg) <- rownames(datExpr.CBL)
colnames(datExpr.CBL.reg) <- rownames(datMeta.CBL)

## Save adjusted FPKM data
save(datExpr.CTX.reg,datMeta.CTX,datExpr.CBL.reg,datMeta.CBL,file="../data/modified/AdjustedExpressionData.Rdata")
