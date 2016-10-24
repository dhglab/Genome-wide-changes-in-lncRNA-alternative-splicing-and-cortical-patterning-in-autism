## Code related to "Genome-wide changes in lncRNA, alternative splicing, and cortical patterning in autism"
## This script will use the samples described in Table S1 to perform a paired analysis to show attenuated cortical patterning in ASD vs CTL
## Note that, for ease of sharing the code, the data are adjusted without bootstrapping, so the final ACP set differs slightly from what was published. However, we show that we get the same results (including SOX5 enrichment in transcription factor binding site analysis) despite this minor variation in how the data are analyzed.

rm(list=ls())
#library(Cairo)
library(WGCNA)
library(biomaRt)
library(NMF)
options(stringsAsFactors=FALSE)

## Start with the adjusted FPKM data and the samples used in the age-matched, paired analysis
load("../data/modified/AdjustedExpressionData.Rdata")
datSupp <- read.csv("../data/provided/SuppTab1C_pairedsamples.csv")
TC.samp.tags <- paste(datSupp[,"Brain.ID"],"ba41-42-22",datSupp[,"TC.RIN"],sep="_")
FC.samp.tags <- paste(datSupp[,"Brain.ID"],"ba9",datSupp[,"FC.RIN"],sep="_")
all.samp.tags <- paste(datMeta.CTX[,"BrainID"],datMeta.CTX[,"RegionID"],datMeta.CTX[,"RIN"],sep="_")

## Extract the paired samples
keepTC <- match(TC.samp.tags,all.samp.tags)
keepFC <- match(FC.samp.tags,all.samp.tags)

## Separate the data out into a mfc and a stc dataset
seq.dat.mfc <- datExpr.CTX.reg[,keepFC]
seq.dat.stc <- datExpr.CTX.reg[,keepTC]
met.dat.mfc <- datMeta.CTX[keepFC,]
met.dat.stc <- datMeta.CTX[keepTC,]

## Run a paired wilcoxon rank-sum test (also can switch this to a t-test, they give similar results)
asd.seq.dat.t.vec <- asd.seq.dat.p.vec <- asd.seq.dat.mean.vec <- rep(NA,nrow(seq.dat.mfc))
ctl.seq.dat.t.vec <- ctl.seq.dat.p.vec <- ctl.seq.dat.mean.vec <- rep(NA,nrow(seq.dat.mfc))

asdSubset <- met.dat.mfc[,"ASD.CTL"]=="ASD"
ctlSubset <- met.dat.mfc[,"ASD.CTL"]=="CTL"

asd.var <- ctl.var <- bartlett.p <- rep(NA,nrow(seq.dat.mfc))

for (i in 1:nrow(seq.dat.mfc)) {
  bout <- bartlett.test(x=c(as.numeric(seq.dat.mfc[i,]),as.numeric(seq.dat.stc[i,])),
                        g=c(as.character(met.dat.mfc[,"ASD.CTL"]),as.character(met.dat.stc[,"ASD.CTL"])))
  bartlett.p[i] <- bout$p.value
  asd.var[i] <- by(data=c(as.numeric(seq.dat.mfc[i,]),as.numeric(seq.dat.stc[i,])),
                   INDICES=c(as.character(met.dat.mfc[,"ASD.CTL"]),as.character(met.dat.stc[,"ASD.CTL"])),
                   FUN=var)[1]
  ctl.var[i] <- by(data=c(as.numeric(seq.dat.mfc[i,]),as.numeric(seq.dat.stc[i,])),
                   INDICES=c(as.character(met.dat.mfc[,"ASD.CTL"]),as.character(met.dat.stc[,"ASD.CTL"])),
                   FUN=var)[2]


  wout <- wilcox.test(as.numeric(seq.dat.mfc[i,asdSubset]),as.numeric(seq.dat.stc[i,asdSubset]),paired=TRUE)
  asd.seq.dat.t.vec[i] <- wout$statistic
  asd.seq.dat.p.vec[i] <- wout$p.value
  asd.seq.dat.mean.vec[i] <- median(as.numeric(seq.dat.mfc[i,asdSubset]))-median(as.numeric(seq.dat.stc[i,asdSubset]))

  wout <- wilcox.test(as.numeric(seq.dat.mfc[i,ctlSubset]),as.numeric(seq.dat.stc[i,ctlSubset]),paired=TRUE)
  ctl.seq.dat.t.vec[i] <- wout$statistic
  ctl.seq.dat.p.vec[i] <- wout$p.value
  ctl.seq.dat.mean.vec[i] <- median(as.numeric(seq.dat.mfc[i,ctlSubset]))-median(as.numeric(seq.dat.stc[i,ctlSubset]))
}

## calculate FDR adjusted p values
asd.qvec <- p.adjust(asd.seq.dat.p.vec,method="BH")
table(asd.qvec < 0.05) ## 34 genes with patterning in ASD

ctl.qvec <- p.adjust(ctl.seq.dat.p.vec,method="BH")
table(ctl.qvec < 0.05) ## 472 genes with patterning in CTL

bart.qvec <- p.adjust(bartlett.p,method="BH")
table(bart.qvec < 0.05) ## many genes have a variance change at q < 0.05

## Plot histograms of the p value distributions
par(mfrow=c(1,3))
hist(ctl.seq.dat.p.vec,ylim=c(0,5000))
hist(asd.seq.dat.p.vec,ylim=c(0,5000))
hist(bartlett.p,ylim=c(0,5000))
hist(bartlett.p[(ctl.qvec<0.05)&!(asd.qvec<0.05)],ylim=c(0,5000),col="red",add=TRUE)
par(mfrow=c(1,2))
plot(ecdf(bartlett.p),main="P value CDF, ACP in red, background in black",cex.main=0.5)
plot(ecdf(bartlett.p[(ctl.qvec<0.05)&!(asd.qvec<0.05)]),add=TRUE,col="red")
boxplot(bartlett.p,bartlett.p[(ctl.qvec<0.05)&!(asd.qvec<0.05)],names=c("background","ACP"),main="Bartlett test P values")
## Genes in the ACP set are less variable between ASD and CTL

## Compute statistics
median(bartlett.p[(ctl.qvec<0.05)&!(asd.qvec<0.05)]) ## median P value = 0.26 for 457 genes in ACP set
median(bartlett.p[!((ctl.qvec<0.05)&!(asd.qvec<0.05))]) ## median P value = 0.19 15853 genes NOT in ACP set
ks.test(bartlett.p[(ctl.qvec<0.05)&!(asd.qvec<0.05)],
        bartlett.p[!((ctl.qvec<0.05)&!(asd.qvec<0.05))])
## In contrast to the results in the publication, there is a significant P value here, but it's that the ACP set distribution shows less variance between regions than non-patterend genes

## Write output
CPAT.dat <- data.frame(ENSGID=rownames(datExpr.CTX.reg),
              ASD.FC.vs.TC.P=asd.seq.dat.p.vec,
              ASD.FC.vs.TC.FDR=asd.qvec,
              CTL.FC.vs.TC.P=ctl.seq.dat.p.vec,
              CTL.FC.vs.TC.FDR=ctl.qvec)
write.csv(CPAT.dat,file="../output/CPATdat.csv")

## TFBS analysis on this ACP set using the pipeline available at https://tfenrichment.semel.ucla.edu/ yielded
##                        Raw_score	    P-value_DB	  P-value_CpG	  P-value_Chro20    GENENAME    MOTIF_ID
##MA0516.1_SP2	          598.902	      0	            0.011	        0.001	            SP2	        MA0516.1
##MA0079.3_SP1	          566.714	      0	            0.01	        0.001	            SP1	        MA0079.3
##MA0599.1_KLF5	          429	          0.001	        0.017	        0.001	            KLF5	      MA0599.1
##MA0472.1_EGR2	          145.851	      0.006	        0.049	        0.005	            EGR2	      MA0472.1
##MA0042.1_FOXI1	        140.215	      0.011	        0.032	        0.014	            FOXI1	      MA0042.1
##MA0136.1_ELF5	          129.226	      0.039	        0.016	        0.009	            ELF5	      MA0136.1
##MA0084.1_SRY	          108.704	      0.04	        0.016	        0.008	            SRY	        MA0084.1
##MA0151.1_ARID3A	        105.037	      0.015	        0.016	        0.014	            ARID3A	    MA0151.1
##MA0442.1_SOX10	        89.0583	      0.046	        0.019	        0.021	            SOX10	      MA0442.1
##MA0087.1_Sox5	          83.1089	      0.019	        0.014	        0.016	            SOX5	      MA0087.1
##MA0074.1_RXRA::VDR	    -6.37162	    0.995	        0.964	        0.967	            RXRA::VDR	  MA0074.1
## Note that SOX5 is in the ACP set and shows enrichment with TFBS.
