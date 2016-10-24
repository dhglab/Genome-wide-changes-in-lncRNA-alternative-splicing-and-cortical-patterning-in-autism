## Code related to "Genome-wide changes in lncRNA, alternative splicing, and cortical patterning in autism"
## This script will run differential alternative splicing analyses for CTX and CBL

rm(list=ls())
options(stringsAsFactors=FALSE)
library(nlme)

runlme <- function(thisdat) {
  lme1 <- eval(parse(text=expression));    
  ##Get the summary of the model
  smodel = summary(lme1);
  return(smodel)
}

##### Major differential alternative splicing (DS) comparisons
## Data structure to keep all of these results
DStables <- spliceDEXolego <- vector(mode="list",length=6)
nameVec <- c("ASDvCTLinCTX","ASDvCTLinCBL",
             "dup15qvCTLinCTX","dup15qvCTLinCBL",
             "AllASDvsCTLinCTX","AllASDvsCTLinCBL")
for (i in 1:length(DStables)) {
    DStables[[i]]$name <- nameVec[i]
}

## Load splicing data, as percent spliced in (PSI) levels processed by MATS.
load("../data/provided/ASD_PSI_Data.Rdata")

for (comp in c(1:6)) {
  if (comp == 1) { ## Idiopathic ASD vs CTL comparison in CTX
    keepsamps <- datMeta.CTX[,"ASD.vs.CTL..CTX"]
    datPSI.sub <- datPSI.CTX[,keepsamps]
    datMeta.sub <- datMeta.CTX[keepsamps,]
  } else if (comp == 2) { ## Idiopathic ASD vs CTL comparison in CBL
    keepsamps <- datMeta.CBL[,"ASD.vs.CTL..CB"]
    datPSI.sub <- datPSI.CBL[,keepsamps]
    datMeta.sub <- datMeta.CBL[keepsamps,]      
  } else if (comp == 3) { ## dup15q vs CTL (independent controls) analysis in CTX
    keepsamps <- datMeta.CTX[,"dup15q.vs.CTL..CTX"]
    datPSI.sub <- datPSI.CTX[,keepsamps]
    datMeta.sub <- datMeta.CTX[keepsamps,]      
  } else if (comp == 4) { ## dup15q vs CTL (independent controls) analysis in CBL
    keepsamps <- datMeta.CBL[,"dup15q.vs.CTL..CB"]
    datPSI.sub <- datPSI.CBL[,keepsamps]
    datMeta.sub <- datMeta.CBL[keepsamps,]      
  } else if (comp == 5) { ## All ASD vs CTL (all samples used in network analysis) analysis in CTX
    keepsamps <- datMeta.CTX[,"Network.analysis..CTX"]
    datPSI.sub <- datPSI.CTX[,keepsamps]
    datMeta.sub <- datMeta.CTX[keepsamps,]      
  } else if (comp == 6) { ## All ASD vs CTL (all samples used in network analysis) analysis in CBL
    keepsamps <- datMeta.CBL[,"Network.Analysis..CB"]
    datPSI.sub <- datPSI.CBL[,keepsamps]
    datMeta.sub <- datMeta.CBL[keepsamps,]      
  }  

    ## First, the MATS PSI values
    print(DStables[[comp]]$name)
    datPSI <- datPSI.sub
    datMeta <- datMeta.sub
    
    ## Calculate the seqSVs
    ## Do a PCA of the sequencing statistics of the full sample
    datSeq <- data.matrix(datMeta[,c(25:43)])
    datSeq <- datSeq[,c(6:19,c(5))]
    datSeqNorm <- t(scale(datSeq,scale=F))
    PC.datSeq <- prcomp(datSeqNorm);
    varexp <- (PC.datSeq$sdev)^2 / sum(PC.datSeq$sdev^2)
    print(varexp[1:5])
    topPC.datSeq <- PC.datSeq$rotation[,1:2]; ## Explains 99% of variance in datSeq
    colnames(topPC.datSeq) <- c("SeqSV1 - Depth","SeqSV2 - GC/Length")
    
    ## Get the metadata
    biolrep <- as.numeric(as.factor(datMeta[,"BrainID"]))
    condition <- 2-as.numeric(as.factor(datMeta[,"ASD.CTL"]))
    age <- as.numeric(datMeta[,"Age"])
    sex <- as.numeric(as.factor(datMeta[,"Sex"]))-1
    region <- as.numeric(as.factor(datMeta[,"RegionID"]))-1
    batch2 <- as.numeric(datMeta[,"SeqBatch"]=="batch2")
    batch3 <- as.numeric(datMeta[,"SeqBatch"]=="batch3")
    RIN <- as.numeric(datMeta[,"RIN"])
    bank <- as.numeric(as.factor(datMeta[,"BrainBank"]))-1
    seqStatPC1 <- as.numeric(topPC.datSeq[,1])
    seqStatPC2 <- as.numeric(topPC.datSeq[,2])
    
    if (comp == 3 | comp == 4) {
      bank <- batch2 <- batch3 <- rep(0,length(batch2)) ## For the dup15q splicing comparison, we only use condition and region as covariates because they are either confounded with diagnosis or introduce singularities into the LMM calculation. We provide an analysis against other control samples to demonstrate that the results are identical when using the set of CTL samples used in the idiopathic ASD analysis.
    }
    
    varnames <- c("condition","age","sex","region","batch2","batch3","RIN","bank","seqStatPC1","seqStatPC2")
    
    varmat <- cbind(condition, age, sex, region,batch2,batch3,RIN,bank,seqStatPC1,seqStatPC2)
    varkeep <- rep(TRUE,length(varnames))
    for (i in 1:length(varnames)) {
      if (length(table(varmat[,i])) == 1) {
        varkeep[i] <- FALSE
      }
    }
    
    Bmat <- SEmat <- Pmat <- matrix(NA,nrow=nrow(datPSI),ncol=length(varnames))[,varkeep]
    rownames(Bmat) <- rownames(SEmat) <- rownames(Pmat) <- rownames(datPSI)
    colnames(Bmat) <- paste("beta",varnames,sep=".")[varkeep]
    colnames(SEmat) <- paste("SE",varnames,sep=".")[varkeep]
    colnames(Pmat) <- paste("p",varnames,sep=".")[varkeep]
    Nsamp <- matrix(NA,nrow=nrow(datPSI),ncol=ncol(datPSI))
   
    for (i in 1:nrow(datPSI)) {
        if (i %% 1000 == 0) {cat(paste("On event ",i,"\n",sep=""))}
        thisPSI <- as.numeric(datPSI[i,])
        keep <- !is.na(thisPSI) ## There exist cases where the PSI is called as NA in a sample
        Nsamp[i,] <- keep ## Keep track of how many PSI values are not NA for each event

        expression <- paste("lme(thisPSI ~ ",paste(colnames(varmat)[varkeep],collapse=" + "),", rand = ~1|biolrep, data = thisdat)",sep="")
        designmatrix <- data.frame(thisPSI, varmat[,varkeep], biolrep)
        lme1.out <- try(runlme(designmatrix),silent=T);
        
        if (substr(lme1.out[1],1,5)!="Error") {
            tabOut <- lme1.out$tTable
            Bmat[i,] <- tabOut[-c(1),"Value"]
            SEmat[i,] <- tabOut[-c(1),"Std.Error"]
            Pmat[i,] <- tabOut[-c(1),"p-value"]
        } else {
            ##cat('Error in LME of event ',rownames(datPSI)[i],"id",'\n')
            ##cat('Setting P-value=NA,Beta value=NA, and SE=NA\n');
            Bmat[i,] <- SEmat[i,] <- Pmat[i,] <- NA;
        }
    }

    ## Get p-values and q-values
    DStables[[comp]]$metadata <- datMeta
    DStables[[comp]]$datPSI <- datPSI
    DStables[[comp]]$datSignifAnno <- cbind(Bmat,SEmat,Pmat)
    DStables[[comp]]$keptSamples <- Nsamp
}
save(DStables,file = "../output/DStables.Rdata")
## Note there are minor differences in results with what is published for CB due to the use of the linear mixed model above while in the study, we removed technical replicates in CB and used a linear model
## CTX DS results are correlated at R = 1
## CBL DS results are correlated at R > 0.976
