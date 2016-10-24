## Code related to "Genome-wide changes in lncRNA, alternative splicing, and cortical patterning in autism"
## This script will run differential gene expression analyses for CTX and CBL

rm(list=ls())
options(stringsAsFactors=FALSE)
library(biomaRt)
library(nlme)

runlme <- function(thisdat,expression) {
  lme1 <- eval(parse(text=expression));
  ##Get the summary of the model
  smodel = summary(lme1);
  return(smodel)
}

## Data structure to keep all of these results
DGEtables <- vector(mode="list",length=12)
nameVec <- c("ASDvCTLinCTX","ASDvCTLinCBL",
             "dup15qvCTLinCTX","dup15qvCTLinCBL",
             "dup15qvCTLinCTX_rep","dup15qvCTLinCBL_rep",
             "AllASDvsCTLinCTX","AllASDvsCTLinCBL",
             "ASDvCTLinCTX_oldonly","ASDvCTLinCBL_oldonly",
             "ASDvCTLinCTX_newonly","ASDvCTLinCBL_newonly")
for (i in 1:length(DGEtables)) {
  DGEtables[[i]]$name <- nameVec[i]
  print(nameVec[i])
}

## RNA-seq data
load("../data/modified/NormalizedExpressionData.Rdata") ## Use the normalized FPKM data, since we utilize covariates in the LMM

for (comp in c(1:12)) {
    if (comp == 1) { ## Idiopathic ASD vs CTL comparison in CTX
      keepsamps <- datMeta.CTX[,"ASD.vs.CTL..CTX"]
      datExpr.sub <- datExpr.CTX[,keepsamps]
      datMeta.sub <- datMeta.CTX[keepsamps,]
    } else if (comp == 2) { ## Idiopathic ASD vs CTL comparison in CBL
      keepsamps <- datMeta.CBL[,"ASD.vs.CTL..CB"]
      datExpr.sub <- datExpr.CBL[,keepsamps]
      datMeta.sub <- datMeta.CBL[keepsamps,]
    } else if (comp == 3) { ## dup15q vs CTL (independent controls) analysis in CTX
      keepsamps <- datMeta.CTX[,"dup15q.vs.CTL..CTX"]
      datExpr.sub <- datExpr.CTX[,keepsamps]
      datMeta.sub <- datMeta.CTX[keepsamps,]
    } else if (comp == 4) { ## dup15q vs CTL (independent controls) analysis in CBL
      keepsamps <- datMeta.CBL[,"dup15q.vs.CTL..CB"]
      datExpr.sub <- datExpr.CBL[,keepsamps]
      datMeta.sub <- datMeta.CBL[keepsamps,]
    } else if (comp == 5) { ## dup15q vs CTL (controls from comparison 1 above) analysis in CTX
      keepsamps1 <- datMeta.CTX[,"dup15q.vs.CTL..CTX"]&datMeta.CTX[,"ASD.CTL"]=="ASD"&datMeta.CTX[,"SeqBatch"]!="batch3"
      keepsamps2 <- datMeta.CTX[,"ASD.vs.CTL..CTX"]&datMeta.CTX[,"ASD.CTL"]=="CTL"
      keepsamps <- keepsamps1|keepsamps2
      datExpr.sub <- datExpr.CTX[,keepsamps]
      datMeta.sub <- datMeta.CTX[keepsamps,]
    } else if (comp == 6) { ## dup15q vs CTL (contorls from comparison 2 above) analysis in CBL
      keepsamps1 <- datMeta.CBL[,"dup15q.vs.CTL..CB"]&datMeta.CBL[,"ASD.CTL"]=="ASD"&datMeta.CBL[,"SeqBatch"]!="batch3"
      keepsamps2 <- datMeta.CBL[,"ASD.vs.CTL..CB"]&datMeta.CBL[,"ASD.CTL"]=="CTL"
      keepsamps <- keepsamps1|keepsamps2
      datExpr.sub <- datExpr.CBL[,keepsamps]
      datMeta.sub <- datMeta.CBL[keepsamps,]
    } else if (comp == 7) { ## All ASD vs CTL (all samples used in network analysis) analysis in CTX
      keepsamps <- datMeta.CTX[,"Network.analysis..CTX"]
      datExpr.sub <- datExpr.CTX[,keepsamps]
      datMeta.sub <- datMeta.CTX[keepsamps,]
    } else if (comp == 8) { ## All ASD vs CTL (all samples used in network analysis) analysis in CBL
      keepsamps <- datMeta.CBL[,"Network.Analysis..CB"]
      datExpr.sub <- datExpr.CBL[,keepsamps]
      datMeta.sub <- datMeta.CBL[keepsamps,]
    } else if (comp == 9) { ## Idiopathic ASD vs CTL comparison in CTX, old samples from matched set only
      keepsamps <- datMeta.CTX[,"ASD.vs.CTL..CTX"]&datMeta.CTX[,"Voineagu.et.al..CTX.sample"]
      datExpr.sub <- datExpr.CTX[,keepsamps]
      datMeta.sub <- datMeta.CTX[keepsamps,]
    } else if (comp == 10) { ## Idiopathic ASD vs CTL comparison in CBL, old samples from matched set only
      keepsamps <- datMeta.CBL[,"ASD.vs.CTL..CB"]&datMeta.CBL[,"Voineagu.et.al..CB.sample"]
      datExpr.sub <- datExpr.CBL[,keepsamps]
      datMeta.sub <- datMeta.CBL[keepsamps,]
    } else if (comp == 11) { ## Idiopathic ASD vs CTL comparison in CTX, new samples from matched set only
      keepsamps <- datMeta.CTX[,"ASD.vs.CTL..CTX"]&!datMeta.CTX[,"Voineagu.et.al..CTX.sample"]
      datExpr.sub <- datExpr.CTX[,keepsamps]
      datMeta.sub <- datMeta.CTX[keepsamps,]
    } else if (comp == 12) { ## Idiopathic ASD vs CTL comparison in CBL, new samples from matched set only
      keepsamps <- datMeta.CBL[,"ASD.vs.CTL..CB"]&!datMeta.CBL[,"Voineagu.et.al..CB.sample"]
      datExpr.sub <- datExpr.CBL[,keepsamps]
      datMeta.sub <- datMeta.CBL[keepsamps,]
    }}

    print(DGEtables[[comp]]$name)
    datMeta <- datMeta.sub
    datExpr <- datExpr.sub

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
      bank <- batch2 <- batch3 <- rep(0,length(batch2)) ## For the dup15q comparison, we don't use batch or brain bank as covariates because they are confounded with diagnosis. We provide an analysis against other control samples to demonstrate that the results are identical when using the set of CTL samples used in the idiopathic ASD analysis.
    }

    varnames <- c("condition","age","sex","region","batch2","batch3","RIN","bank","seqStatPC1","seqStatPC2")

    varmat <- cbind(condition, age, sex, region,batch2,batch3,RIN,bank,seqStatPC1,seqStatPC2)
    varkeep <- rep(TRUE,length(varnames))
    for (i in 1:length(varnames)) {
        if (length(table(varmat[,i])) == 1) {
            varkeep[i] <- FALSE
        }
    }

    Bmat <- SEmat <- Pmat <- matrix(NA,nrow=nrow(datExpr),ncol=length(varnames))[,varkeep]
    rownames(Bmat) <- rownames(SEmat) <- rownames(Pmat) <- rownames(datExpr)
    colnames(Bmat) <- paste("beta",varnames,sep=".")[varkeep]
    colnames(SEmat) <- paste("SE",varnames,sep=".")[varkeep]
    colnames(Pmat) <- paste("p",varnames,sep=".")[varkeep]

    for (i in 1:nrow(datExpr)) {
        if (i %% 1000 == 0) {cat(paste("On gene ",i,"\n",sep=""))}
        thisExpr <- as.numeric(datExpr[i,])
        expression <- paste("lme(thisExpr ~ ",paste(colnames(varmat)[varkeep],collapse=" + "),", rand = ~1|biolrep, data = thisdat)",sep="")
        designmatrix <- data.frame(thisExpr, varmat[,varkeep], biolrep)

        lme1.out <- try(runlme(designmatrix,expression),silent=F);

        if (substr(lme1.out[1],1,5)!="Error") {
            tabOut <- lme1.out$tTable
            Bmat[i,] <- tabOut[-c(1),"Value"]
            SEmat[i,] <- tabOut[-c(1),"Std.Error"]
            Pmat[i,] <- tabOut[-c(1),"p-value"]
        } else {
            cat('Error in LME of Gene',rownames(datExpr)[i],"id",'\n')
            cat('Setting P-value=NA,Beta value=NA, and SE=NA\n');
            Bmat[i,] <- SEmat[i,] <- Pmat[i,] <- NA;
        }
    }

    ## Get gene annotation information
    getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","band","gene_biotype")
    mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="sep2013.archive.ensembl.org") ## Gencode v18, otherwise use ##mart <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")

    geneAnno1 <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values=rownames(Bmat),mart=mart)
    geneAnno2 <- geneAnno1[match(rownames(Bmat),geneAnno1[,1]),]
    datSignifAnno <- cbind(Bmat,SEmat,Pmat,geneAnno2)[!is.na(Bmat[,1]),]

    ## Save the data
    DGEtables[[comp]]$metadata <- datMeta
    DGEtables[[comp]]$datExpr <- datExpr
    DGEtables[[comp]]$datSignifAnno <- datSignifAnno
}
save(DGEtables,file="../output/DGEtables.Rdata")
## Note there are minor differences in results with what is published due to differences in initial normalization.
