########################
## Fabio Morgante
## 7/20/2017
## GBLUP with GO annotation
##
## Modified on 4-24-2018
## Reason: adjust phenotypes, calculate h2 e2 and pool folds and reps to calculate se
##
## Modified on 5-3-2018
## Reason: Get genotype data for only lines with phenotype
## 
## Modified on 11/11/2018
## Reason: add options(stringsAsFactors=FALSE) and remove as.is=T        
##
## Modified on 03/19/2020
## Reason: added options(warn=1) to deal with foreach
##
## Modified on 05/06/2020
## Reason: Use relative paths, use as.numeric(NA),
## 		   set maxcyc=10000, 
##		   improved error/warning handling             
#############################

###Set number of cores for R open
setMKLthreads(1)

###Set global options
options(stringsAsFactors=FALSE)
options(warn=1)

###Set number of cores for foreach
cores <- 10

###Load libraries
library(regress)
library(foreach)
library(doMC)

###Select trait to analyze from command line arguments
args <- commandArgs(TRUE)
trait <- as.numeric(args)

###Source custom functions
source(file="../Software/myRfunctions.R")
load("Data/adjustData.RData")

###Create a function to build the GRM from the genes in a GO, and the GRM with all the other genes
go2GRMs <- function(GO_num=GO_num, GO_data=GO_data, W_data=W_data, phenoLines=phenoLines){
	variants <- unique(unlist(GO_data[[GO_num]]))
    W_data_GO <- W_data[, which(colnames(W_data) %in% variants)]
    GRM_GO <- tcrossprod(W_data_GO)/ncol(W_data_GO)
    GRM_GO <- GRM_GO[phenoLines, phenoLines]
			
    W_data_notGO <- W_data[, -which(colnames(W_data) %in% variants)]
	GRM_notGO <- tcrossprod(W_data_notGO)/ncol(W_data_notGO)
	GRM_notGO <- GRM_notGO[phenoLines, phenoLines]

	return(list(GO=GRM_GO, notGO=GRM_notGO))
}

###Load the data
##Annotations
load("Annotations/go2variants.RData")

##Transcript data
dataF <- as.matrix(read.table("Data/combined_samples_known_novel_fpkm_VR_NoWol_F_line_means_rename_noflag_transp_common.txt",
                              header=T, row.names=1, sep="\t"))
                              
##Centered and scaled genotype matrix (W)
load(file="../DGRPdata/dgrp2_200lines_common_scaled.RData")
                              
##Phenotypes
pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)

###Extract GOs with at least 5 genes
go2variants <- go2variants[sapply(go2variants,length)>=5]

###Select phenotype of interest (discard missing values)
phenoSel <- na.omit(pheno[, c(1, trait)])

###Find lines with both phenotypes and gene expression
common <- intersect(phenoSel[, 1], row.names(dataF))
rm(dataF)

###Select lines with expression in the phenotype data
phenoSel <- phenoSel[which(phenoSel[, 1] %in% common), ] 

###Adjust phenotype for inversion and wolbachia
phenoSel_adj <- adjustPheno(phenoSel, colnames(phenoSel)[2])
yF <- phenoSel_adj[!is.na(phenoSel_adj)]

###Get transcript data for only lines with phenotype
W <- W[names(yF), ]


###Register cores for parallelel execution
registerDoMC(cores)

###Loop across GOs
complete_res <- foreach(i=1:length(go2variants), .combine='rbind', .packages='regress') %dopar% {

	###Compute GRMs and store matrices in a list needed in the prediction function
	Mat <- go2GRMs(GO_num=i, GO_data=go2variants, W_data=W, phenoLines=names(yF))
	
	###Prediction (k-fold CV)
	folds <- 5
	nReps <- 30

	COR.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
	colnames(COR.CV) <- paste('fold=', 1:folds,sep='')

	R2.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
	colnames(R2.CV) <- paste('fold=', 1:folds,sep='')

	Vgo.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
	colnames(Vgo.CV) <- paste('fold=', 1:folds,sep='')

	Vnotgo.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
	colnames(Vnotgo.CV) <- paste('fold=', 1:folds,sep='')

	Ve.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
	colnames(Ve.CV) <- paste('fold=', 1:folds,sep='')
	
	h2go.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
    colnames(h2go.CV) <- paste('fold=', 1:folds,sep='')

	h2notgo.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
    colnames(h2notgo.CV) <- paste('fold=', 1:folds,sep='')
    
    e2.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
    colnames(e2.CV) <- paste('fold=', 1:folds,sep='')



	for(rep in 1:nReps){
  		##Create sets for dividing the lines into folds (if statement needed for when length(yF)/folds) is NOT an integer)
  		set.seed(120+rep)
  		sets1 <- rep(1:folds, (length(yF)/folds))
  		if(length(sets1) != length(yF)){sets2 <- sample(x=1:folds, size=(length(yF)-length(sets1)))
  								  sets <- c(sets1, sets2)} else {sets <- sets1}
  		sets <- sets[order(runif(length(yF)))]
    
  		for(fold in 1:folds){
    		##Assign NAs to lines in the test set
    		yNa <- yF
    		whichNa <- which(sets==fold)
    		yNa[whichNa] <- NA
    
    		##GO-TBLUP
    		GRM_GO <- Mat[[1]]
    		GRM_notGO <- Mat[[2]]
    		
    		fm <- tryCatch(regress(yNa ~ 1, ~GRM_GO+GRM_notGO, pos=c(T, T, T), maxcyc = 10000), error=function(e){message(paste0("Error in GO = ", i, ", CV = ", rep, ", fold = ", fold, ":\n", e)); return(NULL)}, 
                   warning=function(w){message(paste0("Warning in GO = ", i, ", CV = ", rep, ", fold = ", fold, ":\n", w)); return(NULL)})
    		
    		##when fm is NULL assign NA to variance partition and prediction accuracy
    		if(!is.null(fm) && fm$converged){
    
    			##Phenotype prediction
    			yhat <- pred_y(fm, Mat, whichNa)
    			yhatNa <- yhat$VDtot
    
    			##Variance partition in the training set
    			Vgo.CV[rep, fold] <- fm$sigma[1]
    			Vnotgo.CV[rep, fold] <- fm$sigma[2]
    			Ve.CV[rep, fold] <- fm$sigma[3]
    		    h2go.CV[rep, fold] <- fm$sigma[1]/(fm$sigma[1]+fm$sigma[2]+fm$sigma[3])
    		    h2notgo.CV[rep, fold] <- fm$sigma[2]/(fm$sigma[1]+fm$sigma[2]+fm$sigma[3])
                e2.CV[rep, fold] <- fm$sigma[3]/(fm$sigma[1]+fm$sigma[2]+fm$sigma[3])

    
    			##Compute correlation and R2 in the test set
    			COR.CV[rep, fold] <- cor(yF[whichNa], yhatNa)
    			R2.CV[rep, fold] <- (cor(yF[whichNa], yhatNa))^2
    			
    			} else {
    			
    			##Variance partition in the training set
    			Vgo.CV[rep, fold] <- NA
    			Vnotgo.CV[rep, fold] <- NA
    			Ve.CV[rep, fold] <- NA
    		    h2go.CV[rep, fold] <- NA
    		    h2notgo.CV[rep, fold] <- NA
                e2.CV[rep, fold] <- NA
    
    			##Compute correlation and R2 in the test set
    			COR.CV[rep, fold] <- NA
    			R2.CV[rep, fold] <- NA
    			}
    			
    		#print(paste("finished fold =", fold))

  		}

  		#print(paste("finished CV =", rep))
	}

	###Compute mean and se across CV replicates
	meanVgo <- mean(Vgo.CV, na.rm=T)
	se_meanVgo <- sd(Vgo.CV, na.rm=T)/sqrt(sum(!is.na(Vgo.CV)))
	meanVnotgo <- mean(Vnotgo.CV, na.rm=T)
	se_meanVnotgo <- sd(Vnotgo.CV, na.rm=T)/sqrt(sum(!is.na(Vnotgo.CV)))
	meanVe <- mean(Ve.CV, na.rm=T)
	se_meanVe <- sd(Ve.CV, na.rm=T)/sqrt(sum(!is.na(Ve.CV)))
	meanCOR <- mean(COR.CV, na.rm=T)
	se_meanCOR <- sd(COR.CV, na.rm=T)/sqrt(sum(!is.na(COR.CV)))
	meanR2 <- mean(R2.CV, na.rm=T)
	se_meanR2 <- sd(R2.CV, na.rm=T)/sqrt(sum(!is.na(R2.CV)))
	meanh2go <- mean(h2go.CV, na.rm=T)
    se_meanh2go <- sd(h2go.CV, na.rm=T)/sqrt(sum(!is.na(h2go.CV)))
	meanh2notgo <- mean(h2notgo.CV, na.rm=T)
    se_meanh2notgo <- sd(h2notgo.CV, na.rm=T)/sqrt(sum(!is.na(h2notgo.CV)))
    meane2 <- mean(e2.CV, na.rm=T)
    se_meane2 <- sd(e2.CV, na.rm=T)/sqrt(sum(!is.na(e2.CV)))
    nreps_no_NA <- sum(!is.na(R2.CV))


	###Put the results together
	Results <- matrix(c(meanVgo, se_meanVgo, meanVnotgo, se_meanVnotgo, meanVe, se_meanVe, meanCOR, se_meanCOR, meanR2, se_meanR2, meanh2go, se_meanh2go, meanh2notgo, se_meanh2notgo, meane2, se_meane2, nreps_no_NA), ncol=17, byrow=T)		
	GO <- names(go2variants)[i]
	Results_final <- data.frame(GO, Results)
	Results_final
	
	
	#print(paste("finished GO =", i))
}

colnames(complete_res) <- c("GO", "mean Vgo", "se Vgo", "mean Vnotgo", "se Vnotgo", "mean Ve", "se Ve", "mean r", "se r", "mean r2", "se r2", "mean h2go", "se h2go", "mean h2notgo", "se h2notgo", "mean e2", "se e2", "nreps")

write.table(complete_res, paste('Results/Pred_SNPs_GO_regress_adjPheno_', colnames(pheno)[trait],'_females.csv', sep=""), 
		  col.names=T, row.names=F, quote=F, sep=",")





