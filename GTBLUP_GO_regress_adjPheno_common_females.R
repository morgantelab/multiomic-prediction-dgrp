########################
## Fabio Morgante
## 7/13/2017
## GTBLUP with GO annotation
##
## Modified on 4-24-2018
## Reason: adjust phenotypes, calculate h2 e2 and pool folds and reps to calculate se
##
## Modified on 5-3-2018
## Reason: Get genotype and transcript data for only lines with phenotype
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

###Create a function to build the TRM from the genes in a GO, and the TRM with all the other genes
go2TRMs <- function(GO_num=GO_num, GO_data=GO_data, W_transc_data=W_transc_data, phenoLines=phenoLines){
	genes <- GO_data[[GO_num]]
    W_transc_data_GO <- W_transc_data[, which(colnames(W_transc_data) %in% genes)]
    TRM_GO <- tcrossprod(W_transc_data_GO)/ncol(W_transc_data_GO)
    TRM_GO <- TRM_GO[phenoLines, phenoLines]
			
    W_transc_data_notGO <- W_transc_data[, -which(colnames(W_transc_data) %in% genes)]
	TRM_notGO <- tcrossprod(W_transc_data_notGO)/ncol(W_transc_data_notGO)
	TRM_notGO <- TRM_notGO[phenoLines, phenoLines]

	return(list(GO=TRM_GO, notGO=TRM_notGO))
}

###Load the data
##Annotations
#Transcripts
load("Annotations/go2fb_F.RData")
#SNPs
load("Annotations/go2variants.RData")

##Centered and scaled genotype matrix (W)
load(file="../DGRPdata/dgrp2_200lines_common_scaled.RData")

##Centered and scaled transcript matrix (W_F)
load(file="Matrices/dgrp2_scaled_transcripts_common_females.RData")
                              
##Phenotypes
pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)

###Select phenotype of interest (discard missing values)
phenoSel <- na.omit(pheno[, c(1, trait)])

###Find lines with both phenotypes and gene expression
common <- intersect(phenoSel[, 1], row.names(W_F))

###Select lines with expression in the phenotype data
phenoSel <- phenoSel[which(phenoSel[, 1] %in% common), ] 

###Adjust phenotype for inversion and wolbachia
phenoSel_adj <- adjustPheno(phenoSel, colnames(phenoSel)[2])
yF <- phenoSel_adj[!is.na(phenoSel_adj)]

###Get genotype data for only lines with phenotype
W <- W[names(yF), ]

###Get transcript data for only lines with phenotype
W_F <- W_F[names(yF), ]

###Extract GOs-transcripts with at least 5 genes
go2fb_F <- go2fb_F[sapply(go2fb_F,length)>=5]

###Extract GOs-SNPs with at least 5 genes
go2variants <- go2variants[sapply(go2variants,length)>=5]

###Extract common GOs
common_GOs <- intersect(names(go2fb_F), names(go2variants))
go2fb_F <- go2fb_F[common_GOs]
go2variants <- go2variants[common_GOs]


###Register cores for parallelel execution
registerDoMC(cores)

###Loop across GOs
complete_res <- foreach(i=1:length(go2fb_F), .combine='rbind', .packages='regress') %dopar% {

	###Compute TRMs and store matrices in a list needed in the prediction function
	Mat1 <- go2TRMs(GO_num=i, GO_data=go2fb_F, W_transc_data=W_F, phenoLines=names(yF))
	Mat2 <- go2GRMs(GO_num=i, GO_data=go2variants, W_data=W, phenoLines=names(yF))
	Mat <- c(Mat1, Mat2)
	rm(Mat1, Mat2)
	
	###Prediction (k-fold CV)
	folds <- 5
	nReps <- 30

	COR.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
	colnames(COR.CV) <- paste('fold=', 1:folds,sep='')

	R2.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
	colnames(R2.CV) <- paste('fold=', 1:folds,sep='')

	Vt_go.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
	colnames(Vt_go.CV) <- paste('fold=', 1:folds,sep='')

	Vt_notgo.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
	colnames(Vt_notgo.CV) <- paste('fold=', 1:folds,sep='')
	
	Vg_go.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
	colnames(Vg_go.CV) <- paste('fold=', 1:folds,sep='')

	Vg_notgo.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
	colnames(Vg_notgo.CV) <- paste('fold=', 1:folds,sep='')

	Ve.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
	colnames(Ve.CV) <- paste('fold=', 1:folds,sep='')

	h2t_go.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
    colnames(h2t_go.CV) <- paste('fold=', 1:folds,sep='')

	h2t_notgo.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
    colnames(h2t_notgo.CV) <- paste('fold=', 1:folds,sep='')
    
	h2g_go.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
    colnames(h2g_go.CV) <- paste('fold=', 1:folds,sep='')

	h2g_notgo.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
    colnames(h2g_notgo.CV) <- paste('fold=', 1:folds,sep='')
    
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
    
    		##GO-GTBLUP
    		TRM_GO <- Mat[[1]]
    		TRM_notGO <- Mat[[2]]
    		GRM_GO <- Mat[[3]]
    		GRM_notGO <- Mat[[4]]

    		
    		fm <- tryCatch(regress(yNa ~ 1, ~TRM_GO+TRM_notGO+GRM_GO+GRM_notGO, pos=c(T, T, T, T, T), maxcyc = 10000), error=function(e){message(paste0("Error in GO = ", i, ", CV = ", rep, ", fold = ", fold, ":\n", e)); return(NULL)}, 
                   warning=function(w){message(paste0("Warning in GO = ", i, ", CV = ", rep, ", fold = ", fold, ":\n", w)); return(NULL)})
    		
    		##when fm is NULL assign NA to variance partition and prediction accuracy
    		if(!is.null(fm) && fm$converged){
    
    			##Phenotype prediction
    			yhat <- pred_y(fm, Mat, whichNa)
    			yhatNa <- yhat$VDtot
    
    			##Variance partition in the training set
    			Vt_go.CV[rep, fold] <- fm$sigma[1]
    			Vt_notgo.CV[rep, fold] <- fm$sigma[2]
    			Vg_go.CV[rep, fold] <- fm$sigma[3]
    			Vg_notgo.CV[rep, fold] <- fm$sigma[4]
    			Ve.CV[rep, fold] <- fm$sigma[5]
    		    h2t_go.CV[rep, fold] <- fm$sigma[1]/(fm$sigma[1]+fm$sigma[2]+fm$sigma[3]+fm$sigma[4]+fm$sigma[5])
    		    h2t_notgo.CV[rep, fold] <- fm$sigma[2]/(fm$sigma[1]+fm$sigma[2]+fm$sigma[3]+fm$sigma[4]+fm$sigma[5])
    		    h2g_go.CV[rep, fold] <- fm$sigma[3]/(fm$sigma[1]+fm$sigma[2]+fm$sigma[3]+fm$sigma[4]+fm$sigma[5])
    		    h2g_notgo.CV[rep, fold] <- fm$sigma[4]/(fm$sigma[1]+fm$sigma[2]+fm$sigma[3]+fm$sigma[4]+fm$sigma[5])
    		    e2.CV[rep, fold] <- fm$sigma[5]/(fm$sigma[1]+fm$sigma[2]+fm$sigma[3]+fm$sigma[4]+fm$sigma[5])

    
    			##Compute correlation and R2 in the test set
    			COR.CV[rep, fold] <- cor(yF[whichNa], yhatNa)
    			R2.CV[rep, fold] <- (cor(yF[whichNa], yhatNa))^2
    			
    			} else {
    			
    			##Variance partition in the training set
    			Vt_go.CV[rep, fold] <- NA
    			Vt_notgo.CV[rep, fold] <- NA
    			Vg_go.CV[rep, fold] <- NA
    			Vg_notgo.CV[rep, fold] <- NA
    			Ve.CV[rep, fold] <- NA
    		    h2t_go.CV[rep, fold] <- NA
    		    h2t_notgo.CV[rep, fold] <- NA
    		    h2g_go.CV[rep, fold] <- NA
    		    h2g_notgo.CV[rep, fold] <- NA
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
	meanVt_go <- mean(Vt_go.CV, na.rm=T)
	se_meanVt_go <- sd(Vt_go.CV, na.rm=T)/sqrt(sum(!is.na(Vt_go.CV)))
	meanVt_notgo <- mean(Vt_notgo.CV, na.rm=T)
	se_meanVt_notgo <- sd(Vt_notgo.CV, na.rm=T)/sqrt(sum(!is.na(Vt_notgo.CV)))
	meanVg_go <- mean(Vg_go.CV, na.rm=T)
	se_meanVg_go <- sd(Vg_go.CV, na.rm=T)/sqrt(sum(!is.na(Vg_go.CV)))
	meanVg_notgo <- mean(Vg_notgo.CV, na.rm=T)
	se_meanVg_notgo <- sd(Vg_notgo.CV, na.rm=T)/sqrt(sum(!is.na(Vg_notgo.CV)))
	meanVe <- mean(Ve.CV, na.rm=T)
	se_meanVe <- sd(Ve.CV, na.rm=T)/sqrt(sum(!is.na(Ve.CV)))
	meanCOR <- mean(COR.CV, na.rm=T)
	se_meanCOR <- sd(COR.CV, na.rm=T)/sqrt(sum(!is.na(COR.CV)))
	meanR2 <- mean(R2.CV, na.rm=T)
	se_meanR2 <- sd(R2.CV, na.rm=T)/sqrt(sum(!is.na(R2.CV)))
	meanh2t_go <- mean(h2t_go.CV, na.rm=T)
    se_meanh2t_go <- sd(h2t_go.CV, na.rm=T)/sqrt(sum(!is.na(h2t_go.CV)))
	meanh2t_notgo <- mean(h2t_notgo.CV, na.rm=T)
    se_meanh2t_notgo <- sd(h2t_notgo.CV, na.rm=T)/sqrt(sum(!is.na(h2t_notgo.CV)))
	meanh2g_go <- mean(h2g_go.CV, na.rm=T)
    se_meanh2g_go <- sd(h2g_go.CV, na.rm=T)/sqrt(sum(!is.na(h2g_go.CV)))
	meanh2g_notgo <- mean(h2g_notgo.CV, na.rm=T)
    se_meanh2g_notgo <- sd(h2g_notgo.CV, na.rm=T)/sqrt(sum(!is.na(h2g_notgo.CV)))
    meane2 <- mean(e2.CV, na.rm=T)
    se_meane2 <- sd(e2.CV, na.rm=T)/sqrt(sum(!is.na(e2.CV)))
    nreps_no_NA <- sum(!is.na(R2.CV))

	

	###Put the results together
	Results <- matrix(c(meanVt_go, se_meanVt_go, meanVt_notgo, se_meanVt_notgo, meanVg_go, se_meanVg_go, meanVg_notgo, se_meanVg_notgo, meanVe, se_meanVe, meanCOR, se_meanCOR, meanR2, se_meanR2,
						meanh2t_go, se_meanh2t_go, meanh2t_notgo, se_meanh2t_notgo, meanh2g_go, se_meanh2g_go, meanh2g_notgo, se_meanh2g_notgo, meane2, se_meane2, nreps_no_NA), ncol=25, byrow=T)		
	GO <- names(go2fb_F)[i]
	Results_final <- data.frame(GO, Results)
	Results_final
	
	
	#print(paste("finished GO =", i))
}

colnames(complete_res) <- c("GO", "mean Vt_go", "se Vt_go", "mean Vt_notgo", "se Vt_notgo", "mean Vg_go", "se Vg_go", "mean Vg_notgo", "se Vg_notgo", "mean Ve", "se Ve", "mean r", "se r", "mean r2", "se r2",
							"mean h2t_go", "se h2t_go", "mean h2t_notgo", "se h2t_notgo", "mean h2g_go", "se h2g_go", "mean h2g_notgo", "se h2g_notgo", "mean e2", "se e2", "nreps")

write.table(complete_res, paste('Results/Pred_SNPs+Expression_GO_regress_adjPheno_common_', colnames(pheno)[trait],'_females.csv', sep=""), 
		  col.names=T, row.names=F, quote=F, sep=",")





