########################
## Fabio Morgante
## 2/3/2017
## BLUP with expression data
##
## Modified on 4-5-2018
## Reason: adjust phenotypes, use updated GRM
##
## Modified on 4-6-2018
## Reason: calculate h2 and pool folds and reps to calculate se
##
## Modified on 4-17-2018
## Reason: calculate e2
## 
## Modified on 11/11/2018
## Reason: add options(stringsAsFactors=FALSE) and remove as.is=T        
##
## Modified on 05/06/2020
## Reason: Use relative paths, use as.numeric(NA),
## 		   set maxcyc=10000, add options(warn=1)             
##
## Modified on 05/13/2020
## Reason: Add error/warning handling in REML        
#############################

###Set number of cores for R open
setMKLthreads(2)

###Set global options
options(stringsAsFactors=FALSE)
options(warn=1)

###Load library
library(regress)

###Select trait to analyze from command line arguments
args <- commandArgs(TRUE)
trait <- as.numeric(args)


###############
####Females####
###############

#######
##TRM##
#######

###Source custom functions
source(file="../Software/myRfunctions.R")
load("Data/adjustData.RData")

###Load the data
#TRM
TRM <- as.matrix(read.table("Matrices/TRM_common_females.txt", header=TRUE, sep=" ", row.names=1))

#Phenotypes
pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)


###Select phenotype of interest (discard missing values)
phenoSel <- na.omit(pheno[, c(1, trait)])

###Find lines with both phenotypes and gene expression
common <- intersect(phenoSel[, 1], row.names(TRM))

###Select lines with expression in the phenotype data
phenoSel <- phenoSel[which(phenoSel[, 1] %in% common), ] 

###Adjust phenotype for inversion and wolbachia
phenoSel_adj <- adjustPheno(phenoSel, colnames(phenoSel)[2])
yF <- phenoSel_adj[!is.na(phenoSel_adj)]

###Get TRM for only lines with phenotype
TRM <- TRM[names(yF), names(yF)]

###Store matrices in a list needed in the prediction function
Mat <- as.list(NULL)
Mat[[1]] <- TRM

###Prediction (k-fold CV)
folds <- 5
nReps <- 30

COR.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(COR.CV) <- paste('fold=', 1:folds,sep='')

R2.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(R2.CV) <- paste('fold=', 1:folds,sep='')

h2t.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(h2t.CV) <- paste('fold=', 1:folds,sep='')

Vt.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(Vt.CV) <- paste('fold=', 1:folds,sep='')

Ve.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(Ve.CV) <- paste('fold=', 1:folds,sep='')

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
    
    ##TBLUP
    fm <- tryCatch(regress(yNa ~ 1, ~TRM, pos=c(T, T), maxcyc=10000), error=function(e){message(paste0("Error in CV = ", rep, ", fold = ", fold, ":\n", e)); return(NULL)}, 
                   warning=function(w){message(paste0("Warning in CV = ", rep, ", fold = ", fold, ":\n", w)); return(NULL)})
                   
    if(!is.null(fm) && fm$converged){   
    	##Phenotype prediction
    	yhat <- pred_y(fm, Mat, whichNa)
    	yhatNa <- yhat$VDtot  
    
    	##Variance partition in the training set
    	Vt.CV[rep, fold] <- fm$sigma[1]
    	Ve.CV[rep, fold] <- fm$sigma[2]
    	h2t.CV[rep, fold] <- fm$sigma[1]/(fm$sigma[1]+fm$sigma[2])
    	e2.CV[rep, fold] <- fm$sigma[2]/(fm$sigma[1]+fm$sigma[2])
    
    	##Compute correlation and R2 in the test set
    	COR.CV[rep, fold] <- cor(yF[whichNa], yhatNa)
    	R2.CV[rep, fold] <- (cor(yF[whichNa], yhatNa))^2
    } else {
    	##Variance partition in the training set
    	Vt.CV[rep, fold] <- NA
    	Ve.CV[rep, fold] <- NA
    	h2t.CV[rep, fold] <- NA
    	e2.CV[rep, fold] <- NA
    
    	##Compute correlation and R2 in the test set
    	COR.CV[rep, fold] <- NA
    	R2.CV[rep, fold] <- NA  
    }
    
    print(paste("finished fold =", fold))
  }
  
  print(paste("finished CV =", rep))
}

###Compute mean and se across CV replicates
meanVt_TRM <- mean(Vt.CV, na.rm=T)
se_meanVt_TRM <- sd(Vt.CV, na.rm=T)/sqrt(sum(!is.na(Vt.CV)))
meanVe_TRM <- mean(Ve.CV, na.rm=T)
se_meanVe_TRM <- sd(Ve.CV, na.rm=T)/sqrt(sum(!is.na(Ve.CV)))
meanh2t_TRM <- mean(h2t.CV, na.rm=T)
se_meanh2t_TRM <- sd(h2t.CV, na.rm=T)/sqrt(sum(!is.na(h2t.CV)))
meane2_TRM <- mean(e2.CV, na.rm=T)
se_meane2_TRM <- sd(e2.CV, na.rm=T)/sqrt(sum(!is.na(e2.CV)))
meanCOR_TRM <- mean(COR.CV, na.rm=T)
se_meanCOR_TRM <- sd(COR.CV, na.rm=T)/sqrt(sum(!is.na(COR.CV)))
meanR2_TRM <- mean(R2.CV, na.rm=T)
se_meanR2_TRM <- sd(R2.CV, na.rm=T)/sqrt(sum(!is.na(R2.CV)))
nreps_no_NA_TRM <- sum(!is.na(R2.CV))



rm(list=setdiff(ls(), c("trait", "meanVt_TRM", "se_meanVt_TRM", "meanVe_TRM", "se_meanVe_TRM", "meanCOR_TRM", "se_meanCOR_TRM",
						"meanR2_TRM", "se_meanR2_TRM", "meanh2t_TRM", "se_meanh2t_TRM", "meane2_TRM", "se_meane2_TRM", "nreps_no_NA_TRM")))




######################
##Pearson's r matrix##
######################

###Source custom functions
source(file="../Software/myRfunctions.R")
load("Data/adjustData.RData")

###Load the data
#RRM
RRM <- as.matrix(read.table("Matrices/exp_cor_matrix_common_females.txt", header=TRUE, sep=" ", row.names=1))

#Phenotypes
pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)


###Select phenotype of interest (discard missing values)
phenoSel <- na.omit(pheno[, c(1, trait)])

###Find lines with both phenotypes and gene expression
common <- intersect(phenoSel[, 1], row.names(RRM))

###Select lines with expression in the phenotype data
phenoSel <- phenoSel[which(phenoSel[, 1] %in% common), ] 

###Adjust phenotype for inversion and wolbachia
phenoSel_adj <- adjustPheno(phenoSel, colnames(phenoSel)[2])
yF <- phenoSel_adj[!is.na(phenoSel_adj)]

###Get RRM for only lines with phenotype
RRM <- RRM[names(yF), names(yF)]

###Store matrices in a list needed in the prediction function
Mat <- as.list(NULL)
Mat[[1]] <- RRM

###Prediction (k-fold CV)
folds <- 5
nReps <- 30

COR.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(COR.CV) <- paste('fold=', 1:folds,sep='')

R2.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(R2.CV) <- paste('fold=', 1:folds,sep='')

h2t.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(h2t.CV) <- paste('fold=', 1:folds,sep='')

Vt.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(Vt.CV) <- paste('fold=', 1:folds,sep='')

Ve.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(Ve.CV) <- paste('fold=', 1:folds,sep='')

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
    
    ##TBLUP
    fm <- tryCatch(regress(yNa ~ 1, ~RRM, pos=c(T, T), maxcyc=10000), error=function(e){message(paste0("Error in CV = ", rep, ", fold = ", fold, ":\n", e)); return(NULL)}, 
                   warning=function(w){message(paste0("Warning in CV = ", rep, ", fold = ", fold, ":\n", w)); return(NULL)})
                   
    if(!is.null(fm) && fm$converged){   
    	##Phenotype prediction
    	yhat <- pred_y(fm, Mat, whichNa)
    	yhatNa <- yhat$VDtot  
    
    	##Variance partition in the training set
    	Vt.CV[rep, fold] <- fm$sigma[1]
    	Ve.CV[rep, fold] <- fm$sigma[2]
    	h2t.CV[rep, fold] <- fm$sigma[1]/(fm$sigma[1]+fm$sigma[2])
    	e2.CV[rep, fold] <- fm$sigma[2]/(fm$sigma[1]+fm$sigma[2])
    
    	##Compute correlation and R2 in the test set
    	COR.CV[rep, fold] <- cor(yF[whichNa], yhatNa)
    	R2.CV[rep, fold] <- (cor(yF[whichNa], yhatNa))^2
    } else {
    	##Variance partition in the training set
    	Vt.CV[rep, fold] <- NA
    	Ve.CV[rep, fold] <- NA
    	h2t.CV[rep, fold] <- NA
    	e2.CV[rep, fold] <- NA
    
    	##Compute correlation and R2 in the test set
    	COR.CV[rep, fold] <- NA
    	R2.CV[rep, fold] <- NA  
    }
    
    print(paste("finished fold =", fold))
  }
  
  print(paste("finished CV =", rep))
}

###Compute mean and se across CV replicates
meanVt_RRM <- mean(Vt.CV, na.rm=T)
se_meanVt_RRM <- sd(Vt.CV, na.rm=T)/sqrt(sum(!is.na(Vt.CV)))
meanVe_RRM <- mean(Ve.CV, na.rm=T)
se_meanVe_RRM <- sd(Ve.CV, na.rm=T)/sqrt(sum(!is.na(Ve.CV)))
meanh2t_RRM <- mean(h2t.CV, na.rm=T)
se_meanh2t_RRM <- sd(h2t.CV, na.rm=T)/sqrt(sum(!is.na(h2t.CV)))
meane2_RRM <- mean(e2.CV, na.rm=T)
se_meane2_RRM <- sd(e2.CV, na.rm=T)/sqrt(sum(!is.na(e2.CV)))
meanCOR_RRM <- mean(COR.CV, na.rm=T)
se_meanCOR_RRM <- sd(COR.CV, na.rm=T)/sqrt(sum(!is.na(COR.CV)))
meanR2_RRM <- mean(R2.CV, na.rm=T)
se_meanR2_RRM <- sd(R2.CV, na.rm=T)/sqrt(sum(!is.na(R2.CV)))
nreps_no_NA_RRM <- sum(!is.na(R2.CV))



rm(list=setdiff(ls(), c("trait", "meanVt_TRM", "se_meanVt_TRM", "meanVe_TRM", "se_meanVe_TRM", "meanCOR_TRM", "se_meanCOR_TRM",
						"meanR2_TRM", "se_meanR2_TRM", "meanh2t_TRM", "se_meanh2t_TRM", "meane2_TRM", "se_meane2_TRM", "nreps_no_NA_TRM", "meanVt_RRM", "se_meanVt_RRM", "meanVe_RRM", "se_meanVe_RRM", 
						"meanCOR_RRM", "se_meanCOR_RRM", "meanR2_RRM", "se_meanR2_RRM", "meanh2t_RRM", "se_meanh2t_RRM", "meane2_RRM", "se_meane2_RRM", "nreps_no_NA_RRM")))
						
						
						
						
##############################
##RKHS with one kernel h=0.5##
##############################

###Source custom functions
source(file="../Software/myRfunctions.R")
load("Data/adjustData.RData")

###Load the data
#KRM
KRM <- as.matrix(read.table("Matrices/K_matrix_h_0.5_common_females.txt", header=TRUE, sep=" ", row.names=1))

#Phenotypes
pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)


###Select phenotype of interest (discard missing values)
phenoSel <- na.omit(pheno[, c(1, trait)])

###Find lines with both phenotypes and gene expression
common <- intersect(phenoSel[, 1], row.names(KRM))

###Select lines with expression in the phenotype data
phenoSel <- phenoSel[which(phenoSel[, 1] %in% common), ] 

###Adjust phenotype for inversion and wolbachia
phenoSel_adj <- adjustPheno(phenoSel, colnames(phenoSel)[2])
yF <- phenoSel_adj[!is.na(phenoSel_adj)]

###Get KRM for only lines with phenotype
KRM <- KRM[names(yF), names(yF)]

###Store matrices in a list needed in the prediction function
Mat <- as.list(NULL)
Mat[[1]] <- KRM


###Prediction (k-fold CV)
folds <- 5
nReps <- 30

COR.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(COR.CV) <- paste('fold=', 1:folds,sep='')

R2.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(R2.CV) <- paste('fold=', 1:folds,sep='')

h2t.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(h2t.CV) <- paste('fold=', 1:folds,sep='')

Vt.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(Vt.CV) <- paste('fold=', 1:folds,sep='')

Ve.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(Ve.CV) <- paste('fold=', 1:folds,sep='')

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
    
    ##TBLUP
    fm <- tryCatch(regress(yNa ~ 1, ~KRM, pos=c(T, T), maxcyc=10000), error=function(e){message(paste0("Error in CV = ", rep, ", fold = ", fold, ":\n", e)); return(NULL)}, 
                   warning=function(w){message(paste0("Warning in CV = ", rep, ", fold = ", fold, ":\n", w)); return(NULL)})
                   
    if(!is.null(fm) && fm$converged){   
    	##Phenotype prediction
    	yhat <- pred_y(fm, Mat, whichNa)
    	yhatNa <- yhat$VDtot  
    
    	##Variance partition in the training set
    	Vt.CV[rep, fold] <- fm$sigma[1]
    	Ve.CV[rep, fold] <- fm$sigma[2]
    	h2t.CV[rep, fold] <- fm$sigma[1]/(fm$sigma[1]+fm$sigma[2])
    	e2.CV[rep, fold] <- fm$sigma[2]/(fm$sigma[1]+fm$sigma[2])
    
    	##Compute correlation and R2 in the test set
    	COR.CV[rep, fold] <- cor(yF[whichNa], yhatNa)
    	R2.CV[rep, fold] <- (cor(yF[whichNa], yhatNa))^2
    } else {
    	##Variance partition in the training set
    	Vt.CV[rep, fold] <- NA
    	Ve.CV[rep, fold] <- NA
    	h2t.CV[rep, fold] <- NA
    	e2.CV[rep, fold] <- NA
    
    	##Compute correlation and R2 in the test set
    	COR.CV[rep, fold] <- NA
    	R2.CV[rep, fold] <- NA  
    }
    
    print(paste("finished fold =", fold))
  }

  print(paste("finished CV =", rep))
}

###Compute mean and se across CV replicates
meanVt_KRM_0.5 <- mean(Vt.CV, na.rm=T)
se_meanVt_KRM_0.5 <- sd(Vt.CV, na.rm=T)/sqrt(sum(!is.na(Vt.CV)))
meanVe_KRM_0.5 <- mean(Ve.CV, na.rm=T)
se_meanVe_KRM_0.5 <- sd(Ve.CV, na.rm=T)/sqrt(sum(!is.na(Ve.CV)))
meanh2t_KRM_0.5 <- mean(h2t.CV, na.rm=T)
se_meanh2t_KRM_0.5 <- sd(h2t.CV, na.rm=T)/sqrt(sum(!is.na(h2t.CV)))
meane2_KRM_0.5 <- mean(e2.CV, na.rm=T)
se_meane2_KRM_0.5 <- sd(e2.CV, na.rm=T)/sqrt(sum(!is.na(e2.CV)))
meanCOR_KRM_0.5 <- mean(COR.CV, na.rm=T)
se_meanCOR_KRM_0.5 <- sd(COR.CV, na.rm=T)/sqrt(sum(!is.na(COR.CV)))
meanR2_KRM_0.5 <- mean(R2.CV, na.rm=T)
se_meanR2_KRM_0.5 <- sd(R2.CV, na.rm=T)/sqrt(sum(!is.na(R2.CV)))
nreps_no_NA_KRM_0.5 <- sum(!is.na(R2.CV))


rm(list=setdiff(ls(), c("trait", "meanVt_TRM", "se_meanVt_TRM", "meanVe_TRM", "se_meanVe_TRM", "meanCOR_TRM", "se_meanCOR_TRM",
						"meanR2_TRM", "se_meanR2_TRM", "meanh2t_TRM", "se_meanh2t_TRM", "meane2_TRM", "se_meane2_TRM", "nreps_no_NA_TRM", "meanVt_RRM", "se_meanVt_RRM", "meanVe_RRM", "se_meanVe_RRM", 
						"meanCOR_RRM", "se_meanCOR_RRM", "meanR2_RRM", "se_meanR2_RRM", "meanh2t_RRM", "se_meanh2t_RRM", "meane2_RRM", "se_meane2_RRM", "nreps_no_NA_RRM", "meanVt_KRM_0.5", "se_meanVt_KRM_0.5", 
						"meanVe_KRM_0.5", "se_meanVe_KRM_0.5", "meanCOR_KRM_0.5", "se_meanCOR_KRM_0.5", "meanR2_KRM_0.5", "se_meanR2_KRM_0.5", "meanh2t_KRM_0.5", "se_meanh2t_KRM_0.5", "meane2_KRM_0.5", "se_meane2_KRM_0.5", "nreps_no_NA_KRM_0.5")))
	
	
	
	
						
##############################
##RKHS with one kernel h=0.1##
##############################

###Source custom functions
source(file="../Software/myRfunctions.R")
load("Data/adjustData.RData")

###Load the data
#KRM
KRM <- as.matrix(read.table("Matrices/K_matrix_h_0.1_common_females.txt", header=TRUE, sep=" ", row.names=1))

#Phenotypes
pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)


###Select phenotype of interest (discard missing values)
phenoSel <- na.omit(pheno[, c(1, trait)])

###Find lines with both phenotypes and gene expression
common <- intersect(phenoSel[, 1], row.names(KRM))

###Select lines with expression in the phenotype data
phenoSel <- phenoSel[which(phenoSel[, 1] %in% common), ] 

###Adjust phenotype for inversion and wolbachia
phenoSel_adj <- adjustPheno(phenoSel, colnames(phenoSel)[2])
yF <- phenoSel_adj[!is.na(phenoSel_adj)]

###Get KRM for only lines with phenotype
KRM <- KRM[names(yF), names(yF)]

###Store matrices in a list needed in the prediction function
Mat <- as.list(NULL)
Mat[[1]] <- KRM


###Prediction (k-fold CV)
folds <- 5
nReps <- 30

COR.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(COR.CV) <- paste('fold=', 1:folds,sep='')

R2.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(R2.CV) <- paste('fold=', 1:folds,sep='')

h2t.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(h2t.CV) <- paste('fold=', 1:folds,sep='')

Vt.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(Vt.CV) <- paste('fold=', 1:folds,sep='')

Ve.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(Ve.CV) <- paste('fold=', 1:folds,sep='')

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
    
    ##TBLUP
    fm <- tryCatch(regress(yNa ~ 1, ~KRM, pos=c(T, T), maxcyc=10000), error=function(e){message(paste0("Error in CV = ", rep, ", fold = ", fold, ":\n", e)); return(NULL)}, 
                   warning=function(w){message(paste0("Warning in CV = ", rep, ", fold = ", fold, ":\n", w)); return(NULL)})
                   
    if(!is.null(fm) && fm$converged){   
    	##Phenotype prediction
    	yhat <- pred_y(fm, Mat, whichNa)
    	yhatNa <- yhat$VDtot  
    
    	##Variance partition in the training set
    	Vt.CV[rep, fold] <- fm$sigma[1]
    	Ve.CV[rep, fold] <- fm$sigma[2]
    	h2t.CV[rep, fold] <- fm$sigma[1]/(fm$sigma[1]+fm$sigma[2])
    	e2.CV[rep, fold] <- fm$sigma[2]/(fm$sigma[1]+fm$sigma[2])
    
    	##Compute correlation and R2 in the test set
    	COR.CV[rep, fold] <- cor(yF[whichNa], yhatNa)
    	R2.CV[rep, fold] <- (cor(yF[whichNa], yhatNa))^2
    } else {
    	##Variance partition in the training set
    	Vt.CV[rep, fold] <- NA
    	Ve.CV[rep, fold] <- NA
    	h2t.CV[rep, fold] <- NA
    	e2.CV[rep, fold] <- NA
    
    	##Compute correlation and R2 in the test set
    	COR.CV[rep, fold] <- NA
    	R2.CV[rep, fold] <- NA  
    }
    
    print(paste("finished fold =", fold))
  }

  print(paste("finished CV =", rep))
}

###Compute mean and se across CV replicates
meanVt_KRM_0.1 <- mean(Vt.CV, na.rm=T)
se_meanVt_KRM_0.1 <- sd(Vt.CV, na.rm=T)/sqrt(sum(!is.na(Vt.CV)))
meanVe_KRM_0.1 <- mean(Ve.CV, na.rm=T)
se_meanVe_KRM_0.1 <- sd(Ve.CV, na.rm=T)/sqrt(sum(!is.na(Ve.CV)))
meanh2t_KRM_0.1 <- mean(h2t.CV, na.rm=T)
se_meanh2t_KRM_0.1 <- sd(h2t.CV, na.rm=T)/sqrt(sum(!is.na(h2t.CV)))
meane2_KRM_0.1 <- mean(e2.CV, na.rm=T)
se_meane2_KRM_0.1 <- sd(e2.CV, na.rm=T)/sqrt(sum(!is.na(e2.CV)))
meanCOR_KRM_0.1 <- mean(COR.CV, na.rm=T)
se_meanCOR_KRM_0.1 <- sd(COR.CV, na.rm=T)/sqrt(sum(!is.na(COR.CV)))
meanR2_KRM_0.1 <- mean(R2.CV, na.rm=T)
se_meanR2_KRM_0.1 <- sd(R2.CV, na.rm=T)/sqrt(sum(!is.na(R2.CV)))
nreps_no_NA_KRM_0.1 <- sum(!is.na(R2.CV))



rm(list=setdiff(ls(), c("pheno", "trait", "meanVt_TRM", "se_meanVt_TRM", "meanVe_TRM", "se_meanVe_TRM", "meanCOR_TRM", "se_meanCOR_TRM",
						"meanR2_TRM", "se_meanR2_TRM", "meanh2t_TRM", "se_meanh2t_TRM", "meane2_TRM", "se_meane2_TRM", "nreps_no_NA_TRM", "meanVt_RRM", "se_meanVt_RRM", "meanVe_RRM", "se_meanVe_RRM", 
						"meanCOR_RRM", "se_meanCOR_RRM", "meanR2_RRM", "se_meanR2_RRM", "meanh2t_RRM", "se_meanh2t_RRM", "meane2_RRM", "se_meane2_RRM", "nreps_no_NA_RRM", "meanVt_KRM_0.5", "se_meanVt_KRM_0.5", 
						"meanVe_KRM_0.5", "se_meanVe_KRM_0.5", "meanCOR_KRM_0.5", "se_meanCOR_KRM_0.5", "meanR2_KRM_0.5", "se_meanR2_KRM_0.5", "meanh2t_KRM_0.5", "se_meanh2t_KRM_0.5", "meane2_KRM_0.5", "se_meane2_KRM_0.5", "nreps_no_NA_KRM_0.5",
						"meanVt_KRM_0.1", "se_meanVt_KRM_0.1", "meanVe_KRM_0.1", "se_meanVe_KRM_0.1", "meanCOR_KRM_0.1", "se_meanCOR_KRM_0.1", 
						"meanR2_KRM_0.1", "se_meanR2_KRM_0.1", "meanh2t_KRM_0.1", "se_meanh2t_KRM_0.1", "meane2_KRM_0.1", "se_meane2_KRM_0.1", "nreps_no_NA_KRM_0.1")))								


###Put the results together
Results <- matrix(c(meanVt_TRM, se_meanVt_TRM, meanVe_TRM, se_meanVe_TRM, meanCOR_TRM, se_meanCOR_TRM, meanR2_TRM, se_meanR2_TRM, meanh2t_TRM, se_meanh2t_TRM, meane2_TRM, se_meane2_TRM, nreps_no_NA_TRM,
                    meanVt_RRM, se_meanVt_RRM, meanVe_RRM, se_meanVe_RRM, meanCOR_RRM, se_meanCOR_RRM, meanR2_RRM, se_meanR2_RRM, meanh2t_RRM, se_meanh2t_RRM, meane2_RRM, se_meane2_RRM, nreps_no_NA_RRM,
                    meanVt_KRM_0.5, se_meanVt_KRM_0.5, meanVe_KRM_0.5, se_meanVe_KRM_0.5, meanCOR_KRM_0.5, se_meanCOR_KRM_0.5, meanR2_KRM_0.5, se_meanR2_KRM_0.5, meanh2t_KRM_0.5, se_meanh2t_KRM_0.5, meane2_KRM_0.5, se_meane2_KRM_0.5, nreps_no_NA_KRM_0.5,
                    meanVt_KRM_0.1, se_meanVt_KRM_0.1, meanVe_KRM_0.1, se_meanVe_KRM_0.1, meanCOR_KRM_0.1, se_meanCOR_KRM_0.1, meanR2_KRM_0.1, se_meanR2_KRM_0.1, meanh2t_KRM_0.1, se_meanh2t_KRM_0.1, meane2_KRM_0.1, se_meane2_KRM_0.1, nreps_no_NA_KRM_0.1), ncol=13, byrow=T)		
colnames(Results) <- c("mean Vt", "se Vt", "mean Ve", "se Ve", "mean r", "se r", "mean r2", "se r2", "mean h2t", "se h2t", "mean e2", "se e2", "nreps")
Method <- c("TRM", "rRM", "KRM_h=0.5", "KRM_h=0.1")
Results_final <- data.frame(Method, Results)

###Write to a file
write.table(Results_final, paste('Results/Pred_Expression_regress_adjPheno_common_', colnames(pheno)[trait],'_females.csv', sep=""), 
		  col.names=T, row.names=F, quote=F, sep=",")




						
