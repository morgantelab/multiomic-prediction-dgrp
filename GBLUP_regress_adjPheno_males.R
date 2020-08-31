########################
## Fabio Morgante
## 2/6/2017
## BLUP with GRM
##
## Modified on 4-5-2018
## Reason: adjust phenotypes, use updated GRM
##
## Modified on 4-6-2018
## Reason: calculate h2g and pool folds and reps to calculate se
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
####males####
###############

#######
##GRM##
#######

###Source custom functions
source(file="../Software/myRfunctions.R")
load("Data/adjustData.RData")

###Load the data
#TRM
TRM <- as.matrix(read.table("Matrices/TRM_common_males.txt", header=TRUE, sep=" ", row.names=1))

#GRM
GRM <- as.matrix(read.table("Matrices/GRM_200lines_common.txt", header=TRUE, sep="\t", row.names=1))

#Phenotypes
pheno <- read.csv("Data/eQTL_traits_males.csv", header=TRUE)


###Select phenotype of interest (discard missing values)
phenoSel <- na.omit(pheno[, c(1, trait)])

###Find lines with both phenotypes and gene expression
common <- intersect(phenoSel[, 1], row.names(TRM))

###Select lines with expression in the phenotype data
phenoSel <- phenoSel[which(phenoSel[, 1] %in% common), ] 

###Adjust phenotype for inversion and wolbachia
phenoSel_adj <- adjustPheno(phenoSel, colnames(phenoSel)[2])
yM <- phenoSel_adj[!is.na(phenoSel_adj)]

###Get GRM for only lines with phenotype
GRM <- GRM[names(yM), names(yM)]

###Store matrices in a list needed in the prediction function
Mat <- as.list(NULL)
Mat[[1]] <- GRM


###Prediction (k-fold CV)
folds <- 5
nReps <- 30

COR.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(COR.CV) <- paste('fold=', 1:folds,sep='')

R2.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(R2.CV) <- paste('fold=', 1:folds,sep='')

Vg.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(Vg.CV) <- paste('fold=', 1:folds,sep='')

Ve.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(Ve.CV) <- paste('fold=', 1:folds,sep='')

h2g.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(h2g.CV) <- paste('fold=', 1:folds,sep='')

e2.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(e2.CV) <- paste('fold=', 1:folds,sep='')



for(rep in 1:nReps){
  ##Create sets for dividing the lines into folds (if statement needed for when length(yM)/folds) is NOT an integer)
  set.seed(120+rep)
  sets1 <- rep(1:folds, (length(yM)/folds))
  if(length(sets1) != length(yM)){sets2 <- sample(x=1:folds, size=(length(yM)-length(sets1)))
  								  sets <- c(sets1, sets2)} else {sets <- sets1}
  sets <- sets[order(runif(length(yM)))]
    
  for(fold in 1:folds){
    ##Assign NAs to lines in the test set
    yNa <- yM
    whichNa <- which(sets==fold)
    yNa[whichNa] <- NA
    
    ##GBLUP
    fm <- tryCatch(regress(yNa ~ 1, ~GRM, pos=c(T, T), maxcyc=10000), error=function(e){message(paste0("Error in CV = ", rep, ", fold = ", fold, ":\n", e)); return(NULL)}, 
                   warning=function(w){message(paste0("Warning in CV = ", rep, ", fold = ", fold, ":\n", w)); return(NULL)})
    
    if(!is.null(fm) && fm$converged){   
    	##Phenotype prediction
    	yhat <- pred_y(fm, Mat, whichNa)
    	yhatNa <- yhat$VDtot  
    
    	##Variance partition in the training set
    	Vg.CV[rep, fold] <- fm$sigma[1]
    	Ve.CV[rep, fold] <- fm$sigma[2]
    	h2g.CV[rep, fold] <- fm$sigma[1]/(fm$sigma[1]+fm$sigma[2])
    	e2.CV[rep, fold] <- fm$sigma[2]/(fm$sigma[1]+fm$sigma[2])
    
    	##Compute correlation and R2 in the test set
    	COR.CV[rep, fold] <- cor(yM[whichNa], yhatNa)
    	R2.CV[rep, fold] <- (cor(yM[whichNa], yhatNa))^2
    } else {
    	##Variance partition in the training set
    	Vg.CV[rep, fold] <- NA
    	Ve.CV[rep, fold] <- NA
    	h2g.CV[rep, fold] <- NA
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
meanVg_GRM <- mean(Vg.CV, na.rm=T)
se_meanVg_GRM <- sd(Vg.CV, na.rm=T)/sqrt(sum(!is.na(Vg.CV)))
meanVe_GRM <- mean(Ve.CV, na.rm=T)
se_meanVe_GRM <- sd(Ve.CV, na.rm=T)/sqrt(sum(!is.na(Ve.CV)))
meanh2g_GRM <- mean(h2g.CV, na.rm=T)
se_meanh2g_GRM <- sd(h2g.CV, na.rm=T)/sqrt(sum(!is.na(h2g.CV)))
meane2_GRM <- mean(e2.CV, na.rm=T)
se_meane2_GRM <- sd(e2.CV, na.rm=T)/sqrt(sum(!is.na(e2.CV)))
meanCOR_GRM <- mean(COR.CV, na.rm=T)
se_meanCOR_GRM <- sd(COR.CV, na.rm=T)/sqrt(sum(!is.na(COR.CV)))
meanR2_GRM <- mean(R2.CV, na.rm=T)
se_meanR2_GRM <- sd(R2.CV, na.rm=T)/sqrt(sum(!is.na(R2.CV)))
nreps_no_NA_GRM <- sum(!is.na(R2.CV))




rm(list=setdiff(ls(), c("pheno", "trait", "meanVg_GRM", "se_meanVg_GRM", "meanVe_GRM", "se_meanVe_GRM", "meanCOR_GRM", "se_meanCOR_GRM",
						"meanR2_GRM", "se_meanR2_GRM", "meanh2g_GRM", "se_meanh2g_GRM", "meane2_GRM", "se_meane2_GRM", "nreps_no_NA_GRM")))



###Put the results together
Results <- matrix(c(meanVg_GRM, se_meanVg_GRM, meanVe_GRM, se_meanVe_GRM, meanCOR_GRM, se_meanCOR_GRM, meanR2_GRM, se_meanR2_GRM, meanh2g_GRM, se_meanh2g_GRM, meane2_GRM, se_meane2_GRM, nreps_no_NA_GRM), ncol=13, byrow=T)		
colnames(Results) <- c("mean Vg", "se Vg", "mean Ve", "se Ve", "mean r", "se r", "mean r2", "se r2", "mean h2g", "se h2g", "mean e2", "se e2", "nreps")
Method <- c("GRM")
Results_final <- data.frame(Method, Results)

###Write to a file
write.table(Results_final, paste('Results/Pred_SNPs_regress_adjPheno_', colnames(pheno)[trait],'_males.csv', sep=""), 
		  col.names=T, row.names=F, quote=F, sep=",")
