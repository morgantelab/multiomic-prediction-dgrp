########################
## Fabio Morgante
## 3/19/2017
## BLUP with genotype + expression data
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

##Load libraries
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

#GRM
GRM <- as.matrix(read.table("Matrices/GRM_200lines_common.txt", header=TRUE, sep="\t", row.names=1))

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

###Get GRM for only lines with phenotype
GRM <- GRM[names(yF), names(yF)]

###Get TRM for only lines with phenotype
TRM <- TRM[names(yF), names(yF)]

###Store matrices in a list needed in the prediction function
Mat <- as.list(NULL)
Mat[[1]] <- GRM
Mat[[2]] <- TRM


###Prediction (k-fold CV)
folds <- 5
nReps <- 30

COR.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(COR.CV) <- paste('fold=', 1:folds,sep='')

R2.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(R2.CV) <- paste('fold=', 1:folds,sep='')

Vg.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(Vg.CV) <- paste('fold=', 1:folds,sep='')

Vt.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(Vt.CV) <- paste('fold=', 1:folds,sep='')

Ve.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(Ve.CV) <- paste('fold=', 1:folds,sep='')

h2g.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(h2g.CV) <- paste('fold=', 1:folds,sep='')

h2t.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(h2t.CV) <- paste('fold=', 1:folds,sep='')

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
    
    ##GTBLUP
    fm <- tryCatch(regress(yNa ~ 1, ~GRM+TRM, pos=c(T, T, T), maxcyc=10000), error=function(e){message(paste0("Error in CV = ", rep, ", fold = ", fold, ":\n", e)); return(NULL)}, 
                   warning=function(w){message(paste0("Warning in CV = ", rep, ", fold = ", fold, ":\n", w)); return(NULL)})

    if(!is.null(fm) && fm$converged){    
    	##Phenotype prediction
    	yhat <- pred_y(fm, Mat, whichNa)
    	yhatNa <- yhat$VDtot 
    
    	##Variance partition in the training set
    	Vg.CV[rep, fold] <- fm$sigma[1]
    	Vt.CV[rep, fold] <- fm$sigma[2]
    	Ve.CV[rep, fold] <- fm$sigma[3]
    	h2g.CV[rep, fold] <- fm$sigma[1]/(fm$sigma[1]+fm$sigma[2]+fm$sigma[3])
    	h2t.CV[rep, fold] <- fm$sigma[2]/(fm$sigma[1]+fm$sigma[2]+fm$sigma[3])
    	e2.CV[rep, fold] <- fm$sigma[3]/(fm$sigma[1]+fm$sigma[2]+fm$sigma[3])
    
    	##Compute correlation and R2 in the test set
    	COR.CV[rep, fold] <- cor(yF[whichNa], yhatNa)
    	R2.CV[rep, fold] <- (cor(yF[whichNa], yhatNa))^2
    } else {
     	##Variance partition in the training set
    	Vg.CV[rep, fold] <- NA
    	Vt.CV[rep, fold] <- NA
    	Ve.CV[rep, fold] <- NA
    	h2g.CV[rep, fold] <- NA
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
meanVg <- mean(Vg.CV, na.rm=T)
se_meanVg <- sd(Vg.CV, na.rm=T)/sqrt(sum(!is.na(Vg.CV)))
meanVt <- mean(Vt.CV, na.rm=T)
se_meanVt <- sd(Vt.CV, na.rm=T)/sqrt(sum(!is.na(Vt.CV)))
meanVe <- mean(Ve.CV, na.rm=T)
se_meanVe <- sd(Ve.CV, na.rm=T)/sqrt(sum(!is.na(Ve.CV)))
meanh2g <- mean(h2g.CV, na.rm=T)
se_meanh2g <- sd(h2g.CV, na.rm=T)/sqrt(sum(!is.na(h2g.CV)))
meanh2t <- mean(h2t.CV, na.rm=T)
se_meanh2t <- sd(h2t.CV, na.rm=T)/sqrt(sum(!is.na(h2t.CV)))
meane2 <- mean(e2.CV, na.rm=T)
se_meane2 <- sd(e2.CV, na.rm=T)/sqrt(sum(!is.na(e2.CV)))
meanCOR <- mean(COR.CV, na.rm=T)
se_meanCOR <- sd(COR.CV, na.rm=T)/sqrt(sum(!is.na(COR.CV)))
meanR2 <- mean(R2.CV, na.rm=T)
se_meanR2 <- sd(R2.CV, na.rm=T)/sqrt(sum(!is.na(R2.CV)))
nreps_no_NA <- sum(!is.na(R2.CV))



###Put the results together
Results <- matrix(c(meanVg, se_meanVg, meanVt, se_meanVt, meanVe, se_meanVe, meanCOR, se_meanCOR, meanR2, se_meanR2, meanh2g, se_meanh2g, meanh2t, se_meanh2t, meane2, se_meane2, nreps_no_NA), ncol=17, byrow=T)		
colnames(Results) <- c("mean Vg", "se Vg", "mean Vt", "se Vt", "mean Ve", "se Ve", "mean r", "se r", "mean r2", "se r2", "mean h2g", "se h2g", "mean h2t", "se h2t", "mean e2", "se e2", "nreps")
Method <- c("GRM+TRM")
Results_final <- data.frame(Method, Results)

###Write to a file
write.table(Results_final, paste('Results/Pred_SNPs+Expression_regress_adjPheno_common_', colnames(pheno)[trait],'_females.csv', sep=""), 
		  col.names=T, row.names=F, quote=F, sep=",")




						
