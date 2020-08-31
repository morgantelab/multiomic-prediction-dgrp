########################
## Fabio Morgante
## 7/10/2017
## Random Forest with expression data
##
## Modified on 4-5-2018
## Reason: use expression data adjusted ONLY for alignment bias after Logan fixed the issue
##			and adjust pheno
##
## Modified on 4-6-2018
## Reason: calculate h2 and pool folds and reps to calculate se
## 
## Modified on 11/11/2018
## Reason: add options(stringsAsFactors=FALSE) and remove as.is=T        
##
## Modified on 05/06/2020
## Reason: Use relative paths, use as.numeric(NA),
## 		   add options(warn=1)          
##
## Modified on 05/29/2020
## Reason: Use dorng to ensure reproducibility
#############################

###Set global options
options(stringsAsFactors=FALSE)
options(warn=1)

###Select trait to analyze from command line arguments
args <- commandArgs(TRUE)
trait <- as.numeric(args[1])
ntrees <- as.numeric(args[2])
cores <- as.numeric(args[3])

###Set number of cores for R open
setMKLthreads(1)

###Load libraries
library(foreach)
library(doRNG)
library(doMC)
library(randomForest)

###Register cores
registerDoMC(cores)


###############
####Females####
###############

###Source custom functions
load("Data/adjustData.RData")

###Load the data
#Transcript data
dataF <- as.matrix(read.table("Data/combined_samples_known_novel_fpkm_VR_NoWol_F_line_means_rename_noflag_transp_common.txt",
					header=T, row.names=1, sep="\t"))

#Phenotypes
pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)


###Select phenotype of interest (discard missing values)
phenoSel <- na.omit(pheno[, c(1, trait)])

###Find lines with both phenotypes and gene expression
common <- intersect(phenoSel[, 1], row.names(dataF))

###Select lines with expression in the phenotype data
phenoSel <- phenoSel[which(phenoSel[, 1] %in% common), ] 

###Adjust phenotype for inversion and wolbachia
phenoSel_adj <- adjustPheno(phenoSel, colnames(phenoSel)[2])
yF <- phenoSel_adj[!is.na(phenoSel_adj)]

###Get transcript data  for only lines with phenotype
dataFSel <- dataF[names(yF), ]


###Prediction (k-fold CV)
folds <- 5
nReps <- 30

COR.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(COR.CV) <- paste('fold=', 1:folds,sep='')

R2.CV <- matrix(as.numeric(NA), ncol=(folds), nrow=(nReps))
colnames(R2.CV) <- paste('fold=', 1:folds,sep='')


for(rep in 1:nReps){
  ##Create sets for dividing the lines into folds (if statement needed for when length(yF)/folds) is NOT an integer)
  set.seed(120+rep)
  sets1 <- rep(1:folds, (length(yF)/folds))
  if(length(sets1) != length(yF)){sets2 <- sample(x=1:folds, size=(length(yF)-length(sets1)))
  								  sets <- c(sets1, sets2)} else {sets <- sets1}
  sets <- sets[order(runif(length(yF)))]
    
  for(fold in 1:folds){
    ##Assign NAs to lines in the test set
    whichNa <- which(sets==fold)
    yTrain <- yF[-whichNa]
    xTrain <- dataFSel[-whichNa, ]
    yTest <- yF[whichNa]
    xTest <- dataFSel[whichNa, ]
    
    ##Random Forest model
    fit <- foreach(ntree=rep((ntrees/cores), cores), .combine=combine, .multicombine=T, .packages='randomForest') %dorng% {
    				randomForest(xTrain, yTrain, ntree=ntree, importance=F)}
    
    ##Phenotype prediction
    yhatNa <- predict(fit, xTest)  
        
    ##Compute correlation and R2 in the test set
    COR.CV[rep, fold] <- cor(yTest, yhatNa)
    R2.CV[rep, fold] <- (cor(yTest, yhatNa))^2
    
    print(paste("finished fold =", fold))
  }

  print(paste("finished CV =", rep))
}

###Compute mean and se across CV replicates
meanCOR <- mean(COR.CV, na.rm=T)
se_meanCOR <- sd(COR.CV, na.rm=T)/sqrt(sum(!is.na(COR.CV)))
meanR2 <- mean(R2.CV, na.rm=T)
se_meanR2 <- sd(R2.CV, na.rm=T)/sqrt(sum(!is.na(R2.CV)))



rm(list=setdiff(ls(), c("pheno", "trait", "ntrees", "meanCOR", "se_meanCOR", "meanR2", "se_meanR2")))



###Put the results together
Results <- matrix(c(meanCOR, se_meanCOR, meanR2, se_meanR2), ncol=4, byrow=T)		
colnames(Results) <- c("mean r", "se r", "mean r2", "se r2")
Method <- c("TRF")
nTrees <- c(paste(ntrees))
Results_final <- data.frame(Method, nTrees, Results)

###Write to a file
write.table(Results_final, paste('Results/Pred_Expression_RandomForest_adjPheno_common_', colnames(pheno)[trait],'_females.csv', sep=""), 
		  col.names=T, row.names=F, quote=F, sep=",")




						
