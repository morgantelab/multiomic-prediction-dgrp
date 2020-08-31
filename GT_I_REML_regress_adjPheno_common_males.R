########################
## Fabio Morgante
## 3/23/2017
## Variance partition
##
## Modified on 4-5-2018
## Reason: adjust phenotypes, use updated GRM  
## 
## Modified on 11/11/2018
## Reason: add options(stringsAsFactors=FALSE) and remove as.is=T        
##
## Modified on 05/06/2020
## Reason: Use relative paths, use as.numeric(NA),
## 		   set maxcyc=10000, add options(warn=1)             
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
####Males####
###############

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

###Get TRM for only lines with phenotype
TRM <- TRM[names(yM), names(yM)]

###Compute GRM#TRM
IRM <- GRM*TRM


###GREML
fm_G <- regress(yM ~ 1, ~GRM, pos=c(T, T), maxcyc=10000)
Vg_G <- fm_G$sigma[1]
Vt_G <- as.numeric(NA)
Vi_G <- as.numeric(NA)
Ve_G <- fm_G$sigma[2]


###TREML
fm_T <- regress(yM ~ 1, ~TRM, pos=c(T, T), maxcyc=10000)
Vg_T <- as.numeric(NA)
Vt_T <- fm_T$sigma[1]
Vi_T <- as.numeric(NA)
Ve_T <- fm_T$sigma[2]

###GTREML
fm_GT <- regress(yM ~ 1, ~GRM+TRM, pos=c(T, T, T), maxcyc=10000)
Vg_GT <- fm_GT$sigma[1]
Vt_GT <- fm_GT$sigma[2]
Vi_GT <- as.numeric(NA)
Ve_GT <- fm_GT$sigma[3]

###GT-I-REML
fm_GTI <- regress(yM ~ 1, ~GRM+TRM+IRM, pos=c(T, T, T, T), maxcyc=10000)
Vg_GTI <- fm_GTI$sigma[1]
Vt_GTI <- fm_GTI$sigma[2]
Vi_GTI <- fm_GTI$sigma[3]
Ve_GTI <- fm_GTI$sigma[4]


Results <- matrix(c(Vg_G, Vt_G, Vi_G, Ve_G,
					Vg_T, Vt_T, Vi_T, Ve_T,
					Vg_GT, Vt_GT, Vi_GT, Ve_GT,
					Vg_GTI, Vt_GTI, Vi_GTI, Ve_GTI), ncol=4, byrow=T)
colnames(Results) <- c("Vg", "Vt", "Vi", "Ve")
Method <- c("GRM", "TRM", "GRM+TRM", "GRM+TRM+IRM")

Results_final <- data.frame(Method, Results)


###Write to a file
write.table(Results_final, paste('Results/Variance_Partition_SNPs_Expression_Interaction_regress_adjPheno_common_', colnames(pheno)[trait],'_males.csv', sep=""), 
		  col.names=T, row.names=F, quote=F, sep=",")




						
