##################################
## Fabio Morgante		       
## 4-23-2018			       
## Compute centered and scaled 
## genotype matrix from common 
## variants	(200 lines; standard filters)	
##
## Modified on 11/11/2018
## Reason: add options(stringsAsFactors=FALSE) and remove as.is=T		   
##
## Modified on 05/06/2020
## Reason: Use relative paths             
#############################

###Load libraries
library(snpStats)

###Set number of cores
setMKLthreads(5)

###Set global options
options(stringsAsFactors=FALSE)

###Load the data
geno  <- as(read.plink("Data/dgrp2_200lines_common_prediction")$genotypes, "numeric")

###Calculate the mean of each marker
mean_mar <- apply(geno, 2, mean, na.rm=TRUE)

###Substitute NAs with the mean of each marker
for(i in 1:ncol(geno)){
  NAs <- which(is.na(geno[, i]))
  geno[NAs, i] <- mean_mar[i]
  if((i%%10000)==1){print(paste("SNP",i))}
}

###Center and standardize the genotype matrix
W <- scale(geno, center=TRUE, scale=TRUE)

save(W, file="../DGRPdata/dgrp2_200lines_common_scaled.RData")

