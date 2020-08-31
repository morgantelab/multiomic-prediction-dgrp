################################
## Fabio Morgante		     
## 4-23-2018			         
## Compute W for gene expression 
##
## Modified on 11/11/2018
## Reason: add options(stringsAsFactors=FALSE) and remove as.is=T        
##
## Modified on 05/06/2020
## Reason: Use relative paths             
#############################


###Set number of cores for R open
setMKLthreads(5)

###Set global options
options(stringsAsFactors=FALSE)

###############
####Females####
###############

###Load the data
dataF <- as.matrix(read.table("Data/combined_samples_known_novel_fpkm_VR_NoWol_F_line_means_rename_noflag_transp_common.txt",
					header=T, row.names=1, sep="\t"))

###Compute W
W_F <- scale(dataF, center=TRUE, scale=TRUE)

save(W_F, file="Matrices/dgrp2_scaled_transcripts_common_females.RData")


rm(list=ls())


#############
####Males####
#############

###Load the data
dataM <- as.matrix(read.table("Data/combined_samples_known_novel_fpkm_VR_NoWol_M_line_means_rename_noflag_transp_common.txt",
					header=T, row.names=1, sep="\t"))
					
###Compute W
W_M <- scale(dataM, center=TRUE, scale=TRUE)

save(W_M, file="Matrices/dgrp2_scaled_transcripts_common_males.RData")













