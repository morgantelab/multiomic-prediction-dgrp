############################
###Fabio Morgante		 
###5-04-2017			 
###Extract data for genes
###with mean expression  
###> -1.828223 (from Logan's  
###analysis)			 
###
###Modified on 4-2-2018
###Reason: use expression data adjusted ONLY for alignment bias after Logan fixed the issue
###
###Modified on 11/11/2018
###Reason: add options(stringsAsFactors=FALSE) and remove as.is=T        
####
####Modified on 05/06/2020
####Reason: Use relative paths             
#############################

###Set number of cores for R open
setMKLthreads(5)

###Set global options
options(stringsAsFactors=FALSE)

###Source custom functions
source(file="../Software/myRfunctions.R")

###############
####Females####
###############

###Load the data
data_F <- read.table("Data/combined_samples_known_novel_fpkm_VR_NoWol_F_line_means_rename_noflag_transp.txt",
					header=T, row.names=1, sep="\t")

###Calculate summary statistics
mean_F <- apply(data_F, 2, mean)
names(mean_F) <- colnames(data_F)

###Extract genes with mean expression > -1.828223
namesToKeep_F <- names(mean_F[which(mean_F>-1.828223)])
data_Fkeep <- data_F[, colnames(data_F) %in% namesToKeep_F]
line <- row.names(data_F)

###Write to file
###Write the new data, transposed and without flag
write.table(cbind(line, data_Fkeep), "Data/combined_samples_known_novel_fpkm_VR_NoWol_F_line_means_rename_noflag_transp_common.txt", 
			sep="\t", col.names=T, row.names=F, quote=F)


rm(list=ls())



#############
####Males####
#############

###Load the data
data_M <- read.table("Data/combined_samples_known_novel_fpkm_VR_NoWol_M_line_means_rename_noflag_transp.txt",
					header=T, row.names=1, sep="\t")

###Calculate summary statistics
mean_M <- apply(data_M, 2, mean)
names(mean_M) <- colnames(data_M)

###Extract genes with mean expression > -1.828223
namesToKeep_M <- names(mean_M[which(mean_M>-1.828223)])
data_Mkeep <- data_M[, colnames(data_M) %in% namesToKeep_M]
line <- row.names(data_M)

###Write to file
###Write the new data, transposed and without flag
write.table(cbind(line, data_Mkeep), "Data/combined_samples_known_novel_fpkm_VR_NoWol_M_line_means_rename_noflag_transp_common.txt", 
			sep="\t", col.names=T, row.names=F, quote=F)

