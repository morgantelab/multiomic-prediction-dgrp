################################
###Fabio Morgante		     
###4-3-2018			         
###Compute of several kernels
###based on gene expression  
###
###Modified on 11/11/2018
###Reason: add options(stringsAsFactors=FALSE) and remove as.is=T        
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
					
###Correlation matrix
cor_mat_F <- cor(t(dataF))
write.table(cor_mat_F, "Matrices/exp_cor_matrix_common_females.txt", sep=" ", quote=F, row.names=T, col.names=T)

###TRM
X_F <- scale(dataF, center=TRUE, scale=TRUE)
TRM_F <- tcrossprod(X_F)/ncol(X_F)
write.table(TRM_F, "Matrices/TRM_common_females.txt", sep=" ", quote=F, row.names=T, col.names=T)

###Distance matrix (base for RKHS)
D_F <- (as.matrix(dist(X_F, method='euclidean'))^2)/ncol(X_F)
write.table(D_F, "Matrices/D_matrix_common_females.txt", sep=" ", quote=F, row.names=T, col.names=T)



rm(list=ls())


#############
####Males####
#############

###Load the data
dataM <- as.matrix(read.table("Data/combined_samples_known_novel_fpkm_VR_NoWol_M_line_means_rename_noflag_transp_common.txt",
					header=T, row.names=1, sep="\t"))
					
###Correlation matrix
cor_mat_M <- cor(t(dataM))
write.table(cor_mat_M, "Matrices/exp_cor_matrix_common_males.txt", sep=" ", quote=F, row.names=T, col.names=T)

###TRM
X_M <- scale(dataM, center=TRUE, scale=TRUE)
TRM_M <- tcrossprod(X_M)/ncol(X_M)
write.table(TRM_M, "Matrices/TRM_common_males.txt", sep=" ", quote=F, row.names=T, col.names=T)

###Distance matrix (base for RKHS)
D_M <- (as.matrix(dist(X_M, method='euclidean'))^2)/ncol(X_M)
write.table(D_M, "Matrices/D_matrix_common_males.txt", sep=" ", quote=F, row.names=T, col.names=T)















