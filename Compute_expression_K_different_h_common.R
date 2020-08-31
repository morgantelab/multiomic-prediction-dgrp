################################
###Fabio Morgante		     
###4-3-2018			         
###Compute of several kernels
###based on gene expression  
###
###Modified on 11/11/2018
###Reason: add options(stringsAsFactors=FALSE) and remove as.is=T        
###
###Modified on 05/06/2020
###Reason: Use relative paths             
#############################

###Set number of cores for R open
setMKLthreads(5)

###Set global options
options(stringsAsFactors=FALSE)

###############
####Females####
###############

###Load the data
D_F <- as.matrix(read.table("Matrices/D_matrix_common_females.txt",
					header=T, row.names=1, sep=" "))

###Compute kernels
##Choose bandwith parameter h
h <- 0.5 * c(1/5, 1)
##Compute the kernels
for(i in 1:length(h)){
	K_F <- exp(-h[i]*D_F)
	#Write to a file
	write.table(K_F, paste("Matrices/K_matrix_h_", h[i], "_common_females.txt", sep=""), sep=" ", quote=F, row.names=T, col.names=T)
}




rm(list=ls())


#############
####Males####
#############

###Load the data
D_M <- as.matrix(read.table("Matrices/D_matrix_common_males.txt",
					header=T, row.names=1, sep=" "))

###Compute kernels
##Choose bandwith parameter h
h <- 0.5 * c(1/5, 1)
##Compute the kernels
for(i in 1:length(h)){
	K_M <- exp(-h[i]*D_M)
	#Write to a file
	write.table(K_M, paste("Matrices/K_matrix_h_", h[i], "_common_males.txt", sep=""), sep=" ", quote=F, row.names=T, col.names=T)
}















