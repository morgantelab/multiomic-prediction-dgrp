#############################
#### Fabio Morgante      
#### 3-27-2018           
#### Rename expression   
#### data (adjusted ONLY for alignment bias) 
#### columns
####
####Modified on 11/11/2018
####Reason: add options(stringsAsFactors=FALSE) and remove as.is=T
####
####Modified on 04/27/2019
####Reason: use new expression data produced by Wen after first round of review of the Nature MS             
####
####Modified on 05/06/2020
####Reason: Use relative paths             
#############################

setMKLthreads(1)

options(stringsAsFactors=FALSE)

###############
####Females####
###############

###Load the data
dataF <- read.table("/home/whuang9/dgrpRNASeq/female.tx.exp.txt", header=T, sep="\t", check.names=F)

###Remove _F from column names and add line_ as prefix
colnames(dataF) <- gsub("_F", "", colnames(dataF))
colnames(dataF)[3:ncol(dataF)] <- paste("line_", colnames(dataF)[3:ncol(dataF)], sep = "")

###Write the new data
write.table(dataF, "Data/combined_samples_known_novel_fpkm_VR_NoWol_F_line_means_rename.txt", 
			sep="\t", col.names=T, row.names=F, quote=F)
			

###Clear workspace			
rm(list=ls())		


###Load the data
dataF <- read.table("/home/whuang9/dgrpRNASeq/female.tx.exp.txt", 
					header=T, sep="\t", check.names=F, row.names=1)

###Remove _F from column names and add line_ as prefix
colnames(dataF) <- gsub("_F", "", colnames(dataF))
colnames(dataF)[2:ncol(dataF)] <- paste("line_", colnames(dataF)[2:ncol(dataF)], sep = "")

###Transpose and remove flag
dataF_t <- t(dataF[, -1])

###Add a column with line names
dataF_t <- data.frame(rownames(dataF_t), dataF_t)
colnames(dataF_t)[1] <- c("line")

			
###Write the new data, transposed and without flag
write.table(dataF_t, "Data/combined_samples_known_novel_fpkm_VR_NoWol_F_line_means_rename_noflag_transp.txt", 
			sep="\t", col.names=T, row.names=F, quote=F)
	
###Clear workspace			
rm(list=ls())		

	
	
	
	
#############					
####Males####
#############

###Load the data
dataM <- read.table("/home/whuang9/dgrpRNASeq/male.tx.exp.txt", header=T, sep="\t", check.names=F)

###Remove _M from column names and add line_ as prefix
colnames(dataM) <- gsub("_M", "", colnames(dataM))
colnames(dataM)[3:ncol(dataM)] <- paste("line_", colnames(dataM)[3:ncol(dataM)], sep = "")

###Write the new data
write.table(dataM, "Data/combined_samples_known_novel_fpkm_VR_NoWol_M_line_means_rename.txt", 
			sep="\t", col.names=T, row.names=F, quote=F)
			
###Clear workspace			
rm(list=ls())			
			
			
###Load the data
dataM <- read.table("/home/whuang9/dgrpRNASeq/male.tx.exp.txt", 
					header=T, sep="\t", check.names=F, row.names=1)

###Remove _M from column names and add line_ as prefix
colnames(dataM) <- gsub("_M", "", colnames(dataM))
colnames(dataM)[2:ncol(dataM)] <- paste("line_", colnames(dataM)[2:ncol(dataM)], sep = "")

###Transpose and remove flag
dataM_t <- t(dataM[, -1])

###Add a column with line names
dataM_t <- data.frame(rownames(dataM_t), dataM_t)
colnames(dataM_t)[1] <- c("line")

			
###Write the new data, transposed and without flag
write.table(dataM_t, "Data/combined_samples_known_novel_fpkm_VR_NoWol_M_line_means_rename_noflag_transp.txt", 
			sep="\t", col.names=T, row.names=F, quote=F)
				
			
			
			