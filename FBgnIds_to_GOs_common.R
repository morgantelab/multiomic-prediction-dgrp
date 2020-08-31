########################
## Fabio Morgante
## 7/12/2017
## Link FBgn ids to GO terms
## 
## Modified on 4-20-2018
## Reason: use expression data adjusted ONLY for alignment bias after Logan fixed the issue
##
## Modified on 11/11/2018
## Reason: add options(stringsAsFactors=FALSE) and remove as.is=T        
##
## Modified on 05/06/2020
## Reason: Use relative paths             
#############################

###Set number of cores for R open
setMKLthreads(1)

###Set global options
options(stringsAsFactors=FALSE)

###Install and load libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Dm.eg.db")
library(org.Dm.eg.db)

###Link fb gene ids to entrez gene ids.
fb2eg <- org.Dm.egFLYBASE2EG
mapped_genes <- mappedkeys(fb2eg)
fb2eg <- as.list(fb2eg[mapped_genes])

###Link entrez gene ids to fb gene ids.
eg2fb <- org.Dm.egFLYBASE
mapped_genes <- mappedkeys(eg2fb)
eg2fb <- as.list(eg2fb[mapped_genes])

###Link gene GO ids to entrez gene ids.
go2eg <- as.list(org.Dm.egGO2EG)

###Link fb gene ids to GO ids.
##Females
#Load transcript data
dataF <- as.matrix(read.table("Data/combined_samples_known_novel_fpkm_VR_NoWol_F_line_means_rename_noflag_transp_common.txt",
                              header=T, row.names=1, sep="\t"))

#Link data to GO and remove empty elements
go2fb_F <- lapply(go2eg,function(x){ 
  fb <- na.omit(unlist(eg2fb[x]))
  m <- match(fb, colnames(dataF))
  fb <- fb[!is.na(m)]
  fb <- fb[!duplicated(fb)]
  return(fb)
})

go2fb_F <- go2fb_F[sapply(go2fb_F,length)>0]

#Save GO annotations to a file 
save(go2fb_F,file="Annotations/go2fb_F.RData")

#Compute number of genes within each GO and write it out
go2fb_F_length <- sapply(go2fb_F,length)
go2fb_F_length_df <- data.frame(names(go2fb_F_length), go2fb_F_length)
colnames(go2fb_F_length_df) <- c("GO", "length")
write.table(go2fb_F_length_df, "Annotations/go2fb_F_length.txt", sep="\t", row.names=F, col.names=T, quote=F)


##Males
#Load transcript data
dataM <- as.matrix(read.table("Data/combined_samples_known_novel_fpkm_VR_NoWol_M_line_means_rename_noflag_transp_common.txt",
                              header=T, row.names=1, sep="\t"))

#Link data to GO and remove empty elements
go2fb_M <- lapply(go2eg,function(x){ 
  fb <- na.omit(unlist(eg2fb[x]))
  m <- match(fb, colnames(dataM))
  fb <- fb[!is.na(m)]
  fb <- fb[!duplicated(fb)]
  return(fb)
})

go2fb_M <- go2fb_M[sapply(go2fb_M,length)>0]

#Save GO annotations to a file 
save(go2fb_M,file="Annotations/go2fb_M.RData")

#Compute number of genes within each GO and write it out
go2fb_M_length <- sapply(go2fb_M,length)
go2fb_M_length_df <- data.frame(names(go2fb_M_length), go2fb_M_length)
colnames(go2fb_M_length_df) <- c("GO", "length")
write.table(go2fb_M_length_df, "Annotations/go2fb_M_length.txt", sep="\t", row.names=F, col.names=T, quote=F)
