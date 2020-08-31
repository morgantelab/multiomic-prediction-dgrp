##################################
## Fabio Morgante		       
## 4-23-2018			       
## Link SNPs to FBgn ids to GO terms	
##
## Modified on 7-28-2018
## Reason: Keep only unique SNPs within each GO (SNPs may be in multiple genes)	
## Note: Re-ran the last part manually --> creation dates of go2variants.RData and go2variants_length.txt will be different
##
## Modified on 11/11/2018
## Reason: add options(stringsAsFactors=FALSE) and remove as.is=T         
##
## Modified on 05/06/2020
## Reason: Use relative paths             
#############################

###Set number of threads for MKL
setMKLthreads(1)

###Set global options
options(stringsAsFactors=FALSE)

###Load required libraries
library(org.Dm.eg.db)
library(plyr)

###Load centered and scaled genotype matrix
load(file="../DGRPdata/dgrp2_200lines_common_scaled.RData")

###Read in the annotation file and rename rows
annotation <- read.table("../DGRPdata/dgrp.fb557.annot.txt", header=FALSE, colClasses="character")
rownames(annotation) <- as.character(annotation[, 1])

###Create a character vector with only the information we need and rename its elements
variantA0 <- as.character(annotation[, 3])
names(variantA0) <- rownames(annotation)

###Extract only info for site class 
variantA0 <- lapply(variantA0, function(x) {      
  x <- strsplit(x, ",")[[1]][1]
  x <- gsub("SiteClass", "", x)
  x <- gsub("[", "", x, fixed=TRUE)
  x <- gsub("]", "", x, fixed=TRUE)
  x <- strsplit(x, ";")[[1]]
  x
})

###Repeat each variant's name as many times as the number of genes it is in
nA0 <- sapply(variantA0, length)
variantNames <- rep(names(variantA0), times=nA0)
variantA0 <- unlist(variantA0, use.names=FALSE)

###Split records according to |
variantA0 <- lapply(variantA0, function(x) {unlist(strsplit(x, split="|", fixed=TRUE))})

###Store all the information into a matrix for variants that have info in all the 4 fields
variantA <- matrix(NA, nrow=length(variantNames), ncol=4)
rownames(variantA) <- variantNames
nA0 <- sapply(variantA0, length)
variantA[nA0==4, ] <- t(sapply(variantA0[nA0==4], function(x) {unlist(x)}))

###Extract only the variants that are present in W 
variantA <- variantA[which(rownames(variantA) %in% colnames(W)), ]
colnames(variantA) <- c("FBgnId", "GeneName", "SeqOnt", "distance")
#save(variantA, file="/home/fmorgan3/DGRPdata/variantA_W.RData")

###Get rid of variants without annotation and keep only those within genes
variantA <- na.omit(variantA)
inGene <- as.numeric(variantA[, 4])==0

###Create marker sets corresponding to each gene
fbSets <- unstack(rownames(variantA[inGene, ]), rownames(variantA[inGene, ]) ~ as.factor(variantA[inGene, 1]))

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
go2fb <- lapply(go2eg, function(x){ 
  fb <- na.omit(unlist(eg2fb[x]))
  m <- match(fb, names(fbSets)) 
  fb <- fb[!is.na(m)]
  fb <- fb[!duplicated(fb)]
  return(fb)
})

###Keep only GO ids with at least one FBgn from fbSets
go2fb <- go2fb[sapply(go2fb, length)>0]

###Link gene GO ids to variant names
go2variants <- lapply(go2fb, function(x){fbSets[x]})

save(go2variants,file="Annotations/go2variants.RData")

###Compute number of genes within each GO
go2variants_gene_length <- sapply(go2variants, length)

###Unlist the second level
go2variants_variants <- lapply(llply(go2variants, unlist), unname)

###Compute number of variants within each GO
go2variants_variants_unique <- lapply(go2variants_variants, unique)
go2variants_variants_length <- sapply(go2variants_variants_unique, length)

go2variants_length_df <- data.frame(names(go2variants_gene_length), go2variants_gene_length, go2variants_variants_length)
colnames(go2variants_length_df) <- c("GO", "gene.length", "variants.length")
write.table(go2variants_length_df, "Annotations/go2variants_length.txt", sep="\t", row.names=F, col.names=T, quote=F)



