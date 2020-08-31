########################
## Fabio Morgante
## 7/30/2018
## Obtain genes and variants within the most predictive GOs from GO-GBLUP
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

###Load libraries
library(VennDiagram)
library(plyr)

###Get values of arguments
args <- commandArgs(TRUE)
trait <- as.numeric(args)
ntopGOs <- 3

###Load the data
##Annotations
load("Annotations/go2variants.RData")

##Phenotypes (ONLY to get names for the following step)
pheno <- read.csv("Data/eQTL_traits_males.csv", header=TRUE)

##Read in results of GO-TBLUP
accuracies <- read.csv(paste('Results/Pred_SNPs_GO_regress_adjPheno_', colnames(pheno)[trait],'_males.csv', sep=""), header=TRUE)

###Order results by mean r
accuracies_ord <- accuracies[order(accuracies$mean.r, decreasing=TRUE),]

###Get top GO names
name_topGOs <- as.character(accuracies_ord[1:ntopGOs, 1])

###Extract genes within top GOs
topGOs <- go2variants[name_topGOs]

####Gene level overlap####
###Extract gene names within top GOs
genes_topGOs <- list(names(topGOs[[1]]), names(topGOs[[2]]), names(topGOs[[3]]))
names(genes_topGOs) <- names(topGOs)

###Calculate overlap of the 3 sets
overlap <- calculate.overlap(genes_topGOs)
names(overlap) <- c(paste0(names(genes_topGOs)[1], "-", names(genes_topGOs)[2], "-", names(genes_topGOs)[3]), 
					paste0(names(genes_topGOs)[1], "-", names(genes_topGOs)[2]),
					paste0(names(genes_topGOs)[1], "-", names(genes_topGOs)[3]),
					paste0(names(genes_topGOs)[2], "-", names(genes_topGOs)[3]),
					names(genes_topGOs)[1],
					names(genes_topGOs)[2],
					names(genes_topGOs)[3])

###Arrange the overlap list into a matrix (to write to a file)
overlap_matrix <- sapply(overlap, '[', seq(max(sapply(overlap,length))))

###Write results to a file
write.table(overlap_matrix, paste('Results/Overlap_SNPs_gene_level_top3GO_regress_adjPheno_common_', colnames(pheno)[trait],'_males.csv', sep=""), 
		  col.names=T, row.names=F, quote=F, sep=",")
		  
		  
		  
####Variant level overlap####		  
###Unlist the second level to obtain a list of variants instead of genes
variants_topGOs <- lapply(llply(topGOs, unlist), unname)

###Keep only unique variants (variants may be in more than one gene)
variants_topGOs_unique <- lapply(variants_topGOs, unique)

###Calculate overlap of the 3 sets
overlap_variants <- calculate.overlap(variants_topGOs_unique)
names(overlap_variants) <- c(paste0(names(variants_topGOs_unique)[1], "-", names(variants_topGOs_unique)[2], "-", names(variants_topGOs_unique)[3]), 
					paste0(names(variants_topGOs_unique)[1], "-", names(variants_topGOs_unique)[2]),
					paste0(names(variants_topGOs_unique)[1], "-", names(variants_topGOs_unique)[3]),
					paste0(names(variants_topGOs_unique)[2], "-", names(variants_topGOs_unique)[3]),
					names(variants_topGOs_unique)[1],
					names(variants_topGOs_unique)[2],
					names(variants_topGOs_unique)[3])

###Arrange the overlap list into a matrix (to write to a file)
overlap_variants_matrix <- sapply(overlap_variants, '[', seq(max(sapply(overlap_variants,length))))

###Write results to a file
write.table(overlap_variants_matrix, paste('Results/Overlap_SNPs_variant_level_top3GO_regress_adjPheno_common_', colnames(pheno)[trait],'_males.csv', sep=""), 
		  col.names=T, row.names=F, quote=F, sep=",")




