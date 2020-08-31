########################
## Fabio Morgante
## 7/15/2018
## Obtain genes within the most predictive GOs from GO-TBLUP
##
## Modified on 07/30/2018
## Reason: Added header=TRUE, as.is=TRUE when loading accuracies, just to be safe
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

###Get values of arguments
args <- commandArgs(TRUE)
trait <- as.numeric(args)
ntopGOs <- 3

###Load the data
##Annotations
load("Annotations/go2fb_F.RData")

##Phenotypes (ONLY to get names for the following step)
pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)

##Read in results of GO-TBLUP
accuracies <- read.csv(paste('Results/Pred_Expression_GO_regress_adjPheno_common_', colnames(pheno)[trait],'_females.csv', sep=""), header=TRUE)

###Order results by mean r
accuracies_ord <- accuracies[order(accuracies$mean.r, decreasing=TRUE),]

###Get top GO names
name_topGOS <- as.character(accuracies_ord[1:ntopGOs, 1])

###Extract genes within top GOs
topGOs <- go2fb_F[name_topGOS]

###Draw Venn diagram
# venn.diagram(
# x = topGOs,
# category.names = names(topGOs),
# filename = 'prova.png',
#         output = TRUE ,
#         imagetype="png" ,
#         height = 480 , 
#         width = 480 , 
#         resolution = 300,
#         compression = "lzw",
#         lwd = 2,
#         lty = 'blank',
#         fill = c('yellow', 'purple', 'green'),
#         cat.cex = 0.4,
#         rotation = 1
#         )

###Calculate overlap of the 3 sets
overlap <- calculate.overlap(topGOs)
names(overlap) <- c(paste0(names(topGOs)[1], "-", names(topGOs)[2], "-", names(topGOs)[3]), 
					paste0(names(topGOs)[1], "-", names(topGOs)[2]),
					paste0(names(topGOs)[1], "-", names(topGOs)[3]),
					paste0(names(topGOs)[2], "-", names(topGOs)[3]),
					names(topGOs)[1],
					names(topGOs)[2],
					names(topGOs)[3])

###Arrange the overlap list into a matrix (to write to a file)
overlap_matrix <- sapply(overlap, '[', seq(max(sapply(overlap,length))))

###Write results to a file
write.table(overlap_matrix, paste('Results/Overlap_Expression_top3GO_regress_adjPheno_common_', colnames(pheno)[trait],'_females.csv', sep=""), 
		  col.names=T, row.names=F, quote=F, sep=",")

