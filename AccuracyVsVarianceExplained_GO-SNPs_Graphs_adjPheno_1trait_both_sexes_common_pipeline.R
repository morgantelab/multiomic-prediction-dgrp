########################
## Fabio Morgante
## 9/3/2018
## Correlate prediction accuracy with variance explained
## 
## Modified on 11/11/2018
## Reason: add options(stringsAsFactors=FALSE) and remove as.is=T        
########################

###Set options
options(stringsAsFactors=FALSE)

###Set number of cores for R open
setMKLthreads(2)

###Load libraries
library(ggplot2)
library(Rmisc)

###Select trait to analyze from command line arguments
args <- commandArgs(TRUE)
traitNum <- as.numeric(args)

#Load phenotypes (needed only to get trait name!)
pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)

trait <- colnames(pheno)[traitNum]


##Dimension of GOs
go_size <- read.table("Annotations/go2variants_length.txt", sep="\t", header=T)


#################
#####Females#####
#################

###Load the data
##GO-TBLUP results
go_gblup_F <- read.table(paste("Results/Pred_SNPs_GO_regress_adjPheno_", trait, "_females.csv", sep=""), header=T, sep=',')

###Add sex label
go_gblup_F$Sex <- "Females"


###############
#####Males#####
###############

###Load the data
##GO-TBLUP results
go_gblup_M <- read.table(paste("Results/Pred_SNPs_GO_regress_adjPheno_", trait, "_males.csv", sep=""), header=T, sep=',')

###Add sex label
go_gblup_M$Sex <- "Males"



####Merge females and males data
go_gblup_all <- rbind(go_gblup_F, go_gblup_M)



####Scatterplot
p <- ggplot(go_gblup_all, aes(x=mean.h2go, y=mean.r, colour=Sex))

scatter <- p + geom_point() +  geom_smooth(method=lm, colour="black")+
  labs(x = "Mean proportion of variance explained (in training set)", y = expression(paste("mean", sep=" ", italic(r), " (in test set)"))) +
  facet_grid(. ~ Sex) +
  ggtitle(paste(trait)) +
  theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=16), axis.title=element_text(size=18), axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
        strip.text = element_text(size=16), legend.position="none", plot.title=element_text(size=20, hjust=0.5))


##Save to a pdf
pdf(paste('Results/Graphs/Accuracy_vs_VarianceExplained_GO-SNPs_adjPheno_', trait,'.pdf', sep=""), width=15, height=11)
print(scatter)
dev.off()


