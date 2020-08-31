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


#################
#####Females#####
#################

###Load the data
##GO-TBLUP results
go_tblup_F <- read.table(paste("Results/Pred_Expression_GO_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')

###Add sex label
go_tblup_F$Sex <- "Females"


###############
#####Males#####
###############

###Load the data
##GO-TBLUP results
go_tblup_M <- read.table(paste("Results/Pred_Expression_GO_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')

###Add sex label
go_tblup_M$Sex <- "Males"



####Merge females and males data
go_tblup_all <- rbind(go_tblup_F, go_tblup_M)



####Scatterplot
p <- ggplot(go_tblup_all, aes(x=mean.h2go, y=mean.r, colour=Sex))

scatter <- p + geom_point() +  geom_smooth(method=lm, colour="black")+
  labs(x = "Mean proportion of variance explained (in training set)", y = expression(paste("mean", sep=" ", italic(r), " (in test set)"))) +
  facet_grid(. ~ Sex) +
  ggtitle(paste(trait)) +
  theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=16), axis.title=element_text(size=18), axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
        strip.text = element_text(size=16), legend.position="none", plot.title=element_text(size=20, hjust=0.5))


##Save to a pdf
pdf(paste('Results/Graphs/Accuracy_vs_VarianceExplained_GO-transcripts_adjPheno_common_', trait,'.pdf', sep=""), width=15, height=11)
print(scatter)
dev.off()


