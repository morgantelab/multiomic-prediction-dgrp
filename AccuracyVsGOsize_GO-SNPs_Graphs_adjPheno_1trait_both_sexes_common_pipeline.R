########################
## Fabio Morgante
## 5/1/2018
## Correlate prediction accuracy with GO dimension
########################

###Set number of cores for R open
setMKLthreads(2)

###Set global options
options(stringsAsFactors=FALSE)

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

###Merge datasets
go_gblup_size_F <- merge(go_gblup_F, go_size, by="GO", all.x=T)

###Add sex label
go_gblup_size_F$Sex <- "Females"

###Compute correlation between r and size of the GO
print(cor.test(go_gblup_size_F$mean.r, go_gblup_size_F$variants.length, method=c("kendall")))



###############
#####Males#####
###############

###Load the data
##GO-TBLUP results
go_gblup_M <- read.table(paste("Results/Pred_SNPs_GO_regress_adjPheno_", trait, "_males.csv", sep=""), header=T, sep=',')

###Merge datasets
go_gblup_size_M <- merge(go_gblup_M, go_size, by="GO", all.x=T)

###Add sex label
go_gblup_size_M$Sex <- "Males"

###Compute correlation between r and size of the GO
print(cor.test(go_gblup_size_M$mean.r, go_gblup_size_M$variants.length, method=c("kendall")))



####Merge females and males data
go_gblup_size <- rbind(go_gblup_size_F, go_gblup_size_M)



####Scatterplot
p <- ggplot(go_gblup_size, aes(x=variants.length, y=mean.r, colour=Sex))

scatter <- p + geom_point() +  geom_smooth(method=lm, colour="black")+
  labs(x = "GO size (number of variants)", y = expression(paste("mean", sep=" ", italic(r)))) +
  facet_grid(. ~ Sex) +
  ggtitle(paste(trait)) +
  theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=16), axis.title=element_text(size=18), axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
        strip.text = element_text(size=16), legend.position="none", plot.title=element_text(size=20, hjust=0.5))


##Save to a pdf
pdf(paste('Results/Graphs/Accuracy_vs_GOsize_GO-SNPs_adjPheno_', trait,'.pdf', sep=""), width=15, height=11)
print(scatter)
dev.off()


