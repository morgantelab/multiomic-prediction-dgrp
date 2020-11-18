########################
## Fabio Morgante
## 6/1/2020
## Plot var exp vs accuracy
## for GO-GBLUP 
#########################

###Set number of threads
setMKLthreads(1)

###Set global options
options(stringsAsFactors=FALSE)

###Load libraries
library(ggplot2)
library(ggpubr)

###Loop over traits
for(i in 10){
	###Select trait to use
	traitNum <- i

	###Load phenotypes (needed only to get trait name!)
	pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)

	trait <- colnames(pheno)[traitNum]

	###Load dimension of GOs
	go_size <- read.table("Annotations/go2variants_length.txt", sep="\t", header=T)


	###Females
	##Load GO-GBLUP results
	go_gblup_F <- read.table(paste("Results/Pred_SNPs_GO_regress_adjPheno_", trait, "_females.csv", sep=""), header=T, sep=',')

	##Add sex label
	go_gblup_F$Sex <- "Females"


	###Males
	##Load GO-GBLUP results
	go_gblup_M <- read.table(paste("Results/Pred_SNPs_GO_regress_adjPheno_", trait, "_males.csv", sep=""), header=T, sep=',')

	##Add sex label
	go_gblup_M$Sex <- "Males"


	###Merge females and males data
	go_gblup_all <- rbind(go_gblup_F, go_gblup_M)



	###Plot the data
	p <- ggplot(go_gblup_all, aes(x=mean.h2go, y=mean.r, colour=Sex))

	scatter <- p + geom_point(size=1) + geom_smooth(method="lm", colour="black")+
	  stat_cor(method = "pearson", cor.coef.name ="r", label.x=0, label.y=0.3, size=3.7)+
	  labs(x = "mean proportion of variance explained (in training set)", y = expression(paste("mean", sep=" ", italic(r), " (in test set)"))) +
	  facet_wrap(~ Sex, scales="free", ncol=2) +
	  scale_colour_manual(values=c("red", "blue")) +
	  theme_bw() +
	  theme(axis.text.y=element_text(size=10), axis.text.x=element_text(size=10), axis.title=element_text(size=12), axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
			strip.text = element_text(size=10), legend.position="none")


	##Save to a pdf
	jpeg(paste('Results/Graphs/Manuscript/Accuracy_vs_VarianceExplained_GO-SNPs_adjPheno_', trait,'_MS.jpeg', sep=""), width=20, height=15, units="cm", res=300)
	print(scatter)
	dev.off()

	rm(list=ls())
}

