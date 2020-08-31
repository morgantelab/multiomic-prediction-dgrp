########################
## Fabio Morgante
## 6/1/2020
## Plot GO size vs accuracy
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
for(i in c(3, 9, 10)){
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

	##Merge datasets
	go_gblup_size_F <- merge(go_gblup_F, go_size, by="GO", all.x=T)

	##Add sex label
	go_gblup_size_F$Sex <- "Females"


	###Males
	##Load GO-GBLUP results
	go_gblup_M <- read.table(paste("Results/Pred_SNPs_GO_regress_adjPheno_", trait, "_males.csv", sep=""), header=T, sep=',')

	##Merge datasets
	go_gblup_size_M <- merge(go_gblup_M, go_size, by="GO", all.x=T)

	##Add sex label
	go_gblup_size_M$Sex <- "Males"


	###Merge females and males data
	go_gblup_size <- rbind(go_gblup_size_F, go_gblup_size_M)


	###Plot the data
	p <- ggplot(go_gblup_size, aes(x=variants.length, y=mean.r, colour=Sex))

	scatter <- p + geom_point() +  geom_smooth(method="lm", colour="black") +
	  stat_cor(method = "pearson", cor.coef.name ="r", label.x=50000, size=6) +
	  labs(x = "GO size (number of variants)", y = expression(paste("mean", sep=" ", italic(r)))) +
	  facet_wrap(~ Sex, scales="free", ncol=2) +
	  scale_colour_manual(values=c("red", "blue")) +
	  theme_bw() +
	  theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=16), axis.title=element_text(size=18), axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
			strip.text = element_text(size=16), legend.position="none", plot.title=element_text(size=20, hjust=0.5))


	##Save to a pdf
	pdf(paste('Results/Graphs/Manuscript/Accuracy_vs_GOsize_GO-SNPs_adjPheno_', trait,'_MS.pdf', sep=""), width=15, height=11)
	print(scatter)
	dev.off()

	rm(list=ls())
}

