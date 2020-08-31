########################
## Fabio Morgante
## 6/1/2020
## Plot GO size vs accuracy
## for GO-TBLUP 
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


	###Females
	##Load GO-TBLUP results
	go_tblup_F <- read.table(paste("Results/Pred_Expression_GO_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')

	##Load dimension of GOs
	go_size_F <- read.table("Annotations/go2fb_F_length.txt", sep="\t", header=T)

	##Merge datasets
	go_tblup_size_F <- merge(go_tblup_F, go_size_F, by="GO", all.x=T)

	##Add sex label
	go_tblup_size_F$Sex <- "Females"

	###Males
	##Load GO-TBLUP results
	go_tblup_M <- read.table(paste("Results/Pred_Expression_GO_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')

	##Load dimension of GOs
	go_size_M <- read.table("Annotations/go2fb_M_length.txt", sep="\t", header=T)

	##Merge datasets
	go_tblup_size_M <- merge(go_tblup_M, go_size_M, by="GO", all.x=T)

	##Add sex label
	go_tblup_size_M$Sex <- "Males"


	###Merge females and males data
	go_tblup_size <- rbind(go_tblup_size_F, go_tblup_size_M)


	###Plot the data
	p <- ggplot(go_tblup_size, aes(x=length, y=mean.r, colour=Sex))

	scatter <- p + geom_point() + geom_smooth(method="lm", colour="black")+
	stat_cor(method = "pearson", cor.coef.name ="r", label.x=400, size=6)+
	  labs(x = "GO size (number of genes)", y = expression(paste("mean", sep=" ", italic(r)))) +
	  facet_wrap(~ Sex, scales="free", ncol=2) +
	  scale_colour_manual(values=c("red", "blue")) +
	  theme_bw() +
	  theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=16), axis.title=element_text(size=18), axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
			strip.text = element_text(size=16), legend.position="none", plot.title=element_text(size=20, hjust=0.5))


	##Save to a pdf
	pdf(paste('Results/Graphs/Manuscript/Accuracy_vs_GOsize_GO-transcripts_adjPheno_common_', trait,'_MS.pdf', sep=""), width=15, height=11)
	print(scatter)
	dev.off()

	rm(list=ls())
}

