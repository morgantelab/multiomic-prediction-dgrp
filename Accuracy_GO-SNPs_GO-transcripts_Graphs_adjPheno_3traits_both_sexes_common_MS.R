########################
## Fabio Morgante
## 6/1/2020
## Plot accuracy for GO-GTBLUP   
#########################

###Set number of threads
setMKLthreads(1)

###Set global options
options(stringsAsFactors=FALSE)

###Load libraries
library(ggplot2)

###Loop over traits
for(i in c(3, 9, 10)){
	###Select trait to use
	traitNum <- i

	###Load phenotypes (needed only to get trait name!)
	pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)

	trait <- colnames(pheno)[traitNum]

	###Load and reshape the data
	##Females
	GO_trait_F <- read.table(paste("Results/Pred_SNPs+Expression_GO_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	GO_trait_F <- GO_trait_F[, c(1, 12:15)]

	Trait_F <- rep(trait, times=nrow(GO_trait_F))
	Sex_F <- rep('Females', times=nrow(GO_trait_F))
	rank_F <- rank(-GO_trait_F$mean.r)
	data_F <- cbind(GO_trait_F, rank_F, Trait_F, Sex_F)
	colnames(data_F)[6:8] <- c("Rank", "Trait", "Sex")

	##Males
	GO_trait_M <- read.table(paste("Results/Pred_SNPs+Expression_GO_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	GO_trait_M <- GO_trait_M[, c(1, 12:15)]

	Trait_M <- rep(trait, times=nrow(GO_trait_M))
	Sex_M <- rep('Males', times=nrow(GO_trait_M))
	rank_M <- rank(-GO_trait_M$mean.r)
	data_M <- cbind(GO_trait_M, rank_M, Trait_M, Sex_M)
	colnames(data_M)[6:8] <- c("Rank", "Trait", "Sex")

	###Build the final data for plotting
	dat <- rbind(data_F, data_M)


	###Plot the data
	p <- ggplot(data = dat, aes(x=GO, y = mean.r, colour=Sex))

	gPred <- p +  geom_point(shape=16, size=1.5) + 
	  labs(x = "GO", y = expression(paste("mean", sep=" ", italic(r)))) +
	  scale_x_discrete(breaks=seq(from=1, to=(nrow(dat)/2), by=1)) +
	  facet_grid(Sex ~ ., scales="free_y") +
	  geom_text(aes(label=ifelse((Rank==1|Rank==2|Rank==3), as.character(GO), '')), hjust=0.3, vjust=0, size=4) +
	  scale_colour_manual(values=c("red", "blue")) +
	  theme_bw()+
	  theme(axis.text.y=element_text(size=16), axis.title=element_text(size=18), axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
			axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text = element_text(size=16), legend.position="none")

	##Save to a pdf
	pdf(paste('Results/Graphs/Manuscript/Accuracy_GO-SNPs_GO-transcripts_adjPheno_common_', trait,'_MS.pdf', sep=""), width=15, height=11)
	print(gPred)
	dev.off()

	rm(list=ls())

}
