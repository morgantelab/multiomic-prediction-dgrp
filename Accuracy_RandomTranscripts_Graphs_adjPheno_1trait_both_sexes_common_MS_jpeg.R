########################
## Fabio Morgante
## 6/1/2020
## Plot accuracy for TBLUP with
## random transcript selection   
#########################

###Set number of threads
setMKLthreads(1)

###Set global options
options(stringsAsFactors=FALSE)

###Load libraries
library(ggplot2)

###Loop over traits
for(i in 10){
	###Select trait to use
	traitNum <- i

	###Load phenotypes (needed only to get trait name!)
	pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)

	trait <- colnames(pheno)[traitNum]

	###Load and reshape the data
	##Females
	AllTrans_trait_F <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',', nrows=1)
	#colnames(AllTrans_trait_F)[10:11] <- c("mean.h2", "se.h2")
	rnd_5000_trait_F <- read.table(paste("Results/Pred_Expression_Random_5000_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	rnd_1000_trait_F <- read.table(paste("Results/Pred_Expression_Random_1000_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	rnd_500_trait_F <- read.table(paste("Results/Pred_Expression_Random_500_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	rnd_50_trait_F <- read.table(paste("Results/Pred_Expression_Random_50_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	rnd_5_trait_F <- read.table(paste("Results/Pred_Expression_Random_5_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')

	#Add a column to indicate which transcripts were used
	AllTrans_trait_F$Transcripts <- "All"
	rnd_5000_trait_F$Transcripts <- "Random 5000"
	rnd_1000_trait_F$Transcripts <- "Random 1000"
	rnd_500_trait_F$Transcripts <- "Random 500"
	rnd_50_trait_F$Transcripts <- "Random 50"
	rnd_5_trait_F$Transcripts <- "Random 5"

	#Add a column to indicate which trait
	AllTrans_trait_F$Trait <- trait
	rnd_5000_trait_F$Trait <- trait
	rnd_1000_trait_F$Trait <- trait
	rnd_500_trait_F$Trait <- trait
	rnd_50_trait_F$Trait <- trait
	rnd_5_trait_F$Trait <- trait

	#Bind the data
	data_trait_F <- rbind(AllTrans_trait_F, rnd_5000_trait_F, rnd_1000_trait_F, rnd_500_trait_F, rnd_50_trait_F, rnd_5_trait_F)
	data_trait_F$Sex <- "Females"

	##Males
	AllTrans_trait_M <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',', nrows=1)
	#colnames(AllTrans_trait_M)[10:11] <- c("mean.h2", "se.h2")
	rnd_5000_trait_M <- read.table(paste("Results/Pred_Expression_Random_5000_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	rnd_1000_trait_M <- read.table(paste("Results/Pred_Expression_Random_1000_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	rnd_500_trait_M <- read.table(paste("Results/Pred_Expression_Random_500_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	rnd_50_trait_M <- read.table(paste("Results/Pred_Expression_Random_50_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	rnd_5_trait_M <- read.table(paste("Results/Pred_Expression_Random_5_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')

	#Add a column to indicate which transcripts were used
	AllTrans_trait_M$Transcripts <- "All"
	rnd_5000_trait_M$Transcripts <- "Random 5000"
	rnd_1000_trait_M$Transcripts <- "Random 1000"
	rnd_500_trait_M$Transcripts <- "Random 500"
	rnd_50_trait_M$Transcripts <- "Random 50"
	rnd_5_trait_M$Transcripts <- "Random 5"

	#Add a column to indicate which trait
	AllTrans_trait_M$Trait <- trait
	rnd_5000_trait_M$Trait <- trait
	rnd_1000_trait_M$Trait <- trait
	rnd_500_trait_M$Trait <- trait
	rnd_50_trait_M$Trait <- trait
	rnd_5_trait_M$Trait <- trait

	#Bind the data
	data_trait_M <- rbind(AllTrans_trait_M, rnd_5000_trait_M, rnd_1000_trait_M, rnd_500_trait_M, rnd_50_trait_M, rnd_5_trait_M)
	data_trait_M$Sex <- "Males"


	###Create data for final plot
	dat <- rbind(data_trait_F, data_trait_M)


	###Plot the data
	limits <- aes(ymax = dat$mean.r + dat$se.r,
				  ymin = dat$mean.r - dat$se.r)

	p <- ggplot(data = dat, aes(x = Trait, y = mean.r,
								 fill = factor(Transcripts, levels=c("All", "Random 5000", "Random 1000", "Random 500", "Random 50", "Random 5"))))


	gPred <- p + geom_bar(stat = "identity",
						  position = position_dodge(0.9)) +
	  geom_errorbar(limits, position = position_dodge(0.9), width = 0.25) +
	  #geom_text(aes(label=nreps), position=position_dodge(width=0.9), vjust=-5.25)+
	  labs(y = expression(paste("mean", sep=" ", italic(r)))) +
	  facet_grid(. ~ Sex, scales="free_y") +
	  theme_bw()+
	  theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
			axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")), legend.text=element_text(size=10), legend.title=element_text(size=12), strip.text = element_text(size=10), legend.key.size = unit(2, 'lines')) +
	  scale_fill_brewer(name = "Transcripts",palette="Set1")

	##Save to a pdf
	jpeg(paste('Results/Graphs/Manuscript/Accuracy_Expression_Random_common_', trait, '_MS.jpeg', sep=""), width=20, height=15, units="cm", res=300)
	print(gPred)
	dev.off()

	rm(list=ls())
}








