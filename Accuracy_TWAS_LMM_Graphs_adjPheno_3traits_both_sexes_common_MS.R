########################
## Fabio Morgante
## 6/1/2020
## Plot accuracy for TWAS-TBLUP   
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

	#Load phenotypes (needed only to get trait name!)
	pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)

	trait <- colnames(pheno)[traitNum]

	###Load and reshape the data
	##Females
	AllTrans_trait_F <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',', nrows=1)
	#weighted_trait_F <- read.table(paste("Results/Pred_Expression_TWAS_LMM_weighted_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	pvalue0.5_trait_F <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_0.5_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	pvalue0.1_trait_F <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_0.1_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	pvalue0.01_trait_F <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_0.01_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	pvalue0.001_trait_F <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_0.001_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	pvalue1e_4_trait_F <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_1e-04_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	pvalue1e_5_trait_F <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_1e-05_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	pvalue1e_6_trait_F <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_1e-06_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')

	#Add a column to indicate which transcripts were used
	AllTrans_trait_F$Transcripts <- "All"
	#weighted_trait_F$Transcripts <- "Weighted"
	pvalue0.5_trait_F$Transcripts <- "P<0.5"
	pvalue0.1_trait_F$Transcripts <- "P<0.1"
	pvalue0.01_trait_F$Transcripts <- "P<0.01"
	pvalue0.001_trait_F$Transcripts <- "P<0.001"
	pvalue1e_4_trait_F$Transcripts <- "P<0.0001"
	pvalue1e_5_trait_F$Transcripts <- "P<0.00001"
	pvalue1e_6_trait_F$Transcripts<- "P<0.000001"

	#Add a column to indicate which trait
	AllTrans_trait_F$Trait <- trait
	#weighted_trait_F$Trait <- trait
	pvalue0.5_trait_F$Trait <- trait
	pvalue0.1_trait_F$Trait <- trait
	pvalue0.01_trait_F$Trait <- trait
	pvalue0.001_trait_F$Trait <- trait
	pvalue1e_4_trait_F$Trait <- trait
	pvalue1e_5_trait_F$Trait <- trait
	pvalue1e_6_trait_F$Trait<- trait

	#Bind the data
	#data_trait_F <- rbind(AllTrans_trait_F, weighted_trait_F, pvalue0.5_trait_F, pvalue0.1_trait_F, pvalue0.01_trait_F, pvalue0.001_trait_F, pvalue1e_4_trait_F, pvalue1e_5_trait_F, pvalue1e_6_trait_F)
	data_trait_F <- rbind(AllTrans_trait_F, pvalue0.5_trait_F, pvalue0.1_trait_F, pvalue0.01_trait_F, pvalue0.001_trait_F, pvalue1e_4_trait_F, pvalue1e_5_trait_F, pvalue1e_6_trait_F)
	data_trait_F$Sex <- "Females"


	##Males
	AllTrans_trait_M <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',', nrows=1)
	#weighted_trait_M <- read.table(paste("Results/Pred_Expression_TWAS_LMM_weighted_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	pvalue0.5_trait_M <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_0.5_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	pvalue0.1_trait_M <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_0.1_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	pvalue0.01_trait_M <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_0.01_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	pvalue0.001_trait_M <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_0.001_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	pvalue1e_4_trait_M <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_1e-04_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	pvalue1e_5_trait_M <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_1e-05_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	pvalue1e_6_trait_M <- read.table(paste("Results/Pred_Expression_TWAS_LMM_pvalue_1e-06_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')

	#Add a column to indicate which transcripts were used
	AllTrans_trait_M$Transcripts <- "All"
	#weighted_trait_M$Transcripts <- "Weighted"
	pvalue0.5_trait_M$Transcripts <- "P<0.5"
	pvalue0.1_trait_M$Transcripts <- "P<0.1"
	pvalue0.01_trait_M$Transcripts <- "P<0.01"
	pvalue0.001_trait_M$Transcripts <- "P<0.001"
	pvalue1e_4_trait_M$Transcripts <- "P<0.0001"
	pvalue1e_5_trait_M$Transcripts <- "P<0.00001"
	pvalue1e_6_trait_M$Transcripts<- "P<0.000001"

	#Add a column to indicate which trait
	AllTrans_trait_M$Trait <- trait
	#weighted_trait_M$Trait <- trait
	pvalue0.5_trait_M$Trait <- trait
	pvalue0.1_trait_M$Trait <- trait
	pvalue0.01_trait_M$Trait <- trait
	pvalue0.001_trait_M$Trait <- trait
	pvalue1e_4_trait_M$Trait <- trait
	pvalue1e_5_trait_M$Trait <- trait
	pvalue1e_6_trait_M$Trait<- trait

	#Bind the data
	#data_trait_M <- rbind(AllTrans_trait_M, weighted_trait_M, pvalue0.5_trait_M, pvalue0.1_trait_M, pvalue0.01_trait_M, pvalue0.001_trait_M, pvalue1e_4_trait_M, pvalue1e_5_trait_M, pvalue1e_6_trait_M)
	data_trait_M <- rbind(AllTrans_trait_M, pvalue0.5_trait_M, pvalue0.1_trait_M, pvalue0.01_trait_M, pvalue0.001_trait_M, pvalue1e_4_trait_M, pvalue1e_5_trait_M, pvalue1e_6_trait_M)
	data_trait_M$Sex <- "Males"


	###Create data for final plot
	dat <- rbind(data_trait_F, data_trait_M)


	###Plot the data
	limits <- aes(ymax = dat$mean.r + dat$se.r,
				  ymin = dat$mean.r - dat$se.r)

	p <- ggplot(data = dat, aes(x = Trait, y = mean.r,
								 #fill = factor(Transcripts, levels=c("All", "Weighted", "P<0.5", "P<0.1", "P<0.01", "P<0.001", 
								 fill = factor(Transcripts, levels=c("All", "P<0.5", "P<0.1", "P<0.01", "P<0.001",
																	 "P<0.0001", "P<0.00001", "P<0.000001"),
												#labels=c("All", "Weighted", "P<0.5", "P<0.1", "P<1e-2", "P<1e-3",
												labels=c("All", "P<0.5", "P<0.1", "P<1e-2", "P<1e-3",
														 "P<1e-4", "P<1e-5", "P<1e-6"))))


	gPred <- p + geom_bar(stat = "identity",
						  position = position_dodge(0.9)) +
	  geom_errorbar(limits, position = position_dodge(0.9), width = 0.25) +
	  labs(y = expression(paste("mean", sep=" ", italic(r)))) +
	  facet_grid(. ~ Sex, scales="free_y") +
	  theme_bw()+
	  theme(axis.text=element_text(size=16), axis.title=element_text(size=18), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
			axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")), legend.text=element_text(size=16), legend.title=element_text(size=18), strip.text = element_text(size=16), legend.key.size = unit(2, 'lines')) +
	  scale_fill_brewer(name = "Transcripts", palette="Set1")

	##Save to a pdf
	pdf(paste('Results/Graphs/Manuscript/Accuracy_Expression_TWAS_LMM_adjPheno_common_', trait, '_MS.pdf', sep=""), width=15, height=11)
	print(gPred)
	dev.off()

	rm(list=ls())
}









