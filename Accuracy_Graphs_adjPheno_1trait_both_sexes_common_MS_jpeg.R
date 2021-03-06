########################
## Fabio Morgante
## 6/1/2020
## Plot accuracy for [G,T,I]BLUP   
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
	SNPs_trait_F <- read.table(paste("Results/Pred_SNPs_regress_adjPheno_", trait, "_females.csv", sep=""), header=T, sep=',')
	SNPs_trait_F <- SNPs_trait_F[, c(1, 6:9)]
	Expression_trait_F <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	Expression_trait_F <- Expression_trait_F[, c(1, 6:9)]
	SNPs_Expression_trait_F <- read.table(paste("Results/Pred_SNPs+Expression_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	SNPs_Expression_trait_F <- SNPs_Expression_trait_F[, c(1, 8:11)]
	SNPs_Expression_Inter_trait_F <- read.table(paste("Results/Pred_SNPs+Expression+Interaction_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	SNPs_Expression_Inter_trait_F <- SNPs_Expression_Inter_trait_F[, c(1, 10:13)]

	#Bind the genomic and transcriptomic data
	data_trait_F <- rbind(SNPs_trait_F, Expression_trait_F, SNPs_Expression_trait_F, SNPs_Expression_Inter_trait_F)
	Trait_F <- rep(trait, times=nrow(data_trait_F))
	Sex_F <- rep('Females', times=nrow(data_trait_F))
	data_F <- cbind(data_trait_F, Trait_F, Sex_F)
	colnames(data_F)[6:7] <- c("Trait", "Sex")

	#Select only GRM, TRM and GRM+TRM
	data_F <- data_F[which(data_F$Method=='GRM' | data_F$Method=='TRM' | data_F$Method=='GRM+TRM'| data_F$Method=='GRM+TRM+IRM'), ]

	##Males
	SNPs_trait_M <- read.table(paste("Results/Pred_SNPs_regress_adjPheno_", trait, "_males.csv", sep=""), header=T, sep=',')
	SNPs_trait_M <- SNPs_trait_M[, c(1, 6:9)]
	Expression_trait_M <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	Expression_trait_M <- Expression_trait_M[, c(1, 6:9)]
	SNPs_Expression_trait_M <- read.table(paste("Results/Pred_SNPs+Expression_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	SNPs_Expression_trait_M <- SNPs_Expression_trait_M[, c(1, 8:11)]
	SNPs_Expression_Inter_trait_M <- read.table(paste("Results/Pred_SNPs+Expression+Interaction_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	SNPs_Expression_Inter_trait_M <- SNPs_Expression_Inter_trait_M[, c(1, 10:13)]

	#Bind the genomic and transcriptomic data
	data_trait_M <- rbind(SNPs_trait_M, Expression_trait_M, SNPs_Expression_trait_M, SNPs_Expression_Inter_trait_M)
	Trait_M <- rep(trait, times=nrow(data_trait_M))
	Sex_M <- rep('Males', times=nrow(data_trait_M))
	data_M <- cbind(data_trait_M, Trait_M, Sex_M)
	colnames(data_M)[6:7] <- c("Trait", "Sex")

	#Select only GRM, TRM and GRM+TRM
	data_M <- data_M[which(data_M$Method=='GRM' | data_M$Method=='TRM' | data_M$Method=='GRM+TRM'| data_M$Method=='GRM+TRM+IRM'), ]

	###Build the final data for plotting
	dat <- rbind(data_F, data_M)


	###Plot the data
	limits <- aes(ymax = dat$mean.r + dat$se.r,
				  ymin = dat$mean.r - dat$se.r)

	p <- ggplot(data = dat, aes(x = Trait, y = mean.r,
								 fill = factor(Method, levels=c("GRM", "TRM", "GRM+TRM", "GRM+TRM+IRM"), labels=c("GBLUP", "TBLUP", "GTBLUP", "GTIBLUP"))))


	gPred <- p + geom_bar(stat = "identity",
						  position = position_dodge(0.9)) +
	  geom_errorbar(limits, position = position_dodge(0.9), width = 0.25) +
	  labs(y = expression(paste("mean", sep=" ", italic(r)))) +
	  facet_grid(. ~ Sex, scales="free_y") +
	  theme_bw()+
	  theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
			axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")), legend.text=element_text(size=10), legend.title=element_text(size=12), strip.text = element_text(size=10), legend.key.size = unit(2, 'lines')) +
	  scale_fill_brewer(name = "Method",palette="Set1")

	##Save to a pdf
	jpeg(paste('Results/Graphs/Manuscript/Accuracy_GRM_TRM_IRM_adjPheno_common_', trait,'_MS.jpeg', sep=""), width=20, height=15, units="cm", res=300)
	print(gPred)
	dev.off()


	rm(list=ls())
}
