########################
## Fabio Morgante
## 6/1/2020
## Plot variance partition
## in the training set
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

	###Females
	##Load and reshape the data
	SNPs_trait_F <- read.table(paste("Results/Pred_SNPs_regress_adjPheno_", trait, "_females.csv", sep=""), header=T, sep=',')
	SNPs_trait_F <- SNPs_trait_F[, c(1, 10:13)]
	SNPs_trait_F$se.h2i <- SNPs_trait_F$mean.h2i <- SNPs_trait_F$se.h2t <- SNPs_trait_F$mean.h2t <- NA

	Expression_trait_F <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	Expression_trait_F <- Expression_trait_F[, c(1, 10:13)]
	Expression_trait_F$se.h2i <- Expression_trait_F$mean.h2i <- Expression_trait_F$se.h2g <- Expression_trait_F$mean.h2g <- NA

	SNPs_Expression_trait_F <- read.table(paste("Results/Pred_SNPs+Expression_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	SNPs_Expression_trait_F <- SNPs_Expression_trait_F[, c(1, 12:17)]
	SNPs_Expression_trait_F$se.h2i <- SNPs_Expression_trait_F$mean.h2i <- NA

	SNPs_Expression_Inter_trait_F <- read.table(paste("Results/Pred_SNPs+Expression+Interaction_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	SNPs_Expression_Inter_trait_F <- SNPs_Expression_Inter_trait_F[, c(1, 14:21)]

	data_trait_F <- rbind(SNPs_trait_F, Expression_trait_F, SNPs_Expression_trait_F, SNPs_Expression_Inter_trait_F)

	Trait_F <- rep(trait, times=nrow(data_trait_F))
	Sex_F <- rep('Females', times=nrow(data_trait_F))
	data_F <- cbind(data_trait_F, Trait_F, Sex_F)
	colnames(data_F)[10:11] <- c("Trait", "Sex")

	##Select only GRM, TRM and GRM+TRM
	data_F <- data_F[which(data_F$Method=='GRM' | data_F$Method=='TRM' | data_F$Method=='GRM+TRM'| data_F$Method=='GRM+TRM+IRM'), ]


	###Males
	##Load and reshape the data
	SNPs_trait_M <- read.table(paste("Results/Pred_SNPs_regress_adjPheno_", trait, "_males.csv", sep=""), header=T, sep=',')
	SNPs_trait_M <- SNPs_trait_M[, c(1, 10:13)]
	SNPs_trait_M$se.h2i <- SNPs_trait_M$mean.h2i <- SNPs_trait_M$se.h2t <- SNPs_trait_M$mean.h2t <- NA

	Expression_trait_M <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	Expression_trait_M <- Expression_trait_M[, c(1, 10:13)]
	Expression_trait_M$se.h2i <- Expression_trait_M$mean.h2i <- Expression_trait_M$se.h2g <- Expression_trait_M$mean.h2g <- NA

	SNPs_Expression_trait_M <- read.table(paste("Results/Pred_SNPs+Expression_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	SNPs_Expression_trait_M <- SNPs_Expression_trait_M[, c(1, 12:17)]
	SNPs_Expression_trait_M$se.h2i <- SNPs_Expression_trait_M$mean.h2i <- NA

	SNPs_Expression_Inter_trait_M <- read.table(paste("Results/Pred_SNPs+Expression+Interaction_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	SNPs_Expression_Inter_trait_M <- SNPs_Expression_Inter_trait_M[, c(1, 14:21)]

	data_trait_M <- rbind(SNPs_trait_M, Expression_trait_M, SNPs_Expression_trait_M, SNPs_Expression_Inter_trait_M)

	Trait_M <- rep(trait, times=nrow(data_trait_M))
	Sex_M <- rep('Males', times=nrow(data_trait_M))
	data_M <- cbind(data_trait_M, Trait_M, Sex_M)
	colnames(data_M)[10:11] <- c("Trait", "Sex")

	##Select only GRM, TRM and GRM+TRM
	data_M <- data_M[which(data_M$Method=='GRM' | data_M$Method=='TRM' | data_M$Method=='GRM+TRM'| data_M$Method=='GRM+TRM+IRM'), ]


	###Build the final data for plotting
	dat <- rbind(data_F, data_M)

	###Select only the columns of interest
	data_g <- dat[, c(1, 2, 3, 10, 11)]
	data_g$Partition <- "Vg"
	data_t <- dat[, c(1, 6, 7, 10, 11)]
	data_t$Partition <- "Vt"
	data_i <- dat[, c(1, 8, 9, 10, 11)]
	data_i$Partition <- "Vi"
	data_e <- dat[, c(1, 4, 5, 10, 11)]
	data_e$Partition <- "Ve"
	colnames(data_g)[c(2, 3)] <- c("Proportion_Variance", "se_Proportion_Variance")
	colnames(data_t)[c(2, 3)] <- c("Proportion_Variance", "se_Proportion_Variance")
	colnames(data_i)[c(2, 3)] <- c("Proportion_Variance", "se_Proportion_Variance")
	colnames(data_e)[c(2, 3)] <- c("Proportion_Variance", "se_Proportion_Variance")

	dat <- rbind(data_g, data_t, data_i, data_e)

	###Calculate lower and upper bounds for se
	dat$lower <- dat$Proportion_Variance - dat$se_Proportion_Variance
	dat$upper <- dat$Proportion_Variance + dat$se_Proportion_Variance 

	dat[which(dat$Method=='GRM' & dat$Partition=='Ve'), "lower"] <- dat[which(dat$Method=='GRM' & dat$Partition=='Vg'), "Proportion_Variance"] + dat[which(dat$Method=='GRM' & dat$Partition=='Ve'), "lower"]
	dat[which(dat$Method=='GRM' & dat$Partition=='Ve'), "upper"] <- dat[which(dat$Method=='GRM' & dat$Partition=='Vg'), "Proportion_Variance"] + dat[which(dat$Method=='GRM' & dat$Partition=='Ve'), "upper"]

	dat[which(dat$Method=='TRM' & dat$Partition=='Ve'), "lower"] <- dat[which(dat$Method=='TRM' & dat$Partition=='Vt'), "Proportion_Variance"] + dat[which(dat$Method=='TRM' & dat$Partition=='Ve'), "lower"]
	dat[which(dat$Method=='TRM' & dat$Partition=='Ve'), "upper"] <- dat[which(dat$Method=='TRM' & dat$Partition=='Vt'), "Proportion_Variance"] + dat[which(dat$Method=='TRM' & dat$Partition=='Ve'), "upper"]

	dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Vt'), "lower"] <- dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Vg'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Vt'), "lower"]
	dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Vt'), "upper"] <- dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Vg'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Vt'), "upper"]
	dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Ve'), "lower"] <- dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Vg'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Vt'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Ve'), "lower"]
	dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Ve'), "upper"] <- dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Vg'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Vt'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM' & dat$Partition=='Ve'), "upper"]

	dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vt'), "lower"] <- dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vg'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vt'), "lower"]
	dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vt'), "upper"] <- dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vg'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vt'), "upper"]
	dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vi'), "lower"] <- dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vg'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vt'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vi'), "lower"]
	dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vi'), "upper"] <- dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vg'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vt'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vi'), "upper"]
	dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Ve'), "lower"] <- dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vg'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vt'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vi'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Ve'), "lower"]
	dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Ve'), "upper"] <- dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vg'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vt'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Vi'), "Proportion_Variance"] + dat[which(dat$Method=='GRM+TRM+IRM' & dat$Partition=='Ve'), "upper"]


	###Plot the data
	pPropVarPart_starv <- ggplot(data = dat, aes(x = factor(Method, levels=c("GRM", "TRM", "GRM+TRM", "GRM+TRM+IRM"), labels=c("GBLUP", "TBLUP", "GTBLUP", "GTIBLUP")), y = Proportion_Variance,
															  fill = factor(Partition, levels=c("Ve", "Vi", "Vt", "Vg"))))

	gPropVarPart_starv <- pPropVarPart_starv + geom_bar(stat = "identity") + 
	  geom_errorbar(aes(ymax = upper, ymin = lower, width = 0.25), position = "identity") + 
	  facet_grid(. ~ Sex) + 
	  scale_fill_discrete(labels=(expression(italic(sigma[e]^2), italic(sigma[i]^2), italic(sigma[t]^2), italic(sigma[g]^2)))) +
	  labs(x = "Method", y = "Mean proportion of variance explained", fill="Partition") +
	  theme_bw() +
	  theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1),
			axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")), legend.text=element_text(size=10), legend.title=element_text(size=12), strip.text = element_text(size=10), legend.key.size = unit(2, 'lines'))
	
	##Save to a pdf
	jpeg(paste('Results/Graphs/Manuscript/Proportion_Variance_Partition_TrainingSet_GRM_TRM_IRM_adjPheno_common_', trait, '_MS.jpeg', sep=""), width=20, height=15, units="cm", res=300)
	print(gPropVarPart_starv)
	dev.off()
}


rm(list=ls())



