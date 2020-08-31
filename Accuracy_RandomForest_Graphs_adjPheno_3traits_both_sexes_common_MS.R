########################
## Fabio Morgante
## 6/1/2020
## Plot accuracy for TBLUP vs RF  
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
	Expression_trait_F <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	Expression_trait_F <- Expression_trait_F[1, c(1, 6:9)]
	RF_Expression_trait_F <- read.table(paste("Results/Pred_Expression_RandomForest_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
	RF_Expression_trait_F <- RF_Expression_trait_F[, c(1, 3:6)]

	#Bind the genomic and transcriptomic data
	data_trait_F <- rbind(Expression_trait_F, RF_Expression_trait_F)
	Trait_F <- rep(trait, times=nrow(data_trait_F))
	Sex_F <- rep('Females', times=nrow(data_trait_F))
	data_F <- cbind(data_trait_F, Trait_F, Sex_F)
	colnames(data_F)[6:7] <- c("Trait", "Sex")

	##Males
	Expression_trait_M <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	Expression_trait_M <- Expression_trait_M[1, c(1, 6:9)]
	RF_Expression_trait_M <- read.table(paste("Results/Pred_Expression_RandomForest_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
	RF_Expression_trait_M <- RF_Expression_trait_M[, c(1, 3:6)]

	#Bind the genomic and transcriptomic data
	data_trait_M <- rbind(Expression_trait_M, RF_Expression_trait_M)
	Trait_M <- rep(trait, times=nrow(data_trait_M))
	Sex_M <- rep('Males', times=nrow(data_trait_M))
	data_M <- cbind(data_trait_M, Trait_M, Sex_M)
	colnames(data_M)[6:7] <- c("Trait", "Sex")


	###Build the final data for plotting
	dat <- rbind(data_F, data_M)


	###Plot the data
	limits <- aes(ymax = dat$mean.r + dat$se.r,
				  ymin = dat$mean.r - dat$se.r)

	p <- ggplot(data = dat, aes(x = Trait, y = mean.r,
								 fill = factor(Method, levels=c("TRM", "TRF"), labels=c("TBLUP", "Random Forest"))))


	gPred <- p + geom_bar(stat = "identity",
						  position = position_dodge(0.9)) +
	  geom_errorbar(limits, position = position_dodge(0.9), width = 0.25) +
	  facet_grid(. ~ Sex, scales="free_y") +
	  theme_bw()+
	  theme(axis.text=element_text(size=16), axis.title=element_text(size=18), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
			axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")), legend.text=element_text(size=16), legend.title=element_text(size=18), strip.text = element_text(size=16), legend.key.size = unit(2, 'lines')) +
	  scale_fill_brewer(name = "Method",palette="Set1")

	##Save to a pdf
	pdf(paste('Results/Graphs/Manuscript/Accuracy_TRMvsRandomForest_adj_Pheno_common_', trait,'_MS.pdf', sep=""), width=15, height=11)
	print(gPred)
	dev.off()

	rm(list=ls())
}


