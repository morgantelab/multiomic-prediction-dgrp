setMKLthreads(2)

###Set global options
options(stringsAsFactors=FALSE)

###Select trait to analyze from command line arguments
args <- commandArgs(TRUE)
traitNum <- as.numeric(args)

#Load phenotypes (needed only to get trait name!)
pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)

trait <- colnames(pheno)[traitNum]

###Load libraries
library(ggplot2)
library(Rmisc)

####Load the data
###Females
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

####Bind the genomic and transcriptomic data
data_trait_F <- rbind(SNPs_trait_F, Expression_trait_F, SNPs_Expression_trait_F, SNPs_Expression_Inter_trait_F)

Trait_F <- rep(trait, times=nrow(data_trait_F))
Sex_F <- rep('Females', times=nrow(data_trait_F))
data_F <- cbind(data_trait_F, Trait_F, Sex_F)
colnames(data_F)[10:11] <- c("Trait", "Sex")


####Select only GRM, TRM and GRM+TRM
data_F <- data_F[which(data_F$Method=='GRM' | data_F$Method=='TRM' | data_F$Method=='GRM+TRM'| data_F$Method=='GRM+TRM+IRM'), ]



###Males
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

####Bind the genomic and transcriptomic data
data_trait_M <- rbind(SNPs_trait_M, Expression_trait_M, SNPs_Expression_trait_M, SNPs_Expression_Inter_trait_M)

Trait_M <- rep(trait, times=nrow(data_trait_M))
Sex_M <- rep('Males', times=nrow(data_trait_M))
data_M <- cbind(data_trait_M, Trait_M, Sex_M)
colnames(data_M)[10:11] <- c("Trait", "Sex")


####Select only GRM, TRM and GRM+TRM
data_M <- data_M[which(data_M$Method=='GRM' | data_M$Method=='TRM' | data_M$Method=='GRM+TRM'| data_M$Method=='GRM+TRM+IRM'), ]



####Build the final data for plotting
data <- rbind(data_F, data_M)


####Select only GRM and TRM and select only the columns of interest
data_g <- data[, c(1, 2, 3, 10, 11)]
data_g$Partition <- "Vg"
data_t <- data[, c(1, 6, 7, 10, 11)]
data_t$Partition <- "Vt"
data_i <- data[, c(1, 8, 9, 10, 11)]
data_i$Partition <- "Vi"
data_e <- data[, c(1, 4, 5, 10, 11)]
data_e$Partition <- "Ve"
colnames(data_g)[c(2, 3)] <- c("Proportion_Variance", "se_Proportion_Variance")
colnames(data_t)[c(2, 3)] <- c("Proportion_Variance", "se_Proportion_Variance")
colnames(data_i)[c(2, 3)] <- c("Proportion_Variance", "se_Proportion_Variance")
colnames(data_e)[c(2, 3)] <- c("Proportion_Variance", "se_Proportion_Variance")

data <- rbind(data_g, data_t, data_i, data_e)

###Calculate lower and upper bounds for se
data$lower <- data$Proportion_Variance - data$se_Proportion_Variance
data$upper <- data$Proportion_Variance + data$se_Proportion_Variance 

data[which(data$Method=='GRM' & data$Partition=='Ve'), "lower"] <- data[which(data$Method=='GRM' & data$Partition=='Vg'), "Proportion_Variance"] + data[which(data$Method=='GRM' & data$Partition=='Ve'), "lower"]
data[which(data$Method=='GRM' & data$Partition=='Ve'), "upper"] <- data[which(data$Method=='GRM' & data$Partition=='Vg'), "Proportion_Variance"] + data[which(data$Method=='GRM' & data$Partition=='Ve'), "upper"]

data[which(data$Method=='TRM' & data$Partition=='Ve'), "lower"] <- data[which(data$Method=='TRM' & data$Partition=='Vt'), "Proportion_Variance"] + data[which(data$Method=='TRM' & data$Partition=='Ve'), "lower"]
data[which(data$Method=='TRM' & data$Partition=='Ve'), "upper"] <- data[which(data$Method=='TRM' & data$Partition=='Vt'), "Proportion_Variance"] + data[which(data$Method=='TRM' & data$Partition=='Ve'), "upper"]

data[which(data$Method=='GRM+TRM' & data$Partition=='Vt'), "lower"] <- data[which(data$Method=='GRM+TRM' & data$Partition=='Vg'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM' & data$Partition=='Vt'), "lower"]
data[which(data$Method=='GRM+TRM' & data$Partition=='Vt'), "upper"] <- data[which(data$Method=='GRM+TRM' & data$Partition=='Vg'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM' & data$Partition=='Vt'), "upper"]
data[which(data$Method=='GRM+TRM' & data$Partition=='Ve'), "lower"] <- data[which(data$Method=='GRM+TRM' & data$Partition=='Vg'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM' & data$Partition=='Vt'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM' & data$Partition=='Ve'), "lower"]
data[which(data$Method=='GRM+TRM' & data$Partition=='Ve'), "upper"] <- data[which(data$Method=='GRM+TRM' & data$Partition=='Vg'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM' & data$Partition=='Vt'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM' & data$Partition=='Ve'), "upper"]

data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vt'), "lower"] <- data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vg'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vt'), "lower"]
data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vt'), "upper"] <- data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vg'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vt'), "upper"]
data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vi'), "lower"] <- data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vg'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vt'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vi'), "lower"]
data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vi'), "upper"] <- data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vg'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vt'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vi'), "upper"]
data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Ve'), "lower"] <- data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vg'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vt'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vi'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Ve'), "lower"]
data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Ve'), "upper"] <- data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vg'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vt'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Vi'), "Proportion_Variance"] + data[which(data$Method=='GRM+TRM+IRM' & data$Partition=='Ve'), "upper"]



####Plot the data
###Variance partition
##Starvation
#Proportion of variance
pPropVarPart_starv <- ggplot(data = data, aes(x = factor(Method, levels=c("GRM", "TRM", "GRM+TRM", "GRM+TRM+IRM")), y = Proportion_Variance,
                                                          fill = factor(Partition, levels=c("Ve", "Vi", "Vt", "Vg"))))

gPropVarPart_starv <- pPropVarPart_starv + geom_bar(stat = "identity") + geom_errorbar(aes(ymax = upper, ymin = lower, width = 0.25), position = "identity") + facet_grid(. ~ Sex) + ggtitle(trait) + scale_fill_discrete(name="Partition") +
  labs(x = "Method", y = "Mean proportion of variance explained") +
  theme(plot.title=element_text(size=20, hjust=0.5), axis.text=element_text(size=16), axis.title=element_text(size=18), axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), 
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")), legend.text=element_text(size=16), legend.title=element_text(size=18), strip.text = element_text(size=16), 
        axis.text.x = element_text(angle = 90, hjust = 1))

pdf(paste('Results/Graphs/Proportion_Variance_Partition_TrainingSet_GRM_TRM_IRM_adjPheno_common_', trait, '.pdf', sep=""), width=15, height=11)
print(gPropVarPart_starv)
dev.off()



rm(list=ls())



