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
##traitation
SNPs_trait_F <- read.table(paste("Results/Pred_SNPs_regress_adjPheno_", trait, "_females.csv", sep=""), header=T, sep=',')
Expression_trait_F <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
colnames(Expression_trait_F) <- c('Method', 'mean.Vg', 'se.Vg', 'mean.Ve', 'se.Ve', 'mean.r', 'se.r', 'mean.r2', 'se.r2', 'mean.h2g', 'se.h2g', "mean.e2", "se.e2", "nreps")


####Bind the genomic and transcriptomic data
data_trait_F <- rbind(SNPs_trait_F, Expression_trait_F)

Trait_F <- rep(trait, times=nrow(data_trait_F))
Sex_F <- rep('Females', times=nrow(data_trait_F))
data_F <- cbind(data_trait_F, Trait_F, Sex_F)
colnames(data_F)[15:16] <- c("Trait", "Sex")


###Males
##traitation
SNPs_trait_M <- read.table(paste("Results/Pred_SNPs_regress_adjPheno_", trait, "_males.csv", sep=""), header=T, sep=',')
Expression_trait_M <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
colnames(Expression_trait_M) <- c('Method', 'mean.Vg', 'se.Vg', 'mean.Ve', 'se.Ve', 'mean.r', 'se.r', 'mean.r2', 'se.r2', 'mean.h2g', 'se.h2g', "mean.e2", "se.e2", "nreps")


####Bind the genomic and transcriptomic data
data_trait_M <- rbind(SNPs_trait_M, Expression_trait_M)

Trait_M <- rep(trait, times=nrow(data_trait_M))
Sex_M <- rep('Males', times=nrow(data_trait_M))
data_M <- cbind(data_trait_M, Trait_M, Sex_M)
colnames(data_M)[15:16] <- c("Trait", "Sex")

####Build the final data for plotting
data <- rbind(data_F, data_M)



####Plot the data
##Correlation
limits <- aes(ymax = data$mean.r + data$se.r,
              ymin = data$mean.r - data$se.r)

p <- ggplot(data = data, aes(x = Trait, y = mean.r,
                                fill = Method))


gPred <- p + geom_bar(stat = "identity",
                               position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9), width = 0.25) +
  labs(x = "Trait", y = expression(paste("mean", sep=" ", italic(r)))) +
  facet_grid(. ~ Sex) +
  theme(plot.title=element_text(size=20, hjust=0.5), axis.text=element_text(size=16), axis.title=element_text(size=18), axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), 
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")), legend.text=element_text(size=16), legend.title=element_text(size=18), strip.text = element_text(size=16)) +
  scale_fill_brewer(name = "Method",palette="Set1")

##Save to a pdf
pdf(paste('Results/Graphs/Accuracy_SNPvsExpression_different_kernels_adjPheno_common_', trait,'.pdf', sep=""), width=15, height=11)
print(gPred)
dev.off()



#####GRM+TRM+IRM######

####Load the data
###Females
##Trait
SNPs_trait_F <- read.table(paste("Results/Pred_SNPs_regress_adjPheno_", trait, "_females.csv", sep=""), header=T, sep=',')
SNPs_trait_F <- SNPs_trait_F[, c(1, 6:9)]
Expression_trait_F <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
Expression_trait_F <- Expression_trait_F[, c(1, 6:9)]
SNPs_Expression_trait_F <- read.table(paste("Results/Pred_SNPs+Expression_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
SNPs_Expression_trait_F <- SNPs_Expression_trait_F[, c(1, 8:11)]
SNPs_Expression_Inter_trait_F <- read.table(paste("Results/Pred_SNPs+Expression+Interaction_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
SNPs_Expression_Inter_trait_F <- SNPs_Expression_Inter_trait_F[, c(1, 10:13)]


####Bind the genomic and transcriptomic data
data_trait_F <- rbind(SNPs_trait_F, Expression_trait_F, SNPs_Expression_trait_F, SNPs_Expression_Inter_trait_F)
Trait_F <- rep(trait, times=nrow(data_trait_F))
Sex_F <- rep('Females', times=nrow(data_trait_F))
data_F <- cbind(data_trait_F, Trait_F, Sex_F)
colnames(data_F)[6:7] <- c("Trait", "Sex")


####Select only GRM, TRM and GRM+TRM
data_F <- data_F[which(data_F$Method=='GRM' | data_F$Method=='TRM' | data_F$Method=='GRM+TRM'| data_F$Method=='GRM+TRM+IRM'), ]



###Males
##Trait
SNPs_trait_M <- read.table(paste("Results/Pred_SNPs_regress_adjPheno_", trait, "_males.csv", sep=""), header=T, sep=',')
SNPs_trait_M <- SNPs_trait_M[, c(1, 6:9)]
Expression_trait_M <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
Expression_trait_M <- Expression_trait_M[, c(1, 6:9)]
SNPs_Expression_trait_M <- read.table(paste("Results/Pred_SNPs+Expression_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
SNPs_Expression_trait_M <- SNPs_Expression_trait_M[, c(1, 8:11)]
SNPs_Expression_Inter_trait_M <- read.table(paste("Results/Pred_SNPs+Expression+Interaction_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
SNPs_Expression_Inter_trait_M <- SNPs_Expression_Inter_trait_M[, c(1, 10:13)]


####Bind the genomic and transcriptomic data
data_trait_M <- rbind(SNPs_trait_M, Expression_trait_M, SNPs_Expression_trait_M, SNPs_Expression_Inter_trait_M)
Trait_M <- rep(trait, times=nrow(data_trait_M))
Sex_M <- rep('Males', times=nrow(data_trait_M))
data_M <- cbind(data_trait_M, Trait_M, Sex_M)
colnames(data_M)[6:7] <- c("Trait", "Sex")


####Select only GRM, TRM and GRM+TRM
data_M <- data_M[which(data_M$Method=='GRM' | data_M$Method=='TRM' | data_M$Method=='GRM+TRM'| data_M$Method=='GRM+TRM+IRM'), ]



####Build the final data for plotting
data <- rbind(data_F, data_M)




####Plot the data
##Correlation
limits <- aes(ymax = data$mean.r + data$se.r,
              ymin = data$mean.r - data$se.r)

p <- ggplot(data = data, aes(x = Trait, y = mean.r,
                             fill = factor(Method, levels=c("GRM", "TRM", "GRM+TRM", "GRM+TRM+IRM"))))


gPred <- p + geom_bar(stat = "identity",
                      position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9), width = 0.25) +
  labs(x = "Trait", y = expression(paste("mean", sep=" ", italic(r)))) +
  facet_grid(. ~ Sex) +
  theme(plot.title=element_text(size=20, hjust=0.5), axis.text=element_text(size=16), axis.title=element_text(size=18), axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), 
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")), legend.text=element_text(size=16), legend.title=element_text(size=18), strip.text = element_text(size=16)) +
  scale_fill_brewer(name = "Method",palette="Set1")

##Save to a pdf
pdf(paste('Results/Graphs/Accuracy_GRM_TRM_IRM_adjPheno_common_', trait,'.pdf', sep=""), width=15, height=11)
print(gPred)
dev.off()



