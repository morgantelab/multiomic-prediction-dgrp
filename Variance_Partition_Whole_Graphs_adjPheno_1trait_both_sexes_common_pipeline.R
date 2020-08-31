setMKLthreads(2)

###Set global options
options(stringsAsFactors=FALSE)

###Load libraries
library(Rmisc)
library(ggplot2)

###Select trait to analyze from command line arguments
args <- commandArgs(TRUE)
traitNum <- as.numeric(args)

#Load phenotypes (needed only to get trait name!)
pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)

trait <- colnames(pheno)[traitNum]

####Load the data
###Females
data_trait_F <- read.table(paste("Results/Variance_Partition_SNPs_Expression_Interaction_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
Trait_F <- rep(trait, times=nrow(data_trait_F))
Sex_F <- rep('Females', times=nrow(data_trait_F))
data_F <- cbind(data_trait_F, Trait_F, Sex_F)
colnames(data_F)[6:7] <- c("Trait", "Sex")


###Males
data_trait_M <- read.table(paste("Results/Variance_Partition_SNPs_Expression_Interaction_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
Trait_M <- rep(trait, times=nrow(data_trait_M))
Sex_M <- rep('Males', times=nrow(data_trait_M))
data_M <- cbind(data_trait_M, Trait_M, Sex_M)
colnames(data_M)[6:7] <- c("Trait", "Sex")


####Build the final data for plotting
data <- rbind(data_F, data_M)



####Select only GRM and TRM and select only the columns of interest
data$Proportion_Vg <- data$Vg/rowSums(cbind(data$Vg, data$Vt, data$Vi, data$Ve), na.rm=T)
data$Proportion_Vt <- data$Vt/rowSums(cbind(data$Vg, data$Vt, data$Vi, data$Ve), na.rm=T)
data$Proportion_Vi <- data$Vi/rowSums(cbind(data$Vg, data$Vt, data$Vi, data$Ve), na.rm=T)
data$Proportion_Ve <- data$Ve/rowSums(cbind(data$Vg, data$Vt, data$Vi, data$Ve), na.rm=T)
data_g <- data[, c(1, 2, 6, 7, 8)]
data_g$Partition <- "Vg"
data_t <- data[, c(1, 3, 6, 7, 9)]
data_t$Partition <- "Vt"
data_i <- data[, c(1, 4, 6, 7, 10)]
data_i$Partition <- "Vi"
data_e <- data[, c(1, 5, 6, 7, 11)]
data_e$Partition <- "Ve"
colnames(data_g)[c(2, 5)] <- c("Variance", "Proportion_Variance")
colnames(data_t)[c(2, 5)] <- c("Variance", "Proportion_Variance")
colnames(data_i)[c(2, 5)] <- c("Variance", "Proportion_Variance")
colnames(data_e)[c(2, 5)] <- c("Variance", "Proportion_Variance")

data <- rbind(data_g, data_t, data_i, data_e)




####Plot the data
###Variance partition
#Proportion of variance
pPropVarPart <- ggplot(data = data, aes(x = factor(Method, levels=c("GRM", "TRM", "GRM+TRM", "GRM+TRM+IRM")), y = Proportion_Variance,
                                                          fill = Partition))

gPropVarPart <- pPropVarPart + geom_bar(stat = "identity") + facet_grid(. ~ Sex) + ggtitle(trait) +
  labs(x = "Method", y = "Proportion of variance explained") +
  theme(plot.title=element_text(size=20, hjust=0.5), axis.text=element_text(size=16), axis.title=element_text(size=18), axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), 
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")), legend.text=element_text(size=16), legend.title=element_text(size=18), strip.text = element_text(size=16), 
        axis.text.x = element_text(angle = 90, hjust = 1))

pdf(paste('Results/Graphs/Proportion_Variance_Partition_WholeData_GRM_TRM_IRM_adjPheno_common_', trait, '.pdf', sep=""), width=15, height=11)
print(gPropVarPart)
dev.off()



rm(list=ls())



