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
##Trait
GO_trait_F <- read.table(paste("Results/Pred_Expression_GO_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
GO_trait_F <- GO_trait_F[, c(1, 8:11)]

####Bind the genomic and transcriptomic data
Trait_F <- rep(trait, times=nrow(GO_trait_F))
Sex_F <- rep('Females', times=nrow(GO_trait_F))
rank_F <- rank(-GO_trait_F$mean.r)
data_F <- cbind(GO_trait_F, rank_F, Trait_F, Sex_F)
colnames(data_F)[6:8] <- c("Rank", "Trait", "Sex")

All_trait_F <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_females.csv", sep=""), header=T, sep=',')
All_trait_F <- data.frame(mean.r=All_trait_F[1, 6])
All_trait_F$Sex <- c('Females')

###Males
##Trait
GO_trait_M <- read.table(paste("Results/Pred_Expression_GO_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
GO_trait_M <- GO_trait_M[, c(1, 8:11)]


####Bind the genomic and transcriptomic data
Trait_M <- rep(trait, times=nrow(GO_trait_M))
Sex_M <- rep('Males', times=nrow(GO_trait_M))
rank_M <- rank(-GO_trait_M$mean.r)
data_M <- cbind(GO_trait_M, rank_M, Trait_M, Sex_M)
colnames(data_M)[6:8] <- c("Rank", "Trait", "Sex")

All_trait_M <- read.table(paste("Results/Pred_Expression_regress_adjPheno_common_", trait, "_males.csv", sep=""), header=T, sep=',')
All_trait_M <- data.frame(mean.r=All_trait_M[1, 6])
All_trait_M$Sex <- c('Males')



####Build the final data for plotting
data <- rbind(data_F, data_M)
All <- rbind(All_trait_F, All_trait_M)



####Plot the data
##Correlation
p <- ggplot(data = data, aes(x=GO, y = mean.r, colour=Sex))

gPred <- p +  geom_point(shape=16, size=1.5) +
  
  labs(x = "GO", y = expression(paste("mean", sep=" ", italic(r)))) +
  scale_x_discrete(breaks=seq(from=1, to=(nrow(data)/2), by=1)) +
  facet_grid(Sex ~ .) +
  geom_text(aes(label=ifelse((Rank==1|Rank==2|Rank==3), as.character(GO), '')), hjust=0.5, vjust=0) +
  geom_hline(data = All, aes(yintercept = mean.r)) +
  ggtitle(paste(trait)) +
  theme(axis.text.y=element_text(size=16), axis.title=element_text(size=18), axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.text = element_text(size=16), legend.position="none", plot.title=element_text(size=20, hjust=0.5))

##Save to a pdf
pdf(paste('Results/Graphs/Accuracy_GO-transcripts_adjPheno_common_', trait,'.pdf', sep=""), width=15, height=11)
print(gPred)
dev.off()
