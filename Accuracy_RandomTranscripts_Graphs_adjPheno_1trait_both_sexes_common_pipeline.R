setMKLthreads(2)

###Set global options
options(stringsAsFactors=FALSE)


library(ggplot2)
library(Rmisc)


###Select trait to analyze from command line arguments
args <- commandArgs(TRUE)
traitNum <- as.numeric(args)

#Load phenotypes (needed only to get trait name!)
pheno <- read.csv("Data/eQTL_traits_females.csv", header=TRUE)

trait <- colnames(pheno)[traitNum]

###Females
##Starvation
#Load the data
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



###Males
##Starvation
#Load the data
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
data <- rbind(data_trait_F, data_trait_M)
##Replace NA with 0. Useful to plot nreps over absent bars
data[is.na(data)] <- 0


####Plot the data
##Correlation
limits <- aes(ymax = data$mean.r + data$se.r,
              ymin = data$mean.r - data$se.r)

p <- ggplot(data = data, aes(x = Trait, y = mean.r,
                             fill = factor(Transcripts, levels=c("All", "Random 5000", "Random 1000", "Random 500", "Random 50", "Random 5"))))


gPred <- p + geom_bar(stat = "identity",
                      position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9), width = 0.25) +
  #geom_text(aes(label=nreps), position=position_dodge(width=0.9), vjust=-5.25)+
  labs(x = "Trait", y = expression(paste("mean", sep=" ", italic(r)))) +
  facet_grid(. ~ Sex) +
  theme(plot.title=element_text(size=20, hjust=0.5), axis.text=element_text(size=16), axis.title=element_text(size=18), axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), 
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")), legend.text=element_text(size=16), legend.title=element_text(size=18), strip.text = element_text(size=16), legend.key.size = unit(2, 'lines')) +
  scale_fill_brewer(name = "Transcripts",palette="Set1")

##Save to a pdf
pdf(paste('Results/Graphs/Accuracy_Expression_RandomTranscripts_adjPheno_common_', trait, '.pdf', sep=""), width=15, height=11)
gPred
dev.off()









