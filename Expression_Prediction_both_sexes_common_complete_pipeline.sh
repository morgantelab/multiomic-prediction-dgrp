#!/bin/bash


########################
## Fabio Morgante
## 4/6/2017
## Gene expression prediction pipeline
##
## Modified on 4-5-2018
## Reason: use adjusted phenotype
##
## Modified on 11-11-2018
## Reason: Remove steps not to be included in the paper
########################

###Get trait from arguments
traitFem=$1
traitMal=$2


###############
####Females####
###############


###GBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript GBLUP_regress_adjPheno_females.R "$traitFem" > "GBLUP_regress_adjPheno_females_${traitFem}.Rout"
echo Finished GBLUP in females

###TBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript TBLUP_regress_adjPheno_common_females.R "$traitFem" > "TBLUP_regress_adjPheno_common_females_${traitFem}.Rout"
echo Finished TBLUP in females

###GTBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript GTBLUP_regress_adjPheno_common_females.R "$traitFem" > "GTBLUP_regress_adjPheno_common_females_${traitFem}.Rout"
echo Finished GTBLUP in females

###GT-I-BLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript GT_I_BLUP_regress_adjPheno_common_females.R "$traitFem" > "GT_I_BLUP_regress_adjPheno_common_females_${traitFem}.Rout"
echo Finished GTIBLUP in females

###GT-I-REML
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript GT_I_REML_regress_adjPheno_common_females.R "$traitFem" > "GT_I_REML_regress_adjPheno_common_females_${traitFem}.Rout"
echo Finished GTIREML in females

###TBLUP-TWAS
./TBLUP_TWAS_LMM_pvalues_loop_common_females.sh "$traitFem"
echo Finished TWAS-TBLUP in females

###TBLUP_TWAS (weighted TRM)
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript TBLUP_TWAS_LMM_weighted_regress_adjPheno_common_females.R "$traitFem" > "TBLUP_TWAS_LMM_weighted_regress_adjPheno_common_females_${traitFem}.Rout"
echo Finished TWAS-TBLUP weighted in females

###Random transcripts
./TBLUP_random_number_loop_common_females.sh "$traitFem"
echo Finished Random-TBLUP in females

###Random Forest
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript T-RF_adjPheno_common_females.R "$traitFem" 1000 5 > "T-RF_adjPheno_common_females_${traitFem}.Rout"
echo Finished Random Forest in females

###GO-TBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript TBLUP_GO_regress_adjPheno_common_females.R "$traitFem" > "TBLUP_GO_regress_adjPheno_common_females_${traitFem}.Rout"
echo Finished GO-TBLUP in females

###GO-GBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript GBLUP_GO_regress_adjPheno_females.R "$traitFem" > "GBLUP_GO_regress_adjPheno_females_${traitFem}.Rout"
echo Finished GO-GBLUP in females

###GO-GTBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript GTBLUP_GO_regress_adjPheno_common_females.R "$traitFem" > "GTBLUP_GO_regress_adjPheno_common_females_${traitFem}.Rout"
echo Finished GO-GTBLUP in females





#############
####Males####
#############


###GBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript GBLUP_regress_adjPheno_males.R "$traitMal" > "GBLUP_regress_adjPheno_males_${traitMal}.Rout"
echo Finished GBLUP in males

###TBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript TBLUP_regress_adjPheno_common_males.R "$traitMal" > "TBLUP_regress_adjPheno_common_males_${traitMal}.Rout"
echo Finished TBLUP in males

###GTBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript GTBLUP_regress_adjPheno_common_males.R "$traitMal" > "GTBLUP_regress_adjPheno_common_males_${traitMal}.Rout"
echo Finished GTBLUP in males

###GT-I-BLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript GT_I_BLUP_regress_adjPheno_common_males.R "$traitMal" > "GT_I_BLUP_regress_adjPheno_common_males_${traitMal}.Rout"
echo Finished GTIBLUP in males

###GT-I-REML
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript GT_I_REML_regress_adjPheno_common_males.R "$traitMal" > "GT_I_REML_regress_adjPheno_common_males_${traitMal}.Rout"
echo Finished GTIREML in males

###TBLUP-TWAS LMM
./TBLUP_TWAS_LMM_pvalues_loop_common_males.sh "$traitMal"
echo Finished TWAS-TBLUP in males

###TBLUP_TWAS LMM (weighted TRM)
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript TBLUP_TWAS_LMM_weighted_regress_adjPheno_common_males.R "$traitMal" > "TBLUP_TWAS_LMM_weighted_regress_adjPheno_common_males_${traitMal}.Rout"
echo Finished TWAS-TBLUP weighted in males

###Random transcripts
./TBLUP_random_number_loop_common_males.sh "$traitMal"
echo Finished Random-TBLUP in males

###Random Forest
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript T-RF_adjPheno_common_males.R "$traitMal" 1000 5 > "T-RF_adjPheno_common_males_${traitMal}.Rout"
echo Finished Random Forest in males

###GO-TBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript TBLUP_GO_regress_adjPheno_common_males.R "$traitMal" > "TBLUP_GO_regress_adjPheno_common_males_${traitMal}.Rout"
echo Finished GO-TBLUP in males

###GO-GBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript GBLUP_GO_regress_adjPheno_males.R "$traitMal" > "GBLUP_GO_regress_adjPheno_males_${traitMal}.Rout"
echo Finished GO-GBLUP in males

###GO-GTBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript GTBLUP_GO_regress_adjPheno_common_males.R "$traitMal" > "GTBLUP_GO_regress_adjPheno_common_males_${traitMal}.Rout"
echo Finished GO-GTBLUP in males




############################
####Extract useful stuff####
############################

####Females####

###TBLUP - Extract genes in the top 3 most predictive GO terms and calculate the overlap
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript Extract_genes_within_topGOs_TBLUP_GO_mean_r_females.R "$traitFem" > "Extract_genes_within_topGOs_TBLUP_GO_mean_r_females_${traitFem}.Rout"
echo Finished extracting genes in the top 3 most predictive GO terms from TBLUP in females

###GBLUP - Extract genes in the top 3 most predictive GO terms and calculate the overlap
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript Extract_genes_and_variants_within_topGOs_GBLUP_GO_mean_r_females.R "$traitFem" > "Extract_genes_and_variants_within_topGOs_GBLUP_GO_mean_r_females_${traitFem}.Rout"
echo Finished extracting genes in the top 3 most predictive GO terms from GBLUP in females


####Males####

###TBLUP - Extract genes in the top 3 most predictive GO terms and calculate the overlap
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript Extract_genes_within_topGOs_TBLUP_GO_mean_r_males.R "$traitMal" > "Extract_genes_within_topGOs_TBLUP_GO_mean_r_males_${traitMal}.Rout"
echo Finished extracting genes in the top 3 most predictive GO terms from TBLUP in males

###GBLUP - Extract genes in the top 3 most predictive GO terms and calculate the overlap
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript Extract_genes_and_variants_within_topGOs_GBLUP_GO_mean_r_males.R "$traitMal" > "Extract_genes_and_variants_within_topGOs_GBLUP_GO_mean_r_males_${traitMal}.Rout"
echo Finished extracting genes in the top 3 most predictive GO terms from GBLUP in males





#############
####Plots####
#############

###Trait names is retrieved by loading the FEMALE phenotype file in each of the following scripts

###Different kernels + GT and GT_I_BLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript Accuracy_Graphs_adjPheno_1trait_both_sexes_common_pipeline.R "$traitFem" > "Accuracy_Graphs_adjPheno_1trait_both_sexes_common_pipeline_${traitFem}.Rout"
echo Finished plotting BLUPs

###TWAS LMM different pvalues
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript Accuracy_TWAS_LMM_Graphs_adjPheno_1trait_both_sexes_common_pipeline.R "$traitFem" > "Accuracy_TWAS_LMM_Graphs_adjPheno_1trait_both_sexes_common_pipeline_${traitFem}.Rout"
echo Finished plotting TWAS-TBLUP

###Random transcripts
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript Accuracy_RandomTranscripts_Graphs_adjPheno_1trait_both_sexes_common_pipeline.R "$traitFem" > "Accuracy_RandomTranscripts_Graphs_adjPheno_1trait_both_sexes_common_pipeline_${traitFem}.Rout"
echo Finished plotting Random-TBLUP

###Random Forest
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript Accuracy_RandomForest_Graphs_adjPheno_1trait_both_sexes_common_pipeline.R "$traitFem" > "Accuracy_RandomForest_Graphs_adjPheno_1trait_both_sexes_common_pipeline_${traitFem}.Rout"
echo Finished plotting Random Forest

###GO-TBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript Accuracy_GO-transcripts_Graphs_adjPheno_1trait_both_sexes_common_pipeline.R "$traitFem" > "Accuracy_GO-transcripts_Graphs_adjPheno_1trait_both_sexes_common_pipeline_${traitFem}.Rout"
echo Finished plotting GO-TBLUP

###GO-GBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript Accuracy_GO-SNPs_Graphs_adjPheno_1trait_both_sexes_common_pipeline.R "$traitFem" > "Accuracy_GO-SNPs_Graphs_adjPheno_1trait_both_sexes_common_pipeline_${traitFem}.Rout"
echo Finished plotting GO-GBLUP

###GO-GTBLUP
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript Accuracy_GO-SNPs_GO-transcripts_Graphs_adjPheno_1trait_both_sexes_common_pipeline.R "$traitFem" > "Accuracy_GO-SNPs_GO-transcripts_Graphs_adjPheno_1trait_both_sexes_common_pipeline_${traitFem}.Rout"
echo Finished plotting GO-GBLUP

###GO-TBLUP accuracy vs GO size
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript AccuracyVsGOsize_GO-transcripts_Graphs_adjPheno_1trait_both_sexes_common_pipeline.R "$traitFem" > "AccuracyVsGOsize_GO-transcripts_Graphs_adjPheno_1trait_both_sexes_common_pipeline_${traitFem}.Rout"
echo Finished plotting GO-TBLUP accuracy vs GO size

###GO-TBLUP accuracy vs variance explained
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript AccuracyVsVarianceExplained_GO-transcripts_Graphs_adjPheno_1trait_both_sexes_common_pipeline.R "$traitFem" > "AccuracyVsVarianceExplained_GO-transcripts_Graphs_adjPheno_1trait_both_sexes_common_pipeline_${traitFem}.Rout"
echo Finished plotting GO-TBLUP accuracy vs variance explained

###GO-GBLUP accuracy vs GO size
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript AccuracyVsGOsize_GO-SNPs_Graphs_adjPheno_1trait_both_sexes_common_pipeline.R "$traitFem" > "AccuracyVsGOsize_GO-SNPs_Graphs_adjPheno_1trait_both_sexes_common_pipeline_${traitFem}.Rout"
echo Finished plotting GO-GBLUP accuracy vs GO size

###GO-GBLUP accuracy vs variance explained
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript AccuracyVsVarianceExplained_GO-SNPs_Graphs_adjPheno_1trait_both_sexes_common_pipeline.R "$traitFem" > "AccuracyVsVarianceExplained_GO-SNPs_Graphs_adjPheno_1trait_both_sexes_common_pipeline_${traitFem}.Rout"
echo Finished plotting GO-GBLUP accuracy vs variance explained

###Variance partition in the training data set
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript Variance_Partition_Training_Graphs_adjPheno_1trait_both_sexes_common_pipeline.R "$traitFem" > "Variance_Partition_Training_Graphs_adjPheno_1trait_both_sexes_common_pipeline_${traitFem}.Rout"
echo Finished plotting variance partition in the training data

###Variance partition in the whole data set
/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript Variance_Partition_Whole_Graphs_adjPheno_1trait_both_sexes_common_pipeline.R "$traitFem" > "Variance_Partition_Whole_Graphs_adjPheno_1trait_both_sexes_common_pipeline_${traitFem}.Rout"
echo Finished plotting variance partition in the whole data