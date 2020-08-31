#!/bin/bash

########################
## Fabio Morgante
## 2/23/2018
##
## Modified on 4-5-2018
## Reason: use adjusted phenotype
##
## Modified on 11-11-2018
## Reason: use quotes around variables
########################

###Get trait from arguments
trait=$1

###Loop through pvalue threshold for TWAS

for PVAL in 0.000001 0.00001 0.0001 0.001 0.01 0.1 0.5
do
	#Print pvalue to be used
    echo "$PVAL"
    
    #Run the R script
	./Run_R_Open_2args.sh TBLUP_TWAS_LMM_regress_adjPheno_common_females "$trait" "$PVAL"
done