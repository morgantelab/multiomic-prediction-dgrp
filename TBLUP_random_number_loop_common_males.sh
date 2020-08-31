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

###Loop through number of transcripts randomly sampled

for NUM in 5 50 500 1000 5000
do
	#Print pvalue to be used
    echo "$NUM"
    
    #Run the R script
	./Run_R_Open_2args.sh TBLUP_random_regress_adjPheno_common_males "$trait" "$NUM"
done