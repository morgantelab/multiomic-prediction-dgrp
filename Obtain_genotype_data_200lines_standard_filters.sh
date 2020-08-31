#!/bin/bash

########################
## Fabio Morgante
## 4/3/2018
## Extract genotype data for 200 lines with
## expression data, MAF 0.05, call rate 0.8 
##
## Modified on 05/06/2020
## Reason: Use relative paths             
#############################


../Software/plink_linux_x86_64/plink --noweb --bfile ../DGRPdata/dgrp2 --keep Data/lines_expression_plink_format.txt --maf 0.05 --geno 0.2 --make-bed --out Data/dgrp2_200lines_common_prediction
