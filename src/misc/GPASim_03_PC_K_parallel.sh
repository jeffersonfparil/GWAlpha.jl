#!/bin/bash

### PARALLEL execution of:
### PC and K covariates

GEN_PRED_SRC_DIR=$1
fname=$2
covariate=$3

julia ${GEN_PRED_SRC_DIR}/GPASim_03_PC_K_covariates.jl \
  ${fname} \
  $(echo "scale=10; 1/" $(cat ${fname%_GENO.csv*}_PHENO.csv | wc -l) | bc) \
  ${covariate}
