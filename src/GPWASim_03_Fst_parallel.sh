#!/bin/bash

### PARALLEL execution of:
### Fst calculations using Pool-seq data
### across pools or populations
### pairwise Fst per pool or population pair

GEN_PRED_SRC_DIR=$1
window_size=$2
fname_sync=$3
fname_pheno=$4
model=$5

julia ${GEN_PRED_SRC_DIR}/GPWASim_03_Fst.jl \
  ${fname_sync} \
  ${window_size} \
  $(cut -d',' -f 1 ${fname_pheno} | cut -d'.' -f1 | sed ':a;N;$!ba;s/\n/,/g' -) \
  ${model}
