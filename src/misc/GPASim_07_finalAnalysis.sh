#!/bin/bash

###############################################################
###  INTEGRATING ANALYSIS USING THE 2 SUMMAARY STATISTICS   ###
### OR GPAS METRICS (ie. LOG10_RMSD and TRUE_POSITIVE_RATE) ###
###############################################################

### ARGUMENTS
SRC_DIR=$1
DIR=$2
# ### TEST
# SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src
# DIR=/data/Lolium/Quantitative_Genetics/LOLSIM_2019
# # DIR=/data/Lolium/Quantitative_Genetics/LOLSIM_2019_BACKUP_ANALYSIS_20191231/BK_20191230

### SET WORKING DIRECTORY
cd $DIR

### (1) ABC_OPTIM
echo -e "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo -e "ABC Optimization Analysis"
echo -e "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
### prepare merged input file
head -n1 $(ls | grep LOLSIM_ | head -n1)/ABC_RESOURCE_OPTIM*JUST2SUMMSTATS.csv > ABC_MERGED_RESOURCE_OPTIM.csv
for i in $(ls | grep LOLSIM_)
do
  if [ $(ls ${i}/ABC_RESOURCE_OPTIM*JUST2SUMMSTATS.csv | wc -l) -eq 1 ]
  then
    tail -n+2 ${i}/ABC_RESOURCE_OPTIM*JUST2SUMMSTATS.csv >> ABC_MERGED_RESOURCE_OPTIM.csv
  fi
done
### execute script
time Rscript ${SRC_DIR}/GPASim_07_ABC.r ABC_MERGED_RESOURCE_OPTIM.csv

### (2) GPAS performance as a function of number of QTL, migration rates, background selection, and selection gradients
echo -e "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
echo -e "GPAS ~ Poplation + Quantitative Genetic Factors"
echo -e "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
### prepare merged input file
echo REP,QTL,MGR,FGS,BGS,GRAD,$(head -n1 $(ls | grep LOLSIM_ | head -n1)/CROSS_VALIDATION_OUTPUT_MERGED.csv | cut -d, -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,21) > MERGED.csv
### In julia: repeat factors
echo -e "using DelimitedFiles" > repeat_factors.jl
echo -e "factors_temp = DelimitedFiles.readdlm(ARGS[1], ',')" >> repeat_factors.jl
echo -e "nLines = parse(Int64, ARGS[2])" >> repeat_factors.jl
echo -e "# ### test" >> repeat_factors.jl
echo -e "# factors_temp = DelimitedFiles.readdlm(\"factors.temp\", ',')" >> repeat_factors.jl
echo -e "# nLines = parse(Int64, "35774")" >> repeat_factors.jl
echo -e "DelimitedFiles.writedlm(ARGS[1], repeat(factors_temp, nLines), ',')" >> repeat_factors.jl
### In R: remove auto-cross-validation datapoints
echo -e "args = commandArgs(trailingOnly=TRUE)" > remove_autoCV.r
echo -e "# args = c('LOLSIM_5rep_5qtl_0.001mr_0.25fgs_0.25bgs_2grad/CROSS_VALIDATION_OUTPUT_MERGED.csv')" >> remove_autoCV.r
echo -e "dat = read.csv(args[1])" >> remove_autoCV.r
echo -e "idx_noAutoCV = c()" >> remove_autoCV.r
echo -e "for (i in 1:nrow(dat)){" >> remove_autoCV.r
echo -e "  # i = 1" >> remove_autoCV.r
echo -e "  train_set = unlist(strsplit(as.character(dat\$POP_TRAIN[i]), ';'))" >> remove_autoCV.r
echo -e "  test_set = unlist(strsplit(as.character(dat\$POP_TEST[i]), ';'))" >> remove_autoCV.r
echo -e "  if (sum(train_set %in% test_set) + sum(test_set %in% train_set) == 0){" >> remove_autoCV.r
echo -e "    idx_noAutoCV = c(idx_noAutoCV, i)" >> remove_autoCV.r
echo -e "  }" >> remove_autoCV.r
echo -e "}" >> remove_autoCV.r
echo -e "out = dat[idx_noAutoCV, ]" >> remove_autoCV.r
echo -e "out_fname = paste0(dirname(args[1]), '/', gsub('.csv', '_filtered.csv', basename(args[1])))" >> remove_autoCV.r
echo -e "write.table(out, out_fname, row.names=FALSE, sep=',')" >> remove_autoCV.r
### merging parallelizable script
echo -e "i=\$(basename \$1)" > merging_paralel.sh
echo -e "# i=\$(ls | grep LOLSIM_ | head -n1)" >> merging_paralel.sh
echo -e "echo \$i" >> merging_paralel.sh
echo -e "if [ \$(ls \${i}/CROSS_VALIDATION_OUTPUT_MERGED.csv | wc -l) -eq 1 ] && [ \$(ls \${i}/ERRORED | wc -l) -eq 0 ]" >> merging_paralel.sh
echo -e "then" >> merging_paralel.sh
echo -e "  REP=\$(echo \${i} | cut -d_ -f2 | sed 's/rep//g')" >> merging_paralel.sh
echo -e "  QTL=\$(echo \${i} | cut -d_ -f3 | sed 's/qtl//g')" >> merging_paralel.sh
echo -e "  MGR=\$(echo \${i} | cut -d_ -f4 | sed 's/mr//g')" >> merging_paralel.sh
echo -e "  FGS=\$(echo \${i} | cut -d_ -f5 | sed 's/fgs//g')" >> merging_paralel.sh
echo -e "  BGS=\$(echo \${i} | cut -d_ -f6 | sed 's/bgs//g')" >> merging_paralel.sh
echo -e "  GRAD=\$(echo \${i} | cut -d_ -f7 | sed 's/grad//g')" >> merging_paralel.sh
echo -e "  echo -e \${REP},\${QTL},\${MGR},\${FGS},\${BGS},\${GRAD} > \${i}_factors.temp" >> merging_paralel.sh
echo -e "  Rscript remove_autoCV.r \${i}/CROSS_VALIDATION_OUTPUT_MERGED.csv" >> merging_paralel.sh
echo -e "  tail -n+2 \${i}/CROSS_VALIDATION_OUTPUT_MERGED_filtered.csv | cut -d, -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,21 > \${i}_data.temp" >> merging_paralel.sh
echo -e "  julia repeat_factors.jl \${i}_factors.temp \$(cat \${i}_data.temp | wc -l)" >> merging_paralel.sh
echo -e "  paste -d, \${i}_factors.temp \${i}_data.temp > \${i}_MERGED.csv" >> merging_paralel.sh
echo -e "  rm \${i}_*.temp" >> merging_paralel.sh
echo -e "fi" >> merging_paralel.sh
chmod +x merging_paralel.sh
time parallel -j 5 ./merging_paralel.sh {} ::: $(ls | grep LOLSIM_)
cat *_MERGED.csv >> MERGED.csv
rm *_MERGED.csv
### (2.1) Semi-Exhaustive pairwise mixed models
Rscript GPASim_07_exhaustive_pairwise_mixed_models.r MERGED.csv
### (2.2)  Sensible mixed modes ...
### ... where we are interested in the effects of:
###   - the number of training populations and
###   - the different GPAS models
###   - as well as the number of simulated QTL for within population data sets
### while controlling for the random (assumed random) effects of:
###   - the number of  simulated QTL (applicable only for across population datasets),
###   - migration rate,
###   - background selection intensity and
###   - resistance allele emergence gradient.
### We are trying to peer more closely with eyes for sensible models (create an Rscript ("sensible_mixed_modelling.R") with the following contents:)
Rscript GPASim_07_sensible_mixed_models.r GPASim_data.rds
