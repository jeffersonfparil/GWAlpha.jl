#!/bin/bash

#######################################################################
### GENOME-WIDE ASSOCIATION AND GENOMIC PREDICTION CROSS-VALIDATION ###
#######################################################################

#############
### INPUT ###
#############
OUTDIR=${1}             #output directory containing the parsed genotype and phenotype data
QTL_SPEC_FNAME=${2}     #full path and file name of the QTL specification file
GEN_PRED_SRC_DIR=${3}   #full path to the genomic_prediction/src directory
nCores=${4}							#total number of parallel jobs

########################
### SAMPLE EXECUTION ###
########################
# DIR=/data/Lolium/Quantitative_Genetics/LOLSIM_2019_TEST
# rep=1
# nQTL=10
# migration=0.00
# selection=0.25
# bg_selection=0.00
# GRADIENT=0
# BASEDIR=${DIR}/LOLSIM_${rep}rep_${nQTL}qtl_${migration}mr_${selection}fgs_${bg_selection}bgs_${GRADIENT}grad
# OUTPREFIX=LOLSIM_${rep}rep_${nQTL}qtl_${migration}mr_${selection}fgs_${bg_selection}bgs_${GRADIENT}grad
# OUTDIR=${BASEDIR}/${OUTPREFIX}*/
# QTL_SPEC_FNAME=${BASEDIR}/QTL_SPEC.csv
# GEN_PRED_SRC_DIR=${DIR}/Softwares/genomic_prediction/src
# nCore=12
#
# time \
# ${GEN_PRED_SRC_DIR}/GPASim_04_GWAS_GP.sh \
#   ${OUTDIR} \
#   ${QTL_SPEC_FNAME} \
#   ${GEN_PRED_SRC_DIR} \
#		${nCores}
########################

#############################
### SET WORKING DIRECTORY ###
#############################
cd ${OUTDIR}

##################################################
### EXECUTE GWAS & GP MODELS CROSS-VALIDATIONS ###
##################################################
echo -e "#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#"
echo -e "EXECUTE OF GWAS/GP:"
echo -e "#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#"

### Randomize then split into sqrt(nrow) files if nrow > 10 (arbitrary number >10 - seems big enough to warrant splitting eh?!)
ls POP_FNAME_*.csv > FNAMES_LIST.temp
for f in $(cat FNAMES_LIST.temp)
do
	NMERGED=$(echo $(cat $f | wc -l) -1 | bc)
	if [ $NMERGED -gt 10 ]
	then
		head -n1 $f > HEADER.temp
		tail -n+2 $f | shuf | split -d -l $(echo "sqrt($NMERGED)" | bc) - ${f%.csv*}_
		for split_files in $(ls ${f%.csv*}_*)
		do
			cat HEADER.temp > ${split_files}.csv
			cat $split_files >> ${split_files}.csv
			rm $split_files
		done
		rm $f
	fi
done
ls POP_FNAME_*.csv > FNAMES_LIST.temp

### In parallel or iteratively
parallel --jobs ${nCores} --link julia ${GEN_PRED_SRC_DIR}/GPASim_04_GWAS_GP.jl \
	${OUTDIR} \
	{1} \
	{2} \
	${QTL_SPEC_FNAME} \
	1 \
	Bonferroni \
	0.01 ::: $(cat FNAMES_LIST.temp) ::: $(cut -d_ -f4 FNAMES_LIST.temp) || for i in $(grep -vf <<<$(ls CROSS_VALIDATION*) FNAMES_LIST.temp)
do
  echo "#############################################"
  echo $i
  j=$(cut -d_ -f4 <<<$i)
  julia ${GEN_PRED_SRC_DIR}/GPASim_04_GWAS_GP.jl \
    ${OUTDIR} \
    ${i} \
    ${j} \
    ${QTL_SPEC_FNAME} \
    1 \
    Bonferroni \
    0.01
done

  ### INPUT ###
  ### (1) directory where the parsing output are written into
  ### (2) csv file of lists of individual/pool data per/across population/s
  ### ### HEADERS:
  ### ### POP,GENO,PHENO,COVARIATE_PC,COVARIATE_K for: POP_FNAME_WITHIN_INDI_DF.csv and POP_FNAME_ACROSS_INDI_DF.csv
  ### ### POP,GENO_SYNC,GENO_CSV,PHENO_PY,PHENO_CSV,COVARIATE_NPSTAT,COVARIATE_FST for: POP_FNAME_WITHIN_POOL_DF.csv or POP_FNAME_ACROSS_POOL_DF.csv
  ### (3) cross-validations to perform: {"INDI", "POOL"}
  ### (4) QTL identities (CHROM,POS,ALLELE,EFFECT)
  ### (5) size of linkage block in kilobases
  ### (6) p-value adjustment method: "Bonferroni", "BenjaminiHochberg" or "BenjaminiYekutieli"
  ### (7) significance level (decimal): prbability of incorrectly  rejecting the null hypothesis

  ### OUTPUT ###
  ### (1) CROSS_VALIDATION_OUTPUT_POP_FNAME_WITHIN_INDI_DF.csv
  ### (2) CROSS_VALIDATION_OUTPUT_POP_FNAME_WITHIN_POOL_DF.csv
  ### (3) CROSS_VALIDATION_OUTPUT_POP_FNAME_ACROSS_INDI_DF.csv
  ### (4) CROSS_VALIDATION_OUTPUT_POP_FNAME_ACROSS_POOL_DF.csv


####################
### MERGE OUTPUT ###
####################
echo -e "#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#"
echo -e "OUTPUT MERGING:"
echo -e "#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#"
echo -e "GROUPING,ALGORITHM,POP_TRAIN,POP_TEST,MODEL_ITERATION,MODEL_COVARIATE,MODEL_MODEL,PREDICTORS,NON_ZERO_PREDICTORS,MEAN_DEVIANCE,VAR_DEVIANCE,CORRELATION,INTERCEPT,SLOPE,R2,RMSD,TRUE_POSITIVE_RATE,QTL_FRAC_EFF,TRUE_POSITIVE_ID,FALSE_POSITIVE_ID,FALSE_POSITIVE_RATE,QTL_FREQS_TRAINING,QTL_FREQS_VALIDATION" \
  > CROSS_VALIDATION_OUTPUT_MERGED.csv
ls CROSS_VALIDATION_OUTPUT_*_*.csv > OUT_LIST.temp
for f in $(cat OUT_LIST.temp)
do
  GROUP=$(echo $f | cut -d_ -f6)
  ALGOR=$(echo $f | cut -d_ -f7)
  touch COL1_COL2.temp
  for i in $(seq 1 $(cat $f | wc -l)); do echo -e "$GROUP,$ALGOR" >> COL1_COL2.temp; done
  paste -d, COL1_COL2.temp $f > ${f%.*}_WITH_COL1_COL2.csv
  tail -n+2 ${f%.*}_WITH_COL1_COL2.csv >> CROSS_VALIDATION_OUTPUT_MERGED.csv
  rm COL1_COL2.temp
done

echo "#########################################################################################"
