#!/bin/bash

###################################################
### CALCULATE SUMMARY STATISTICS AND COVARIATES ###
###################################################

#############
### INPUT ###
#############
OUTDIR=${1}             #output directory containing the parsed genotype and phenotype data
OUTPREFIX=${2}          #prefix of output data (same name as the INI file and the prefix of the output directory)
POP_ID_txt=${3}         #full path and file name of the text file containing the population ID (strings: e.g. 1 or 01 or 16) of each individual across populations - one ID per line
nGen=${4}               #number of generations
GENOME_SPEC_FNAME=${5}  #full path and file name of the genome specification file
QTL_SPEC_FNAME=${6}     #full path and file name of the QTL specification file
GEN_PRED_SRC_DIR=${7}   #full path to the genomic_prediction/src directory
NPSTAT_DIR=${8}         #full path to the npstat directory
N=${9}                  #maximum number of simulated libraries sequenced (e.g. N=500)
nCores=${10}            #number of cores to use for parallel computations of PC and K covariate (and not for NPSTAT summary statistics and Fst: these use all available cores)

########################
### SAMPLE EXECUTION ###
########################
# DIR=/data/Lolium/Quantitative_Genetics/LOLSIM_2019_TEST
# rep=1
# nQTL=10
# migration=0.00
# selection=0.25
# bg_selection=0.00
# VAR_gradient=0
# BASEDIR=${DIR}/LOLSIM_${rep}rep_${nQTL}qtl_${migration}mr_${selection}fgs_${bg_selection}bgs_${VAR_gradient}grad
# OUTPREFIX=LOLSIM_${rep}rep_${nQTL}qtl_${migration}mr_${selection}fgs_${bg_selection}bgs_${VAR_gradient}grad
# OUTDIR=${BASEDIR}/${OUTPREFIX}*/
# POP_ID_txt=${OUTDIR}/POP_ID.txt
# nGen=100
# GENOME_SPEC_FNAME=${BASEDIR}/Lperenne_genome.spec
# QTL_SPEC_FNAME=${BASEDIR}/QTL_SPEC.csv
# GEN_PRED_SRC_DIR=${DIR}/Softwares/genomic_prediction/src
# NPSTAT_DIR=${DIR}/Softwares/npstat
# N=100
# nCores=$(echo $(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l) / 2 | bc) #half the number of parallel threads for RAM-intensive PC and K covariate caluculations
#
# time \
# ${GEN_PRED_SRC_DIR}/GPASim_03_summaryStats.sh \
#   ${OUTDIR} \
#   ${OUTPREFIX} \
#   ${POP_ID_txt} \
#   ${nGen} \
#   ${GENOME_SPEC_FNAME} \
#   ${QTL_SPEC_FNAME} \
#   ${GEN_PRED_SRC_DIR} \
#   ${NPSTAT_DIR} \
#   ${N} \
#   ${nCores}
########################

#############################
### SET WORKING DIRECTORY ###
#############################
cd ${OUTDIR}

###################################################################
### LIST GENOTYPE AND PHENOTYPE FILES (WITHIN POPULATION TESTS) ###
###################################################################
echo -e "###############################################################################"
echo -e "LISTING INPUT FILENAMES FOR WITHIN POPULATION GWAS/GP:"
echo -e "###############################################################################"
nPop=$(uniq $POP_ID_txt | wc -l)
echo -e "POP,GENO,PHENO,COVARIATE_PC,COVARIATE_K" > POP_FNAME_WITHIN_INDI_DF.csv
echo -e "POP,GENO_SYNC,GENO_CSV,PHENO_PY,PHENO_CSV,COVARIATE_NPSTAT,COVARIATE_FST" > POP_FNAME_WITHIN_POOL_DF.csv

for i in $(uniq $POP_ID_txt)
do
  id=p${i} ### adding "p" prefix to the population ID to avoid misidentification as Int instead of String!

  fname_geno_indi=${OUTPREFIX}_g${nGen}_p${i}_GENO.csv
  fname_pheno_indi=${OUTPREFIX}_g${nGen}_p${i}_PHENO.csv
  fname_covar_PC_indi=${OUTPREFIX}_g${nGen}_p${i}_GENO_COVARIATE_PC.csv
  fname_covar_K_indi=${OUTPREFIX}_g${nGen}_p${i}_GENO_COVARIATE_K.csv

  fname_geno_pool_sync=${OUTPREFIX}_g${nGen}_p${i}_POOLS_GENO.sync
  fname_geno_pool_csv=${OUTPREFIX}_g${nGen}_p${i}_POOLS_GENO_ALLELEFREQ.csv ### HAVE NOT BEEN GENERATED YET --> similar to ${OUTPREFIX}_g${nGen}_p${i}_POOLS_GENO.csv --> but missing the 'N' column AND! generated separately by PoolGPAS_module for standallones in the future 20190924
  fname_pheno_pool_py=${OUTPREFIX}_g${nGen}_p${i}_POOLS_PHENO.py
  fname_pheno_pool_csv=${OUTPREFIX}_g${nGen}_p${i}_POOLS_PHENO.csv
  fname_covar_npstat_pool=${OUTPREFIX}_g${nGen}_p${i}_POOLS_GENO_COVARIATE_NPSTAT.csv
  fname_covar_fst_pool=${OUTPREFIX}_g${nGen}_p${i}_POOLS_GENO_COVARIATE_FST.csv
  n=$(awk -F, '{print NF}' <(head -n1 $fname_geno_indi))
  echo -e "${id},${fname_geno_indi},${fname_pheno_indi},${fname_covar_PC_indi},${fname_covar_K_indi}" >> POP_FNAME_WITHIN_INDI_DF.csv
  echo -e "${id},${fname_geno_pool_sync},${fname_geno_pool_csv},${fname_pheno_pool_py},${fname_pheno_pool_csv},${fname_covar_npstat_pool},${fname_covar_fst_pool}" >> POP_FNAME_WITHIN_POOL_DF.csv
done
### OUTPUTS
# (1) POP_FNAME_WITHIN_INDI_DF.csv (nPop merged populations)
# (2) POP_FNAME_WITHIN_POOL_DF.csv (nPop merged populations)

####################################################################
### LIST GENOTYPE AND PHENOTYPE FILES (ACROSS POPULATIONS TESTS) ###
####################################################################
echo -e "###############################################################################"
echo -e "LISTING INPUT FILENAMES FOR ACROSS POPULATIONS GWAS/GP:"
echo -e "###############################################################################"
julia ${GEN_PRED_SRC_DIR}/GPASim_03_GWAS_GP_pop_grouping.jl \
  POP_FNAME_WITHIN_INDI_DF.csv \
  ${OUTPREFIX}_ALLPOP_GENO.sync \
  ${OUTPREFIX}_ALLPOP_PHENO.pool \
  ${N}
### OUTPUTS
# (1) POP_FNAME_ACROSS_INDI_RANDUNIF_DF.csv (if mod(nPop,2)==0 then {sqrt(nPop)/2} + 1 else {sqrt(nPop)+1}/2)
# (2) POP_FNAME_ACROSS_POOL_RANDUNIF_DF.csv (if mod(nPop,2)==0 then {sqrt(nPop)/2} else [{sqrt(nPop)+1}/2] - 1)
# (3) POP_FNAME_ACROSS_INDI_PSEUDOPT_DF.csv (if mod(nPop,2)==0 then {sqrt(nPop)/2} + 1 else {sqrt(nPop)+1}/2)
# (4) POP_FNAME_ACROSS_POOL_PSEUDOPT_DF.csv (if mod(nPop,2)==0 then {sqrt(nPop)/2} else [{sqrt(nPop)+1}/2] - 1)
# (n) *MERGED* genotype and phenotype files

###############################################
### COVARIATES FOR INDIVIDUAL GENOTYPE DATA ###
###############################################
#####################
### WITHIN POPULATION
#####################
echo -e "################################################################################################"
echo -e "BUILDING COVARIATES (INDIVIDUAL-GENOTYPE-DERIVED: PC & K) FOR WITHIN POPULATION GWAS/GP:"
echo -e "################################################################################################"
parallel -j $nCores ${GEN_PRED_SRC_DIR}/GPASim_03_PC_K_parallel.sh \
  ${GEN_PRED_SRC_DIR} \
  ${OUTPREFIX}_g${nGen}_p{1}_GENO.csv \
  {2} ::: $(uniq $POP_ID_txt) ::: PC K
### using just the first PC: --->  which therefore means we need to remove PC COVARIATES from the mixed (RR, GLMNET, LASSO) and EMMAX models
for i in $(uniq $POP_ID_txt)
do
  echo $i
  cp ${OUTPREFIX}_g${nGen}_p${i}_GENO_COVARIATE_PC.csv ${OUTPREFIX}_g${nGen}_p${i}_GENO_COVARIATE_PC.bk
  awk -F, '{print $1}' ${OUTPREFIX}_g${nGen}_p${i}_GENO_COVARIATE_PC.csv > ${OUTPREFIX}_g${nGen}_p${i}_GENO_COVARIATE_PC.temp
  mv ${OUTPREFIX}_g${nGen}_p${i}_GENO_COVARIATE_PC.temp ${OUTPREFIX}_g${nGen}_p${i}_GENO_COVARIATE_PC.csv
done
#####################
### ACROSS POPULATION
#####################
echo -e "################################################################################################"
echo -e "BUILDING COVARIATES (INDIVIDUAL-GENOTYPE-DERIVED: PC & K) FOR ACROSS POPULATIONS GWAS/GP:"
echo -e "################################################################################################"
ls *_MERGED_GENO.csv > fnames_merged_geno_indi.temp
parallel -j $nCores ${GEN_PRED_SRC_DIR}/GPASim_03_PC_K_parallel.sh \
  ${GEN_PRED_SRC_DIR} \
  {1} \
  {2} ::: $(cat fnames_merged_geno_indi.temp) ::: PC K
for i in $(cat fnames_merged_geno_indi.temp)
do
  echo ${i%.csv*}_COVARIATE_PC.csv
  cp ${i%.csv*}_COVARIATE_PC.csv ${i%.csv*}_COVARIATE_PC.bk
  awk -F, '{print $1}' ${i%.csv*}_COVARIATE_PC.csv > ${i%.csv*}_COVARIATE_PC.temp
  mv ${i%.csv*}_COVARIATE_PC.temp ${i%.csv*}_COVARIATE_PC.csv
done
##############
### OUTPUT ###
##############
# (1) *_COVARIATE_PC.csv ### principal components: nx1
# (2) *_COVARIATE_K.csv  ### kinship matrix: nxn

#########################################################
### NPSTAT SUMMARY STATISTICS PER POOL PER POPULATION ###
#########################################################
echo -e "################################################################################################"
echo -e "BUILDING COVARIATES (POOL-SEQ-DERIVED: NPSTAT) FOR WITHIN POPULATION GWAS/GP:"
echo -e "################################################################################################"
echo -e "args = commandArgs(trailing=TRUE)
fname = args[1]
stats = read.table(fname, header=TRUE)
stats\$Scoverage = stats\$S / stats\$coverage
stats = data.frame(pool=stats\$pool, Scoverage=stats\$Scoverage, theta=stats\$Watterson_estimator, Pi=stats\$Tajimas_Pi, D=stats\$Tajimas_D, H=stats\$Fay_Wu_H_norm)
out_MEAN = aggregate(. ~ pool, data=stats, FUN=mean)[,2:ncol(stats)]
out_MIN = aggregate(. ~ pool, data=stats, FUN=min)[,2:ncol(stats)]
out_MAX = aggregate(. ~ pool, data=stats, FUN=max)[,2:ncol(stats)]
out_VAR = aggregate(. ~ pool, data=stats, FUN=var)[,2:ncol(stats)]
# ### testing Tajima's D parameter estimation (D ~ beta(a1, a2))
# par = optim(par=c(0.1, 0.1), fn=function(par) sum(dbeta(par[1], par[2], stats$Tajimas_D), na.rm=TRUE), lower=c(0.0, 0.0), upper=c(Inf, Inf))$par
# mean_D = par[1] / (par[1] + par[2])
OUT = cbind(out_MEAN, out_MIN, out_MAX, out_VAR)
colnames(OUT) = paste(rep(c('MEAN', 'MIN', 'MAX', 'VAR'), each=5), rep(c('Scoverage', 'theta', 'Pi', 'D', 'H'), times=4), sep='_')
outname = paste0(unlist(strsplit(fname, '_NPSTAT.stats'))[1], '_COVARIATE_NPSTAT.csv')
write.table(OUT, file=outname, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=',')
" > aggregate_npstats.r ### summarize npstat output per pool
nPools=$(cat $(ls ${OUTPREFIX}_g${nGen}_p*_POOLS_PHENO.csv | head -n1) | wc -l) ### same number of pools in all populations
nChrom=$(echo $(cat ${GENOME_SPEC_FNAME} | wc -l) - 1 | bc)
window_size=10000000 #10 Mb
for i in $(uniq $POP_ID_txt)
do
  echo POP_${i}
  for j in $(seq 1 $nPools)
  do
    # popool_size=$(cut -d, -f1 ${OUTPREFIX}_g${nGen}_p${i}_POOLS_PHENO.csv | head -n${j} | tail -n1)
    popool_size=10 #fixed to 10 for quick estimations though highly over-estimated summary statistics
    parallel \
    ${GEN_PRED_SRC_DIR}/GPASim_03_NPSTAT_parallel.sh \
            {1} \
            ${OUTPREFIX}_g${nGen}_p${i}_POOL${j}.pileup \
            ${popool_size} \
            ${window_size} \
            ${NPSTAT_DIR} ::: $(seq 1 $nChrom)
    ##############
    ### INPUT ####
    ##############
    # (1) chrom=1
    # (2) pileup_fname=QUANTI_g140_p*_POOL*.pileup
    # (4) n=popool_size
    # (5) NPSTAT_DIR=${NPSTAT_DIR}
    ##############
    ### OUTPUT ###
    ##############
    # (1) ${OUTPREFIX}_g${nGen}_p*_POOL*_CHR*_NPSTAT.stats (NO HEADER)
    #############
    ### MERGE ###
    #############
    # echo -e "chrom\twindow\tcoverage\toutgroup_coverage\tdepth\tS\tWatterson_estimator\tTajimas_Pi\tTajimas_D\tFay_Wu_H_unorm\tFay_Wu_H_norm\tS_var\tV_theta\toutgroup_divergence\tnonsynonimous_polymorph\tsynonymous_polymorph\tnonsynonimous_divergence\tsynonymous_divergence\talpha" > ${OUTPREFIX}_g${nGen}_p${i}_POOL${j}_NPSTAT.stats
    # cat ${OUTPREFIX}_g${nGen}_p${i}_POOL${j}_CHR*_NPSTAT.stats >> ${OUTPREFIX}_g${nGen}_p${i}_POOL${j}_NPSTAT.stats
    cat ${OUTPREFIX}_g${nGen}_p${i}_POOL${j}_CHR*_NPSTAT.stats > ${OUTPREFIX}_g${nGen}_p${i}_POOL${j}_NPSTAT.CHROM_MERGED
    touch ${OUTPREFIX}_g${nGen}_p${i}_POOL${j}_COL1_POOL_ID.temp
    for k in $(seq 1 $(cat ${OUTPREFIX}_g${nGen}_p${i}_POOL${j}_NPSTAT.CHROM_MERGED | wc -l)); do echo -e "$j" >> ${OUTPREFIX}_g${nGen}_p${i}_POOL${j}_COL1_POOL_ID.temp; done
    paste ${OUTPREFIX}_g${nGen}_p${i}_POOL${j}_COL1_POOL_ID.temp ${OUTPREFIX}_g${nGen}_p${i}_POOL${j}_NPSTAT.CHROM_MERGED > ${OUTPREFIX}_g${nGen}_p${i}_POOL${j}_NPSTAT.WITH_POOL_ID ### add pool ID
    rm *.temp *.CHROM_MERGED
  done
  echo -e "pool\tchrom\twindow\tcoverage\toutgroup_coverage\tdepth\tS\tWatterson_estimator\tTajimas_Pi\tTajimas_D\tFay_Wu_H_unorm\tFay_Wu_H_norm\tS_var\tV_theta\toutgroup_divergence\tnonsynonimous_polymorph\tsynonymous_polymorph\tnonsynonimous_divergence\tsynonymous_divergence\talpha" \
    > ${OUTPREFIX}_g${nGen}_p${i}_POOLS_GENO_NPSTAT.stats
  cat ${OUTPREFIX}_g${nGen}_p${i}_POOL*_NPSTAT.WITH_POOL_ID >> ${OUTPREFIX}_g${nGen}_p${i}_POOLS_GENO_NPSTAT.stats
  Rscript aggregate_npstats.r ${OUTPREFIX}_g${nGen}_p${i}_POOLS_GENO_NPSTAT.stats
  # output: ${OUTPREFIX}_g${nGen}_p${i}_POOLS_GENO_COVARIATE_NPSTAT.csv
  rm *.WITH_POOL_ID
  ##############
  ### OUTPUT ###
  ##############
  # (1) ${OUTPREFIX}_g${nGen}_p${i}_POOLS_GENO_COVARIATE_NPSTAT.csv
  # NO HEADER: COL1:S/COVERAGE, COL2:theta_w, COL3:Pi, COL4:D, COL5:Fay_WU_H_norm
done

################################################
### NPSTAT SUMMARY STATISTICS PER POPULATION ###
################################################
nChrom=$(echo $(cat ${GENOME_SPEC_FNAME} | wc -l) - 1 | bc)
window_size=10000000 #10 Mb
for i in $(uniq $POP_ID_txt)
do
  # pop_size=$(grep "^$i" $POP_ID_txt | wc -l)
  pop_size=10 #fixed to 10 for quick estimations though highly over-estimated summary statistics
  parallel \
  ${GEN_PRED_SRC_DIR}/GPASim_03_NPSTAT_parallel.sh \
          {1} \
          ${OUTPREFIX}_g${nGen}_p${i}_POPULATION.pileup \
          ${pop_size} \
          ${window_size} \
          ${NPSTAT_DIR} ::: $(seq 1 $nChrom)
  ##############
  ### INPUT ####
  ##############
  # (1) chrom=1
  # (2) pileup_fname=QUANTI_g140_p*_POPULATION.pileup
  # (4) n=pop_size
  # (5) NPSTAT_DIR=${NPSTAT_DIR}
  ##############
  ### OUTPUT ###
  ##############
  # (1) ${OUTPREFIX}_g${nGen}_p*_POPULATION_CHR*_NPSTAT.stats (NO HEADER)
  #############
  ### MERGE ###
  #############
  echo -e "chrom\twindow\tcoverage\toutgroup_coverage\tdepth\tS\tWatterson_estimator\tTajimas_Pi\tTajimas_D\tFay_Wu_H_unorm\tFay_Wu_H_norm\tS_var\tV_theta\toutgroup_divergence\tnonsynonimous_polymorph\tsynonymous_polymorph\tnonsynonimous_divergence\tsynonymous_divergence\talpha" > ${OUTPREFIX}_g${nGen}_p${i}_POPULATION_NPSTAT.stats
  cat ${OUTPREFIX}_g${nGen}_p${i}_POPULATION_CHR*_NPSTAT.stats >> ${OUTPREFIX}_g${nGen}_p${i}_POPULATION_NPSTAT.stats
  rm ${OUTPREFIX}_g${nGen}_p${i}_POPULATION_CHR*_NPSTAT.stats
done

############################################ ### NOTE: TESTING UP TO THE LINE ABOVE 20190928
### ACROSS POPULATIONS NPSTAT OUTPUT MERGING
############################################
echo -e "################################################################################################"
echo -e "BUILDING COVARIATES (POOL-SEQ-DERIVED: NPSTAT) FOR ACROSS POPULATIONS GWAS/GP:"
echo -e "################################################################################################"
ls *_MERGED_POOLS_GENO.sync > fnames_merged_geno_pool.temp
for i in $(cat fnames_merged_geno_pool.temp)
do
  # i=$(head -n1 fnames_merged_geno_pool.temp)
  echo $i
  # echo ${i%_MERGED_POOLS_GENO.sync*} > pops.temp
  nmatches=$(grep -c $(echo $i) POP_FNAME_ACROSS_POOL_RANDUNIF_DF.csv )
  if [ $nmatches != 0 ]
  then
    grep $(echo $i) POP_FNAME_ACROSS_POOL_RANDUNIF_DF.csv | cut -d, -f1 -  > pops.temp
  else
    grep $(echo $i) POP_FNAME_ACROSS_POOL_PSEUDOPT_DF.csv | cut -d, -f1 -  > pops.temp
  fi
  # sed -i 's/_/\n/g' pops.temp
  sed -i 's/;/\n/g' pops.temp
  echo -e "pool\tchrom\twindow\tcoverage\toutgroup_coverage\tdepth\tS\tWatterson_estimator\tTajimas_Pi\tTajimas_D\tFay_Wu_H_unorm\tFay_Wu_H_norm\tS_var\tV_theta\toutgroup_divergence\tnonsynonimous_polymorph\tsynonymous_polymorph\tnonsynonimous_divergence\tsynonymous_divergence\talpha" \
    > ${i%.sync*}_NPSTAT.stats
  for j in $(cat pops.temp)
  do
    # j=$(head -n1 pops.temp)
    tail -n+2 ${OUTPREFIX}_g${nGen}_${j}_POPULATION_NPSTAT.stats > ${OUTPREFIX}_g${nGen}_${j}_POPULATION_NPSTAT.noheader.temp
    touch ${j}_COL1_POOL_ID.temp
    for k in $(seq 1 $(cat ${OUTPREFIX}_g${nGen}_${j}_POPULATION_NPSTAT.noheader.temp | wc -l)); do echo -e "$j" >> ${j}_COL1_POOL_ID.temp; done
    paste -d'\t' ${j}_COL1_POOL_ID.temp ${OUTPREFIX}_g${nGen}_${j}_POPULATION_NPSTAT.noheader.temp >> ${i%.sync*}_NPSTAT.stats
    rm *.temp
  done
  Rscript aggregate_npstats.r ${i%.sync*}_NPSTAT.stats
  # output: ${i%.sync*}_COVARIATE_NPSTAT.csv
  rm ${i%.sync*}_NPSTAT.stats
done
##############
### OUTPUT ###
##############
# (1) *MERGED_GENO_COVARIATE_NPSTAT.csv


###############################################
###                                         ###
###             FST ACROSS POOLS            ###
###                                         ###
###############################################
#################################################
### PAIRWISE FST BETWEEN POOLS PER POPULATION ###
#################################################
echo -e "################################################################################################"
echo -e "BUILDING COVARIATES (POOL-SEQ-DERIVED: FST) FOR WITHIN POPULATION GWAS/GP:"
echo -e "################################################################################################"
window_size=1000000 #1Mb
parallel ${GEN_PRED_SRC_DIR}/GPASim_03_Fst_parallel.sh \
  ${GEN_PRED_SRC_DIR} \
  ${window_size} \
  ${OUTPREFIX}_g${nGen}_p{1}_POOLS_GENO.sync \
  ${OUTPREFIX}_g${nGen}_p{1}_POOLS_PHENO.csv \
  WeirCock ::: $(uniq $POP_ID_txt)
######################################################
### PAIRWISE FST BETWEEN POPULATIONS ACROSS GROUPS ###
######################################################
echo -e "################################################################################################"
echo -e "BUILDING COVARIATES (POOL-SEQ-DERIVED: FST) FOR ACROSS POPULATIONS GWAS/GP:"
echo -e "################################################################################################"
window_size=1000000 #1Mb
parallel --link ${GEN_PRED_SRC_DIR}/GPASim_03_Fst_parallel.sh \
  ${GEN_PRED_SRC_DIR} \
  ${window_size} \
  {1} \
  {2} \
  WeirCock ::: $(ls *_MERGED_POOLS_GENO.sync) ::: $(ls *_MERGED_POOLS_PHENO.csv)
  #############
  ### INPUT ###
  #############
  # (1) directory where the parallelization script is located
  # (2) window size in bases
  # (3) pool sync file
  # (4) pool phenotype file
  # (5) Fst calculation method to use: WeirCock or Hivert
  ##############
  ### OUTPUT ###
  ##############
  # (1) whole population Fst per window = string(split(sync_fname, ".")[1], "_window_", window_size, "bp_Fst_data.csv")
  # (2) whole population Fst per pool = string(split(sync_fname, ".")[1], "_window_", window_size, "bp_Fst_sumstats.csv")
  # (3) Pairwise Fst = string(split(sync_fname, ".")[1], "_COVARIATE_FST.csv")
  echo "#########################################################################################"
