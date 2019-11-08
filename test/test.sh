#!/bin/bash

#############################
### TEST GPASim WORKFLOW ###
#############################

########################
### DEFINE VARIABLES ###
########################
DIR=/data/Lolium/Quantitative_Genetics/LOLSIM_2019_TEST
OUTDIR=${DIR}/Output/                                     # existing output directory to dump all the output files and folder
QUANTINEMO_DIR=${DIR}/Softwares/quantinemo_linux          # location of the quantinemo executable
nIndividuals=1000                                         # number of individuals to simulate
nLoci=1000                                                # number of loci to simulate (neutral + QTL)
nQTL=10                                                   # number of QTL among the all the loci simulated
nBGS=100                                                  # number of background selection loci (should be highly polygenic >= 100)
nAlleles=5                                                # number of alleles per loci, e.g. 5 for A,T,C,G, and DEL
allele_eff_model=CHISQ                                    # probability distribution to sample the allele effects from (i.e. CHISQ, NORM, & UNIF)
nGen=500                                                  # number of generations to simulate
# nPop=100                                                  # number of populations/subpopulations to simulate (NOTE: must have a natural number square-root)
nPop=16                                                  # number of populations/subpopulations to simulate (NOTE: must have a natural number square-root)
migration=0.001                                           # migration rate across the populations (surrently using the 1D stepping stone model see line 117)
selection=0.25                                             # selection intensity (trait of interest) defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
bg_selection=0.00                                         # background selection intensity defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
GRADIENT=1                                                # diffusion gradient of non-wildtype alleles
rep=1
prefix=LOLIUM
OUTPREFIX=${prefix}_${rep}rep_${nQTL}QTL_${migration}mr_${selection}fgs_${bg_selection}bgs_${GRADIENT}grad   # prefix for the quantiNemo2 initiation (*.ini) file and the output folder
fname_dat=${OUTPREFIX}_g${nGen}.dat                       # filename of the genotype (*.dat) file within the subdirectory of $OUTDIR
nPools=5                                                  # number of pools to subdivided each population
GEN_PRED_SRC_DIR=${DIR}/Softwares/genomic_prediction/src  # full path to the location of the julia 1.0.5 binary
INI_FNAME=${OUTPREFIX}.ini                                # full path and filename of the quantinemo initiation (*.ini) file
nCores=5                                                  # number of cores to use for parallelisation
OUT_SUBDIR=${OUTDIR}/${OUTPREFIX}*/                       # output directory containing the parsed genotype and phenotype data
POP_ID_txt=${OUT_SUBDIR}/POP_ID.txt                       # full path and file name of the text file containing the population ID (strings: e.g. 1 or 01 or 16) of each individual across populations - one ID per line
GENOME_SPEC_FNAME=${OUTDIR}/Lperenne_genome.spec          # full path and file name of the genome specification file
QTL_SPEC_FNAME=${OUTDIR}/QTL_SPEC.csv                     # full path and file name of the QTL specification file
NPSTAT_DIR=${DIR}/Softwares/npstat                        # full path to the npstat directory
N=500                                                     # maximum number of simulated libraries sequenced (e.g. N=500)


##################################
### GPASim_00_installation.sh ###
##################################
mkdir $DIR
mkdir $OUTDIR
time \
${GEN_PRED_SRC_DIR}/GPASim_00_installation.sh $DIR

##############################
### GPASim_01_simulate.sh ###
##############################
time \
${GEN_PRED_SRC_DIR}/GPASim_01_simulate.sh \
      $QUANTINEMO_DIR \
      $GEN_PRED_SRC_DIR \
      $OUTDIR \
      $OUTPREFIX \
      $nIndividuals \
      $nLoci \
      $nQTL \
      $nBGS \
      $nAlleles \
      $allele_eff_model \
      $nGen \
      $nPop \
      $migration \
      $selection \
      $bg_selection \
      $GRADIENT

###########################
### GPASim_02_parse.sh ###
###########################
time \
${GEN_PRED_SRC_DIR}/GPASim_02_parse.sh \
      ${OUTDIR} \
      ${fname_dat} \
      ${nPools} \
      ${GEN_PRED_SRC_DIR} \
      ${INI_FNAME} \
      ${nCores}

##################################
### GPASim_03_summaryStats.sh ###
##################################
time \
${GEN_PRED_SRC_DIR}/GPASim_03_summaryStats.sh \
  ${OUT_SUBDIR} \
  ${OUTPREFIX} \
  ${POP_ID_txt} \
  ${nGen} \
  ${GENOME_SPEC_FNAME} \
  ${QTL_SPEC_FNAME} \
  ${GEN_PRED_SRC_DIR} \
  ${NPSTAT_DIR} \
  ${N} \
  ${nCores}

#############################
### GPASim_04_GWAS_GP.sh ###
#############################
time \
${GEN_PRED_SRC_DIR}/GPASim_04_GWAS_GP.sh \
  ${OUT_SUBDIR} \
  ${QTL_SPEC_FNAME} \
  ${GEN_PRED_SRC_DIR} \
  ${nCores}
