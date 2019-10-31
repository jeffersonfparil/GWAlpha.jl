#!/bin/bash

##############
### INPUTS ###
##############
DIR=${1}
VAR_nQTL=${2}
VAR_migration=${3}
VAR_foreground_selection=${4}
VAR_background_selection=${5}
VAR_gradient=${6}
VAR_rep=${7}
VAR_partition=${8}
VAR_prefix=${9}
# DIR=/data/cephfs/punim0543/jparil
# VAR_nQTL=10                     # 5, 10, 100
# VAR_migration=0.001             # 0.01, 0.001, 0.0001
# VAR_foreground_selection=0.25   # 0.0, 0.25
# VAR_background_selection=-0.25  # 0.0, 0.25, -0.25
# VAR_gradient=0                  # 0, 1, 2
# VAR_rep=1                       # 1, 2, 3, 4, 5
# VAR_partition=mig               # cloud, physical, snowy, ..., etc

########################
### SAMPLE EXECUTION ###
# ./GPASim_05_slurmer.sh /data/cephfs/punim0543/jparil 100 0.001 0.25 -0.25 0 1 mig LOLSIM
########################

########################################
### AUXILLARY NON-USER-FACING INPUTS ###
########################################
VAR_indiv=1000                                                  ### number of individuals per poplation to simulate
VAR_nLoci=10000                                                 ### number of loci per individual to simulate
VAR_nBGS=100                                                    ### number of background selection QTL
VAR_nGen=500                                                    ### number of generations to simulate
VAR_nPop=100                                                    ### number of populations (m x m: should have a natural number square-root)
VAR_nPools=5                                                    ### number of pools for within population Pool-seq
VAR_maxlib=500                                                  ### maximum number of genotyping libraries
VAR_account=punim0543                                           ### Spartan account name
VAR_reqmem=$(echo "${VAR_nLoci} * ${VAR_indiv} / 1000000" | bc) ### estimated RAM requirement in gigabytes per parallel job during parsing, summary statistics including PC and K estimations and GPAS CV
if [ $VAR_partition == cloud ]
then
  VAR_ncores=12                                                 ### number of cores avaible given the partition
  VAR_mem=95 ### 100 really but I don't want to push it
elif [ $VAR_partition == longcloud ]
then
  VAR_ncores=12
  VAR_mem=95
elif [ $VAR_partition == physical ]
then
  VAR_ncores=32
  VAR_mem=250
elif [ $VAR_partition == bigmem ]
then
  VAR_ncores=36
  VAR_mem=1535
elif [ $VAR_partition == snowy ]
then
  VAR_ncores=32
  VAR_mem=120
elif [ $VAR_partition == mig ]
then
  VAR_ncores=32
  VAR_mem=120
  VAR_account=punim0594
else
  echo -e "Unrecognised Spartan high performance computing partition."
  echo -e "Select from:"
  echo -e "\t - cloud"
  echo -e "\t - longcloud"
  echo -e "\t - physical"
  echo -e "\t - bigmem"
  echo -e "\t - snowy"
  echo -e "\t - mig"
  exit
fi

echo -e '#!/bin/bash' > ${VAR_prefix}_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.slurm
echo -e "
# Partition for the job:
#SBATCH --partition=${VAR_partition}
# Multithreaded (SMP) job: must run on one node and the cloud partition
#SBATCH --nodes=1
# The name of the job:
#SBATCH --job-name=${VAR_prefix}
# The project ID which this job should run under:
#SBATCH --account=${VAR_account}
# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${VAR_ncores}
# The amount of memory in megabytes per process in the job:
#SBATCH --mem=${VAR_mem}GB
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=5-0:0:00
# Send yourself an email when the job:
# aborts abnormally (fails)
# #SBATCH --mail-type=FAIL
# begins
# #SBATCH --mail-type=BEGIN
# ends successfully
# #SBATCH --mail-type=END
# Use this email address:
# #SBATCH --mail-user=jparil@student.unimelb.edu.au

echo '#######################################################'
echo 'DEFINE VARIABLES'
echo '#######################################################'
QUANTINEMO_DIR=${DIR}/Softwares/quantinemo_linux          # location of the quantinemo executable
GEN_PRED_SRC_DIR=${DIR}/Softwares/genomic_prediction/src  # full path to the genomic_prediction/src directory
nIndividuals=${VAR_indiv}                                 # number of individuals to simulate
nLoci=${VAR_nLoci}                                        # number of loci to simulate (neutral + QTL)
nQTL=${VAR_nQTL}                                          # number of QTL among the all the loci simulated
nBGS=${VAR_nBGS}                                          # number of background selection loci (should be highly polygenic >= 100)
nAlleles=5                                                # number of alleles per loci, e.g. 5 for A,T,C,G, and DEL
allele_eff_model=CHISQ                                    # probability distribution to sample the allele effects from (i.e. CHISQ, NORM, & UNIF)
nGen=${VAR_nGen}                                          # number of generations to simulate
nPop=${VAR_nPop}                                          # number of populations/subpopulations to simulate (NOTE: must have a natural number square-root)
migration=${VAR_migration}                                # migration rate across the populations (surrently using the 1D stepping stone model see line 117)
selection=${VAR_foreground_selection}                     # selection intensity (trait of interest) defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
bg_selection=${VAR_background_selection}                  # background selection intensity defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
GRADIENT=${VAR_gradient}                                  # diffusion gradient of non-wildtype alleles
OUTPREFIX=${VAR_prefix}_${VAR_rep}rep_${VAR_nQTL}qtl_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad  # prefix for the quantiNemo2 initiation (*.ini) file and the output folder
OUTDIR=${DIR}/\${OUTPREFIX}; if [ ! -d \${OUTDIR} ]; then mkdir \${OUTDIR}; fi                                         # create the output directory
fname_dat=\${OUTPREFIX}_g${VAR_nGen}.dat                       # filename of the genotype (*.dat) file within the subdirectory of $OUTDIR
nPools=${VAR_nPools}                                      # number of pools to subdivided each population
INI_FNAME=\${OUTPREFIX}.ini                                # full path and filename of the quantinemo initiation (*.ini) file
nCores=$(echo "${VAR_mem} / ${VAR_reqmem}" | bc)          # the number of cores for parallelizing RAM-intensive PC and K covariate caluculations is restricted by the total memory available
OUT_SUBDIR=\${OUTDIR}/\${OUTPREFIX}*/                       # output directory containing the parsed genotype and phenotype data
POP_ID_txt=\${OUT_SUBDIR}/POP_ID.txt                       # full path and file name of the text file containing the population ID (strings: e.g. 1 or 01 or 16) of each individual across populations - one ID per line
GENOME_SPEC_FNAME=\${OUTDIR}/Lperenne_genome.spec          # full path and file name of the genome specification file
QTL_SPEC_FNAME=\${OUTDIR}/QTL_SPEC.csv                     # full path and file name of the QTL specification file
NPSTAT_DIR=${DIR}/Softwares/npstat                        # full path to the npstat directory
N=${VAR_maxlib}                                           # maximum number of simulated libraries sequenced (e.g. N=500)

echo '#######################################################'
echo 'SIMULATE'
echo '#######################################################'
### Add module loading on top of the simulation script and write into a temporary script
temp_SCRIPT=SIMULATE_${VAR_prefix}_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
head -n1 \${GEN_PRED_SRC_DIR}/GPASim_01_simulate.sh > \${temp_SCRIPT}
echo -e 'module load R/3.5.2-GCC-6.2.0' >> \${temp_SCRIPT}
tail -n+2 \${GEN_PRED_SRC_DIR}/GPASim_01_simulate.sh >> \${temp_SCRIPT}
chmod +x \${temp_SCRIPT}
### Execute simulation script
time ./\${temp_SCRIPT} \
      \$QUANTINEMO_DIR \
      \$GEN_PRED_SRC_DIR \
      \$OUTDIR \
      \$OUTPREFIX \
      \$nIndividuals \
      \$nLoci \
      \$nQTL \
      \$nBGS \
      \$nAlleles \
      \$allele_eff_model \
      \$nGen \
      \$nPop \
      \$migration \
      \$selection \
      \$bg_selection \
      \$GRADIENT

echo '#######################################################'
echo 'PARSE'
echo '#######################################################'
### Add module loading on top of the simulation script and write into a temporary script
temp_SCRIPT=PARSE_${VAR_prefix}_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
head -n1 \${GEN_PRED_SRC_DIR}/GPASim_02_parse.sh > \${temp_SCRIPT}
echo -e 'module load parallel/20181222-spartan_gcc-6.2.0.lua' >> \${temp_SCRIPT}
echo -e 'module load Julia/1.1.1-spartan_gcc-6.2.0.lua' >> \${temp_SCRIPT}
echo -e 'module load R/3.5.2-GCC-6.2.0' >> \${temp_SCRIPT}
tail -n+2 \${GEN_PRED_SRC_DIR}/GPASim_02_parse.sh >> \${temp_SCRIPT}
chmod +x \${temp_SCRIPT}
### Execute parsing script
time ./\${temp_SCRIPT} \
      \${OUTDIR} \
      \${fname_dat} \
      \${nPools} \
      \${GEN_PRED_SRC_DIR} \
      \${INI_FNAME} \
      \${nCores}

echo '#######################################################'
echo 'SUMMSTATS'
echo '#######################################################'
### Add module loading on top of the simulation script and write into a temporary script
temp_SCRIPT=SUMMSTATS_${VAR_prefix}_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
head -n1 \${GEN_PRED_SRC_DIR}/GPASim_03_summaryStats.sh > \${temp_SCRIPT}
echo -e 'module load parallel/20181222-spartan_gcc-6.2.0.lua' >> \${temp_SCRIPT}
echo -e 'module load Julia/1.1.1-spartan_gcc-6.2.0.lua' >> \${temp_SCRIPT}
echo -e 'module load R/3.5.2-GCC-6.2.0' >> \${temp_SCRIPT}
tail -n+2 \${GEN_PRED_SRC_DIR}/GPASim_03_summaryStats.sh >> \${temp_SCRIPT}
chmod +x \${temp_SCRIPT}
### Execute summary statistics script
time ./\${temp_SCRIPT} \
      \${OUT_SUBDIR} \
      \${OUTPREFIX} \
      \${POP_ID_txt} \
      \${nGen} \
      \${GENOME_SPEC_FNAME} \
      \${QTL_SPEC_FNAME} \
      \${GEN_PRED_SRC_DIR} \
      \${NPSTAT_DIR} \
      \${N} \
      \${nCores}

echo '#######################################################'
echo 'GPAS'
echo '#######################################################'
### Add module loading on top of the simulation script and write into a temporary script
temp_SCRIPT=GPAS_${VAR_prefix}_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
head -n1 \${GEN_PRED_SRC_DIR}/GPASim_04_GWAS_GP.sh > \${temp_SCRIPT}
echo -e 'module load parallel/20181222-spartan_gcc-6.2.0.lua' >> \${temp_SCRIPT}
echo -e 'module load Julia/1.1.1-spartan_gcc-6.2.0.lua' >> \${temp_SCRIPT}
echo -e 'module load R/3.5.2-GCC-6.2.0' >> \${temp_SCRIPT}
tail -n+2 \${GEN_PRED_SRC_DIR}/GPASim_04_GWAS_GP.sh >> \${temp_SCRIPT}
chmod +x \${temp_SCRIPT}
### Execute genomic prediction and association analysis cross-validation script
time ./\${temp_SCRIPT} \
      \${OUT_SUBDIR} \
      \${QTL_SPEC_FNAME} \
      \${GEN_PRED_SRC_DIR} \
      \${nCores}

" >> ${VAR_prefix}_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.slurm

#########################
### SUBMIT INTO QUEUE ###
#########################
sbatch ${VAR_prefix}_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.slurm
