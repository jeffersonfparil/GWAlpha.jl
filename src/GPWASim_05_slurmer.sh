#!/bin/bash
DIR=${1}
VAR_nQTL=${2}
VAR_migration=${3}
VAR_foreground_selection=${4}
VAR_background_selection=${5}
VAR_gradient=${6}
VAR_rep=${7}
VAR_ncores=${8}
# DIR=/data/cephfs/punim0543/jparil
# VAR_nQTL=10                     # 5, 10, 100
# VAR_migration=0.001             # 0.01, 0.001, 0.0001
# VAR_foreground_selection=0.25   # 0.0, 0.25
# VAR_background_selection=-0.25  # 0.0, 0.25, -0.25
# VAR_gradient=0                  # 0, 1, 2
# VAR_rep=1                       # 1, 2, 3, 4, 5
# VAR_ncores=32                   # 32

########################
### SAMPLE EXECUTION ###
# ./GPWASim_05_slurmer.sh /data/cephfs/punim0543/jparil 10 0.001 0.25 -0.25 1 32
########################

echo -e '#!/bin/bash' > LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.slurm
echo -e "
# Created by the University of Melbourne job script generator for SLURM
# Tue Aug 13 2019 18:03:03 GMT+1000 (Australian Eastern Standard Time)

# Partition for the job:
#SBATCH --partition=physical
# #SBATCH --partition=mig

# # Use the high memory nodes only
# #SBATCH --constraint=physg4

# Multithreaded (SMP) job: must run on one node and the cloud partition
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name=\"LOLSIM\"

# The project ID which this job should run under:
#SBATCH --account=\"punim0543\"
# #SBATCH --account=\"punim0594\"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${VAR_ncores}

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=95GB

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
#SBATCH --mail-user=jparil@student.unimelb.edu.au

# Run the job from the directory where it was launched (default)

# ################
# ### SIMULATE ###
# ################
echo '#######################################################'
echo SIMULATE
echo '#######################################################'
# The modules to load:
module load R/3.5.2-GCC-6.2.0

# The job command(s):
DIR=/data/cephfs/punim0543/jparil
QUANTINEMO_DIR=\${DIR}/Softwares/quantinemo_linux                                                       # location of the quantinemo executable
GEN_PRED_SRC_DIR=\${DIR}/Softwares/genomic_prediction/src                                               # full path to the genomic_prediction/src directory
nIndividuals=1000                                                                                       # number of individuals to simulate
nLoci=10000                                                                                             # number of loci to simulate (neutral + QTL)
nQTL=$VAR_nQTL                                                                                          # number of QTL among the all the loci simulated
nBGS=100                                                                                                # number of background selection loci (should be highly polygenic >= 100)
nAlleles=5                                                                                              # number of alleles per loci, e.g. 5 for A,T,C,G, and DEL
allele_eff_model=CHISQ                                                                                  # probability distribution to sample the allele effects from (i.e. CHISQ, NORM, & UNIF)
nGen=500                                                                                                # number of generations to simulate
nPop=100                                                                                                # number of populations/subpopulations to simulate (NOTE: must have a natural number square-root)
migration=$VAR_migration                                                                                # migration rate across the populations (surrently using the 1D stepping stone model see line 117)
selection=$VAR_foreground_selection                                                                     # selection intensity (trait of interest) defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
bg_selection=$VAR_background_selection                                                                  # background selection intensity defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
GRADIENT=$VAR_gradient
OUTPREFIX=LOLIUM_${VAR_rep}rep_\${nQTL}QTL_\${migration}mr_\${selection}fgs_\${bg_selection}bgs            # prefix for the quantiNemo2 initiation (*.ini) file and the output folder
OUTDIR=\${DIR}/Output_${VAR_rep}rep_\${nQTL}QTL_\${migration}mr_\${selection}fgs_\${bg_selection}bgs/  # existing output directory to dump all the output files and folder
mkdir \${OUTDIR}                                                                                        # create the output directory

# Simulate populations
time \
\${GEN_PRED_SRC_DIR}/GPWASim_01_simulate.sh \
      \${QUANTINEMO_DIR} \
      \${GEN_PRED_SRC_DIR} \
      \${OUTDIR} \
      \${OUTPREFIX} \
      \${nIndividuals} \
      \${nLoci} \
      \${nQTL} \
      \${nBGS} \
      \${nAlleles} \
      \${allele_eff_model} \
      \${nGen} \
      \${nPop} \
      \${migration} \
      \${selection} \
      \${bg_selection} \
      \${GRADIENT}

# #############
# ### PARSE ###
# #############
echo '#######################################################'
echo PARSE
echo '#######################################################'
# The modules to load:
module load parallel/20181222-spartan_gcc-6.2.0.lua
module load R/3.5.2-GCC-6.2.0
module load Julia/1.1.1-spartan_gcc-6.2.0.lua

# Parameter definitons
nGen=500
nQTL=$VAR_nQTL
migration=$VAR_migration
selection=$VAR_foreground_selection
bg_selection=$VAR_background_selection
OUTPREFIX=LOLIUM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad            # prefix for the quantiNemo2 initiation (*.ini) file and the output folder
OUTDIR=${DIR}/Output_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad/  # existing output directory to dump all the output files and folder
fname_dat=\${OUTPREFIX}_g\${nGen}.dat
nPools=5
GEN_PRED_SRC_DIR=${DIR}/Softwares/genomic_prediction/src
INI_FNAME=\${OUTPREFIX}.ini
nCores=$(echo "(${VAR_ncores} / 3) + 3" | bc) #(1/3 + 3) the number of parallel threads for RAM-intensive parsing

# Add module loading
head -n1 \${GEN_PRED_SRC_DIR}/GPWASim_02_parse.sh > temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
echo -e 'module load parallel/20181222-spartan_gcc-6.2.0.lua' >> temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
echo -e 'module load Julia/1.1.1-spartan_gcc-6.2.0.lua' >> temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
echo -e 'module load R/3.5.2-GCC-6.2.0' >> temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
tail -n+2 \${GEN_PRED_SRC_DIR}/GPWASim_02_parse.sh >> temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh

chmod +x temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh

# Execute parsing
time \
./temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh \
\${OUTDIR} \
\${fname_dat} \
\${nPools} \
\${GEN_PRED_SRC_DIR} \
\${INI_FNAME} \
\${nCores}

# #################
# ### SUMMSTATS ###
# #################
echo '#######################################################'
echo SUMMSTATS
echo '#######################################################'
# The modules to load:
module load parallel/20181222-spartan_gcc-6.2.0.lua
module load R/3.5.2-GCC-6.2.0
module load Julia/1.1.1-spartan_gcc-6.2.0.lua
module load GSL/2.5-intel-2018.u4

# Parameter definitons
nGen=500
nQTL=$VAR_nQTL
migration=$VAR_migration
selection=$VAR_foreground_selection
bg_selection=$VAR_background_selection
OUTPREFIX=LOLIUM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad
BASEDIR=${DIR}/Output_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad/
OUTDIR=\${BASEDIR}/\${OUTPREFIX}*/
POP_ID_txt=\${OUTDIR}/POP_ID.txt
GENOME_SPEC_FNAME=\${BASEDIR}/Lperenne_genome.spec
QTL_SPEC_FNAME=\${BASEDIR}/QTL_SPEC.csv
GEN_PRED_SRC_DIR=${DIR}/Softwares/genomic_prediction/src
NPSTAT_DIR=${DIR}/Softwares/npstat
N=500 #maximum number of simulated libraries sequenced (e.g. N=500)
nCores=$(echo "${VAR_ncores} / 4" | bc) #(1/4) the number of parallel threads for RAM-intensive PC and K covariate caluculations

# Add module loading
head -n1 \${GEN_PRED_SRC_DIR}/GPWASim_03_summaryStats.sh > temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
echo -e 'module load parallel/20181222-spartan_gcc-6.2.0.lua' >> temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
echo -e 'module load Julia/1.1.1-spartan_gcc-6.2.0.lua' >> temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
echo -e 'module load R/3.5.2-GCC-6.2.0' >> temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
tail -n+2 \${GEN_PRED_SRC_DIR}/GPWASim_03_summaryStats.sh >> temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh

chmod +x temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh

# Execute summstats
time \
./temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh \
  \${OUTDIR} \
  \${OUTPREFIX} \
  \${POP_ID_txt} \
  \${nGen} \
  \${GENOME_SPEC_FNAME} \
  \${QTL_SPEC_FNAME} \
  \${GEN_PRED_SRC_DIR} \
  \${NPSTAT_DIR} \
  \${N} \
  \${nCores}

# #############
# ### GPWAS ###
# #############
echo '#######################################################'
echo GPWAS
echo '#######################################################'
# The modules to load:
module load parallel/20181222-spartan_gcc-6.2.0.lua
module load R/3.5.2-GCC-6.2.0
module load Python/2.7.13-GCC-6.2.0-bare
module load Julia/1.1.1-spartan_gcc-6.2.0.lua
module load GSL/2.5-intel-2018.u4

# Parameter definitons
nGen=500
nQTL=$VAR_nQTL
migration=$VAR_migration
selection=$VAR_foreground_selection
bg_selection=$VAR_background_selection
OUTPREFIX=LOLIUM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad
BASEDIR=${DIR}/Output_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad/
OUTDIR=\${BASEDIR}/\${OUTPREFIX}*/
QTL_SPEC_FNAME=\${BASEDIR}/QTL_SPEC.csv
GEN_PRED_SRC_DIR=${DIR}/Softwares/genomic_prediction/src

# Add module loading
head -n1 \${GEN_PRED_SRC_DIR}/GPWASim_04_GWAS_GP.sh > temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
echo -e 'module load parallel/20181222-spartan_gcc-6.2.0.lua' >> temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
echo -e 'module load R/3.5.2-GCC-6.2.0' >> temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
echo -e 'module load Python/2.7.13-GCC-6.2.0-bare' >> temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
echo -e 'module load Julia/1.1.1-spartan_gcc-6.2.0.lua' >> temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh
tail -n+2 \${GEN_PRED_SRC_DIR}/GPWASim_04_GWAS_GP.sh >> temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh

chmod +x temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh

# Execute GPWAS
time \
./temp_LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.sh \
  \${OUTDIR} \
  \${QTL_SPEC_FNAME} \
  \${GEN_PRED_SRC_DIR}

" >> LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.slurm
sbatch LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.slurm
