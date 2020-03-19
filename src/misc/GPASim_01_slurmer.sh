#!/bin/bash
VAR_nQTL=$1
VAR_migration=$2
VAR_foreground_selection=$3
VAR_background_selection=$4
VAR_gradient=$5
VAR_rep=$6
# VAR_nQTL=5                  # 5, 10, 100
# VAR_migration=0.001         # 0.01, 0.001, 0.0001
# VAR_foreground_selection=1  # 1, 0
# VAR_background_selection=0  # 0, 1
# VAR_gradient=0              # 0, 1, 2
# VAR_rep=1                   # 1, 2, 3, 4, 5

########################
### SAMPLE EXECUTION ###
# ./GPASim_01_slurmer.sh 10 0.001 0.25 -0.25 0 1
# or
# for nqtl in 5 10 100; do
#   for mr in 0.01 0.001 0.0001; do
#     for fgs in 0.0 0.25; do
#       for bgs in 0.0 0.25 -0.25; do
#         for grad in 0 1 2; do
#            for rep in 1 2 3 4 5; dor
#             ./GPASim_01_slurmer.sh ${nqtl} ${mr} ${fgs} ${bgs} ${grad} ${rep}
#            done
#         done
#       done
#     done
#   done
# done
########################

echo -e '#!/bin/bash' > LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.slurm
echo -e "
# Created by the University of Melbourne job script generator for SLURM
# Tue Aug 13 2019 18:03:03 GMT+1000 (Australian Eastern Standard Time)

# Partition for the job:
#SBATCH --partition=physical

# # Use the high memory nodes only
# #SBATCH --constraint=physg4

# Multithreaded (SMP) job: must run on one node and the cloud partition
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name=\"LOLSIM\"

# The project ID which this job should run under:
#SBATCH --account=\"punim0543\"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=95GB

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-8:0:00

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

# The modules to load:
module load R/3.5.2-GCC-6.2.0

# The job command(s):
DIR=/data/cephfs/punim0543/jparil
QUANTINEMO_DIR=\${DIR}/Softwares/quantinemo_linux                                                                       # location of the quantinemo executable
GEN_PRED_SRC_DIR=\${DIR}/Softwares/genomic_prediction/src                                                               # full path to the genomic_prediction/src directory
nIndividuals=1000                                                                                                       # number of individuals to simulate
nLoci=10000                                                                                                             # number of loci to simulate (neutral + QTL)
nQTL=$VAR_nQTL                                                                                                          # number of QTL among the all the loci simulated
nBGS=100                                                                                                                # number of background selection loci (should be highly polygenic >= 100)
nAlleles=5                                                                                                              # number of alleles per loci, e.g. 5 for A,T,C,G, and DEL
allele_eff_model=CHISQ                                                                                                  # probability distribution to sample the allele effects from (i.e. CHISQ, NORM, & UNIF)
nGen=500                                                                                                                # number of generations to simulate
nPop=100                                                                                                                # number of populations/subpopulations to simulate (NOTE: must have a natural number square-root)
migration=$VAR_migration                                                                                                # migration rate across the populations (surrently using the 1D stepping stone model see line 117)
selection=$VAR_foreground_selection                                                                                     # selection intensity (trait of interest) defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
bg_selection=$VAR_background_selection                                                                                  # background selection intensity defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
GRADIENT=$VAR_gradient                                                                                                  # gradient of diffusing non-wildtype alleles
OUTPREFIX=LOLIUM_${VAR_rep}rep_\${nQTL}QTL_\${migration}mr_\${selection}fgs_\${bg_selection}bgs_\${GRADIENT}grad        # prefix for the quantiNemo2 initiation (*.ini) file and the output folder
OUTDIR=\${DIR}/Output_${VAR_rep}rep_\${nQTL}QTL_\${migration}mr_\${selection}fgs_\${bg_selection}bgs_\${GRADIENT}grad/  # existing output directory to dump all the output files and folder
mkdir \${OUTDIR}                                                                                                        # create the output directory

# Simulate populations
time \
\${GEN_PRED_SRC_DIR}/GPASim_01_simulate.sh \
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
" >> LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.slurm
sbatch LOLSIM_${VAR_rep}rep_${VAR_nQTL}QTL_${VAR_migration}mr_${VAR_foreground_selection}fgs_${VAR_background_selection}bgs_${VAR_gradient}grad.slurm
