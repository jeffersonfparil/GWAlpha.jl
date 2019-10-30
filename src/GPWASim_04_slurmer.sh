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
# TEST WITH sinteractive
# sinteractive -T13 -t120
# EXAMPLE 1:
# ./GPWASim_04_slurmer.sh /data/cephfs/punim0543/jparil 10 0.001 0.25 -0.25 1 32
# EXAMPLE 2:
# for nqtl in 5 10 100; do
#   for mr in 0.01 0.001 0.0001; do
#     for fgs in 0.0 0.25; do
#       for bgs in 0.0 0.25 -0.25; do
#         for grad in 0 1 2; do
#            for rep in 1 2 3 4 5; dor
#             ./GPWASim_04_slurmer.sh /data/cephfs/punim0543/jparil ${nqtl} ${mr} ${fgs} ${bgs} ${grad} ${rep} 32
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
#SBATCH --mem=100GB

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=2-0:0:00

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
