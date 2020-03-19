#!/bin/bash

##############
### INPUTS ###
##############
GEN_PRED_SRC_DIR=${1}
DIR=${2}
VAR_partition=${3}

# GEN_PRED_SRC_DIR=/data/cephfs/punim0543/jparil/Softwares/genomic_prediction/src
# DIR=/data/cephfs/punim0543/jparil/LOLSIM_5rep_5qtl_0.001mr_0.25fgs_0.25bgs_2grad
# VAR_partition=snowy               # cloud, physical, snowy, ..., etc

########################
### SAMPLE EXECUTION ###
### (1) single execution
# cd /data/cephfs/punim0543/jparil/Softwares/genomic_prediction/src
# ./GPASim_06_slurmer.sh \
#                         /data/cephfs/punim0543/jparil/Softwares/genomic_prediction/src \
#                         /data/cephfs/punim0543/jparil/LOLSIM_5rep_5qtl_0.001mr_0.25fgs_0.25bgs_2grad \
#                         snowy
### (2) across all landscapes
# cd /data/cephfs/punim0543/jparil/Softwares/genomic_prediction/src
# for DIR in $(ls /data/cephfs/punim0543/jparil/ | grep LOLSIM_)
# do
#   echo $DIR
#   ./GPASim_06_slurmer.sh \
#         /data/cephfs/punim0543/jparil/Softwares/genomic_prediction/src \
#         /data/cephfs/punim0543/jparil/${DIR} \
#         snowy
# done

########################

########################################
### AUXILLARY NON-USER-FACING INPUTS ###
########################################
VAR_account=punim0543 ### Spartan account name
VAR_timelim="0-08:00:00" ### maximum computing time limit in days-hours:minutes:seconds

if [ $VAR_partition == cloud ]
then
  VAR_ncores=12 ### number of cores avaible given the partition
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

SLURM_NAME=$(basename $DIR)
echo -e '#!/bin/bash' > ${SLURM_NAME}.slurm
echo -e "
# Partition for the job:
#SBATCH --partition=${VAR_partition}
# Multithreaded (SMP) job: must run on one node and the cloud partition
#SBATCH --nodes=1
# The name of the job:
#SBATCH --job-name=${DIR}
# The project ID which this job should run under:
#SBATCH --account=${VAR_account}
# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${VAR_ncores}
# The amount of memory in megabytes per process in the job:
#SBATCH --mem=${VAR_mem}GB
# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=${VAR_timelim}
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
echo 'ANALYSIS: ${DIR}'
echo '#######################################################'
### Add module loading on top of the simulation script and write into a temporary script
temp_SCRIPT=${SLURM_NAME}.sh
head -n1 ${GEN_PRED_SRC_DIR}/GPASim_06_analysis.sh > \${temp_SCRIPT}
echo -e 'module load parallel/20181222-spartan_gcc-6.2.0.lua' >> \${temp_SCRIPT}
echo -e 'module load Julia/1.1.1-spartan_gcc-6.2.0.lua' >> \${temp_SCRIPT}
echo -e 'module load R/3.5.2-GCC-6.2.0' >> \${temp_SCRIPT}
echo -e 'module load GSL/2.5-intel-2018.u4' >> \${temp_SCRIPT}
tail -n+2 ${GEN_PRED_SRC_DIR}/GPASim_06_analysis.sh >> \${temp_SCRIPT}
chmod +x \${temp_SCRIPT}
### Execute analysis script
time ./\${temp_SCRIPT} \
      ${GEN_PRED_SRC_DIR} \
      ${DIR} \
      100 \
      10
rm \${temp_SCRIPT}

echo '#######################################################'
echo 'ANALYSIS FINISHED FOR ${DIR} ::: ' $(date)
echo '#######################################################'
" >> ${SLURM_NAME}.slurm

#########################
### SUBMIT INTO QUEUE ###
#########################
sbatch ${SLURM_NAME}.slurm

# ##################
# ### MONITORING ###
# ##################
# VAR_user="jparil"
# squeue -u ${VAR_user}
# ### monitor via htop in the node (testing: considering only the first job)
# NODE_ID=$(squeue -u ${VAR_user} | grep R | tail -n+2 | head -n1 | sed -z "s/ \+/\t/g" | cut -f2- | rev | cut -f1 | rev)
# ssh $NODE_ID
# htop
# logout
# ### ONE-LINER:
# squeue -u jparil; echo RUNNING: $(echo $(squeue -u jparil | grep "R" | wc -l) - 1 | bc); echo PENDING: $(squeue -u jparil | grep -v "R" | wc -l); echo FINISHED: $(grep "FINISHED" slurm-* | wc -l); echo ERRORED: $(grep "Exiting now" slurm-* | wc -l)

# ###########################
# ### INTERACTIVE SESSION ###
# ###########################
# sinteractive --partition=snowy -c 32 --time=0-2:0:0
# module load parallel/20181222-spartan_gcc-6.2.0.lua
# module load Julia/1.1.1-spartan_gcc-6.2.0.lua
# module load R/3.5.2-GCC-6.2.0
# module load GSL/2.5-intel-2018.u4
