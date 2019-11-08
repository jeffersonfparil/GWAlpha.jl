### start an interactive session
sinteractive --account=punim0594 --partition=mig --threads=13 --time=60

# The modules to load:
module load parallel/20181222-spartan_gcc-6.2.0.lua
module load R/3.5.2-GCC-6.2.0
module load Julia/1.1.1-spartan_gcc-6.2.0.lua

# Parameter definitons
nGen=500
nQTL=10
migration=0.001
selection=0.25
bg_selection=-0.25
OUTPREFIX=LOLIUM_1rep_10QTL_0.001mr_0.25fgs_-0.25bgs            # prefix for the quantiNemo2 initiation (*.ini) file and the output folder
OUTDIR=/data/cephfs/punim0543/jparil/Output_1rep_10QTL_0.001mr_0.25fgs_-0.25bgs/  # existing output directory to dump all the output files and folder
fname_dat=${OUTPREFIX}_g${nGen}.dat
nPools=5
GEN_PRED_SRC_DIR=/data/cephfs/punim0543/jparil/Softwares/genomic_prediction/src
JULIA_DIR=/data/cephfs/punim0543/jparil/Softwares/julia-1.0.5/bin
JULIA_DEPOT=/data/cephfs/punim0543/jparil/Softwares/julia_packages
INI_FNAME=${OUTPREFIX}.ini
nCores=16

# Execute parsing
head -n1 ${GEN_PRED_SRC_DIR}/GPWASim_02_parse.sh > test.sh
echo -e "module load parallel/20181222-spartan_gcc-6.2.0.lua" >> test.sh
echo -e "module load R/3.5.2-GCC-6.2.0" >> test.sh
echo -e "module load Julia/1.1.1-spartan_gcc-6.2.0.lua" >> test.sh
tail -n+2 ${GEN_PRED_SRC_DIR}/GPWASim_02_parse.sh >> test.sh
time ${GEN_PRED_SRC_DIR}/GPWASim_02_parse.sh ${OUTDIR} ${fname_dat} ${nPools} ${JULIA_DIR} ${JULIA_DEPOT} ${GEN_PRED_SRC_DIR} ${INI_FNAME} ${nCores}
