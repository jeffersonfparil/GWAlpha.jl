##############################
### ANALYSIS PER LANDSCAPE ### Shell executioner for modules 1, 2 and 3
##############################

### INPUT
SRC_DIR=${1} ### SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src
DIR=${2} ### DIR=$(pwd)
PLOT=TRUE
bname=$(basename $DIR)
nQTL=$(echo ${bname%qtl_*} | rev | cut -d_ -f1)
nREP=${3} ### nREP=100

### SAMPLE EXECUTION
# SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src
# time \
# parallel -j 10 ${SRC_DIR}/GPASim_06_analysis.sh ${SRC_DIR} {} 100 ::: $(find /data/Lolium/Quantitative_Genetics/LOLSIM_2019_TEST/ -type d -name "LOLSIM_*rep*_*qtl*_*mr*_*fgs*_*bgs*_*grad")

### SET WQRKING DIRECTORY
cd $DIR ### need to move to the directory for @distributed to load the neccessary files across all cores

### EXECUTE
# (1) prelim stats
time Rscript  ${SRC_DIR}/GPASim_06_prelim.r $DIR $PLOT
# (2) prep for ABC
echo $nQTL > nQTL.temp
time julia    ${SRC_DIR}/GPASim_06_ABC.jl $DIR $nQTL $nREP
rm nQTL.temp
# (3) ABC resource optim
time Rscript  ${SRC_DIR}/GPASim_06_ABC.r $DIR $PLOT

### CLEANUP
mkdir PRELIM_STATS
mv MEAN_COMPARISON-*.csv PRELIM_STATS
mv $(basename $DIR)-METRICS* PRELIM_STATS
mv $(basename $DIR)-BOXPLOTS* PRELIM_STATS
mv $(basename $DIR)-HSD* PRELIM_STATS
mkdir ABC_OPTIM
mv ACROSS_*.csv ABC_OPTIM
mv WITHIN_*.csv ABC_OPTIM
mv *ABC*.svg ABC_OPTIM
mv ABC_SAMPLING_OUTPUT_100REPS.csv ABC_OPTIM

### MAIN OUTPUT
# (1) CV_OUTPUT_SUMMARY-${REP}rep-${QTL}qtl-${MGR}mgr-${FGS}fgs-${BGS}bgs-${GRAD}grad.csv
# (2) ABC_RESOURCE_OPTIM-${REP}rep-${QTL}qtl-${MGR}mgr-${FGS}fgs-${BGS}bgs-${GRAD}grad.csv
