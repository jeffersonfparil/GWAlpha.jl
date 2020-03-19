##############################
### ANALYSIS PER LANDSCAPE ### Shell executioner for modules 1, 2 and 3
##############################

### INPUT
SRC_DIR=${1} ### SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src ### SRC_DIR=/data/cephfs/punim0543/jparil/Softwares/genomic_prediction/src
DIR=${2} ### DIR="/data/Lolium/Quantitative_Genetics/LOLSIM_2019/LOLSIM_1rep_100qtl_0.001mr_0.25fgs_0.00bgs_0grad" ### DIR="/data/cephfs/punim0543/jparil/LOLSIM_1rep_100qtl_0.001mr_0.25fgs_0.00bgs_0grad"
PLOT=TRUE
bname=$(basename $DIR)
nQTL=$(echo ${bname%qtl_*} | rev | cut -d_ -f1 | rev)
nREP=${3} ### nREP=100
nLIB=${4} ### nLIB=10

### SAMPLE EXECUTION
# SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src
# time \
# parallel -j 10 ${SRC_DIR}/GPASim_06_analysis.sh ${SRC_DIR} {} 100 10 ::: $(find /data/Lolium/Quantitative_Genetics/LOLSIM_2019_TEST/ -type d -name "LOLSIM_*rep*_*qtl*_*mr*_*fgs*_*bgs*_*grad")

### SET WQRKING DIRECTORY
cd $DIR ### need to move to the directory for @distributed to load the neccessary files across all cores

# ### MERGE SIMLATION DATA
# head -n1 ${DIR}/$(ls $DIR | grep LOLSIM_ | head -n1)/CROSS_VALIDATION_OUTPUT_MERGED.csv > OUTPUT.csv
# for i in $(ls $DIR | grep LOLSIM_)
# do
#   echo $i
#   tail -n+2 ${DIR}/$(ls $DIR | grep LOLSIM_ | head -n1)/CROSS_VALIDATION_OUTPUT_MERGED.csv >> OUTPUT.csv
# done
# cut -d, -f1,2,3,4,5,6,7,9,12,16,17,18,21 OUTPUT.csv > OUTPUT_METRICS_ONLY.csv
# sed -i 's/FALSE_POSITIVE/FALSE_DISCOVERY/g' OUTPUT.csv
# sed -i 's/FALSE_POSITIVE/FALSE_DISCOVERY/g' OUTPUT_METRICS_ONLY.csv

### EXECUTE
# (1) prelim stats
echo -e "@@@@@@@@@@@@@@@@@@@@@@@"
echo -e "Preliminary Analysis"
echo -e "@@@@@@@@@@@@@@@@@@@@@@@"
time Rscript  ${SRC_DIR}/GPASim_06_prelim.r $DIR CROSS_VALIDATION_OUTPUT_MERGED.csv $PLOT
# time Rscript  ${SRC_DIR}/GPASim_06_prelim.r $DIR OUTPUT_METRICS_ONLY.csv $PLOT
# (1.5) test if we hit at least one error during simulation which means stopping analysis here
if [ $(cat ACROSS_INDI.csv | wc -l) -eq 0 ] || [ $(cat ACROSS_POOL.csv | wc -l) -eq 0 ] || [ $(cat WITHIN_INDI.csv | wc -l) -eq 0 ] || [ $(cat WITHIN_POOL.csv | wc -l) -eq 0 ]
then
  echo "ERROR! Exiting now!"
  mkdir ERRORED/
  mv *.csv ERRORED/
  mv *.svg ERRORED/
  mv ERRORED/CROSS_VALIDATION_OUTPUT_MERGED.csv .
  exit 1
fi
# (2) prep for ABC
echo -e "@@@@@@@@@@@@@@@@@@@@@@@"
echo -e "Prepare ABC input"
echo -e "@@@@@@@@@@@@@@@@@@@@@@@"
echo $nQTL > nQTL.temp
time julia    ${SRC_DIR}/GPASim_06_ABC.jl $DIR $nREP $nLIB
rm nQTL.temp
# (2.5) test if we hit a dataset with the number of training populations less than the number of libraries we need
if [ $(ls ABC_SAMPLING_OUTPUT_* | wc -l) -eq 0 ]
then
  echo "ERROR! Exiting now!"
  mkdir ERRORED/
  mv *.csv ERRORED/
  mv *.svg ERRORED/
  mv ERRORED/CROSS_VALIDATION_OUTPUT_MERGED.csv .
  exit 1
fi
# (3) ABC resource optim
echo -e "@@@@@@@@@@@@@@@@@@@@@@@"
echo -e "ABC optimization"
echo -e "@@@@@@@@@@@@@@@@@@@@@@@"
time Rscript  ${SRC_DIR}/GPASim_06_ABC.r $DIR ABC_SAMPLING_OUTPUT_${nREP}REPS.csv $nLIB $PLOT

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
mv *ABC*.png ABC_OPTIM
mv ABC_SAMPLING_OUTPUT_100REPS.csv ABC_OPTIM

### MAIN OUTPUT
# (1) CV_OUTPUT_SUMMARY-${REP}rep-${QTL}qtl-${MGR}mgr-${FGS}fgs-${BGS}bgs-${GRAD}grad.csv
# (2) ABC_RESOURCE_OPTIM-${REP}rep-${QTL}qtl-${MGR}mgr-${FGS}fgs-${BGS}bgs-${GRAD}grad.csv
