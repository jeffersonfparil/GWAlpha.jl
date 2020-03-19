#!/bin/bash

####################################
### Parallel execution of npstat ###
#################################### NOTE: chromosome-wise parallellization

#############
### INPUT ###
#############
chrom=$1
pileup_fname=$2
n=$3 #population size
window_size=$4
NPSTAT_DIR=$5
# chrom=1
# pileup_fname=test_1kloci_g1000_p01_POOL1.pileup
# n=91.0
# window_size=10000000
# NPSTAT_DIR=/data/Lolium/Softwares/npstat

##############
### OUTPUT ###
##############
# ${pileup_fname%.pileup*}_CHR${chrom}_NPSTAT.stats (NO HEADER but each column is supposed to be:
# chrom\twindow\tcoverage\toutgroup_coverage\tdepth\tS\tWatterson_estimator\tTajimas_Pi\tTajimas_D\tFay_Wu_H_unorm\tFay_Wu_H_norm\tS_var\tV_theta\toutgroup_divergence\tnonsynonimous_polymorph\tsynonymous_polymorph\tnonsynonimous_divergence\tsynonymous_divergence\talpha)

########################
### SAMPLE EXECUTION ###
########################
# ./GPASim_03_NPSTAT_parallel.sh \
# 1 \
# QUANTI_g140_p1_POPULATION.pileup \
# 100 \
# /data/Lolium/Softwares/npstat

########################
### SETUP PARAMETERS ###
########################
grep ^${chrom} ${pileup_fname} > ${pileup_fname%.pileup*}_CHR${chrom}.pileup
### Set mincov to at most 50 to minimize combinatorics calculations
if [ ${n%.*} -gt 50 ]
then
  n=50
else
  n=${n%.*}
fi
### Set mincov to n and maxcov to 10*n
mincov=$n
maxcov=$(echo ${n%.*} "* 100" | bc)

######################
### EXECUTE NPSTAT ###
######################
echo "###############################"
echo "CHROM = $chrom"
${NPSTAT_DIR}/npstat  -n ${n%.*} \
                      -l ${window_size%.*} \
                      -mincov ${mincov} \
                      -maxcov ${maxcov} \
                      ${pileup_fname%.pileup*}_CHR${chrom}.pileup
tail -n+2 ${pileup_fname%.pileup*}_CHR${chrom}.pileup.stats > ${pileup_fname%.pileup*}_CHR${chrom}_NPSTAT.no_header ### remove header
for i in $(seq 1 $(cat ${pileup_fname%.pileup*}_CHR${chrom}_NPSTAT.no_header | wc -l)); do echo -e "$chrom" >> ${pileup_fname%.pileup*}_CHR${chrom}_NPSTAT_COL1_POOL_ID.temp; done
paste ${pileup_fname%.pileup*}_CHR${chrom}_NPSTAT_COL1_POOL_ID.temp ${pileup_fname%.pileup*}_CHR${chrom}_NPSTAT.no_header > ${pileup_fname%.pileup*}_CHR${chrom}_NPSTAT.stats ### add chrom ID
rm ${pileup_fname%.pileup*}_CHR${chrom}.pileup ${pileup_fname%.pileup*}_CHR${chrom}.pileup.stats ${pileup_fname%.pileup*}_CHR${chrom}_NPSTAT.no_header ${pileup_fname%.pileup*}_CHR${chrom}_NPSTAT_COL1_POOL_ID.temp
echo "#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#"
