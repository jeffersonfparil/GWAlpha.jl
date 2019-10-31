#!/bin/bash
#######################################################################################
###                                                                                 ###
### Genomic prediction modeling and association analyses using Pool-sequencing data ###
###                                   (PoolGPAS)                                    ###
###                                                                                 ###
#######################################################################################

###################
### USER INPUTS ###
###################
DIR=${1}              ### working directory where the fastq files are located
RGRP=${2}             ### text file indicating the identiying labels for read1 and read2 (e.g. line1:_R1 and line2:_R2) ### NOTE: It's nice to ad whatever delimiter/s flanking these identifiers
PGRP=${3}             ### text file indicating the identiying labels for pool groups (e.g. line1:IF, line2:IG, ...)
MAPQ=${4}             ### minimum mapping quality threshold
BASQ=${5}             ### minimum base quality threshold
REFLINK=${6}          ### download link to the reference genome in fasta format
POPOOLATION_DIR=${7}  ### directory of the popoolation2 software
SRC_DIR=${8}          ### location of genomic_prediction/src
MAF=${9}              ### miimum allele frequency threshold
DEPTH=${10}           ### minimum depth threshold
POOLSIZES=${11}       ### pool sizes per pool grouping: e.g. IT\t100,100,100,100,100\nUG\t100,100,100,100,100\n
# ### TEST:
# DIR=/data/Lolium/Quantitative_Genetics/GWAS_GP_2018_Inverleigh_Urana
# echo -e "_R1" > ${DIR}/read_grouping.txt
# echo -e "_R2" >> ${DIR}/read_grouping.txt
# RGRP=${DIR}/read_grouping.txt
# echo -e "IF" > ${DIR}/pool_grouping.txt
# echo -e "IG" >> ${DIR}/pool_grouping.txt
# echo -e "IS" >> ${DIR}/pool_grouping.txt
# echo -e "IT" >> ${DIR}/pool_grouping.txt
# echo -e "UF" >> ${DIR}/pool_grouping.txt
# echo -e "UG" >> ${DIR}/pool_grouping.txt
# echo -e "US" >> ${DIR}/pool_grouping.txt
# echo -e "UT" >> ${DIR}/pool_grouping.txt
# PGRP=${DIR}/pool_grouping.txt
# MAPQ=20 # minimum mapping quality at 1% error rate
# BASQ=20 # minimum base calling quality at 1% error rate
# REFLINK=http://185.45.23.197:5080/ryegrassdata/GENOMIC_ASSEMBLY/lope_V1.0.fasta
# POPOOLATION_DIR=/data/Lolium/Softwares/popoolation2_1201
# SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src
# MAF=0.001
# DEPTH=10
# echo -e "IF\t74,74,74,74,74" > ${DIR}/pool_sizes_per_grouping.txt
# echo -e "IG\t84,84,85,84,84" >> ${DIR}/pool_sizes_per_grouping.txt
# echo -e "IS\t92,92,95,92,92" >> ${DIR}/pool_sizes_per_grouping.txt
# echo -e "IT\t89,89,91,89,89" >> ${DIR}/pool_sizes_per_grouping.txt
# echo -e "UF\t66,66,67,66,66" >> ${DIR}/pool_sizes_per_grouping.txt
# echo -e "UG\t71,71,71,71,71" >> ${DIR}/pool_sizes_per_grouping.txt
# echo -e "US\t80,80,84,80,80" >> ${DIR}/pool_sizes_per_grouping.txt
# echo -e "UT\t80,80,84,80,80" >> ${DIR}/pool_sizes_per_grouping.txt
# POOLSIZES=${DIR}/pool_sizes_per_grouping.txt

#####################################
### (0) SETUP WORKING DIRECTORIES ###
#####################################
echo -e "################################################"
echo -e "(0) SETUP WORKING DIRECTORIES"
echo -e "################################################"
# (0a) test if there are paired-end fastq files (qzipped or ungzipped) in the working directory
if [ $(ls ${DIR}/*.fastq* | wc -l) -lt 2 ]
then
  echo -e "The paired-end reads (*.fastq or *.fastq.gz) files not found in the directory:"
  echo -e "${DIR}"
  echo -e "Exiting now!"
  exit
fi
# (0b) set working directories
mkdir ${DIR}/FASTQ/; mv ${DIR}/*.fastq* ${DIR}/FASTQ/ ### moving fastq files into FASTQ/ directory
mkdir ${DIR}/REF/   ### dowload reference genome here
mkdir ${DIR}/SAM/   ### for storing BAM alignments
mkdir ${DIR}/VCF/   ### for storing mpileups and sync files
mkdir ${DIR}/GPAS/   ### for storing mpileups and sync files
# (0c) dowload the fasta reference genome and fix the format
if [ $(ls ${DIR}/REF/$(basename $REFLINK) | wc -l) == 0  ]
then
  cd ${DIR}/REF/
  wget ${REFLINK}
  echo -e "import os, sys" > fix_REF.py
  echo -e "from Bio import SeqIO" >> fix_REF.py
  echo -e "INPUT_fasta = sys.argv[1]" >> fix_REF.py
  echo -e "OUTPUT_fasta = sys.argv[2]" >> fix_REF.py
  echo -e "SeqIO.convert(INPUT_fasta, 'fasta', OUTPUT_fasta, 'fasta')" >> fix_REF.py
  python fix_REF.py ${DIR}/REF/$(basename $REFLINK) ${DIR}/REF/Reference.fasta || \
  echo -e "Please check your reference fasta.\nExiting now!" | exit
  bwa index -p Reference -a bwtsw Reference.fasta # index the reference genome for bwa
  samtools faidx ${DIR}/REF/Reference.fasta # index the reference genome for samtools
  # mv REF/Reference.fasta.fai REF/Reference_genome.fai
  cd -
fi

#################
### (1) ALIGN ###
#################
echo -e "################################################"
echo -e "(1) ALIGN"
echo -e "################################################"
READ1=$(head -n1 $RGRP)
READ2=$(tail -n1 $RGRP)
echo -e "bwa mem ${DIR}/REF/Reference ${DIR}/FASTQ/\${1} ${DIR}/FASTQ/\${2} | \\" > align.sh
echo -e "samtools view -q ${MAPQ} -b | \\" >> align.sh
echo -e "samtools sort > ${DIR}/SAM/\${3%.fastq*}.bam" >> align.sh
chmod +x align.sh
time \
parallel --link ./align.sh {1} {2} {3} ::: \
$(ls ${DIR}/FASTQ/ | grep "${READ1}") ::: \
$(ls ${DIR}/FASTQ/ | grep "${READ2}") ::: \
$(ls ${DIR}/FASTQ/ | grep "${READ1}" | sed s/"$READ1"/""/g)
rm align.sh

##################
### (2) PILEUP ###
##################
echo -e "################################################"
echo -e "(2) PILEUP"
echo -e "################################################"
# (2a) generate bam list
for i in $(cat $PGRP)
do
  echo $i
  ls ${DIR}/SAM/${i}*.bam > ${DIR}/bam_list_${i}.txt
done
# (2b) pileup alignments with minimum mapping and base qualities set to PHRED13 or 5% error rate
echo -e "samtools mpileup \\" > mpileup.sh
echo -e "--max-depth 10000 \\" >> mpileup.sh
echo -e "--min-MQ ${MAPQ} \\" >> mpileup.sh
echo -e "--min-BQ ${BASQ} \\" >> mpileup.sh
echo -e "--fasta-ref ${DIR}/REF/Reference.fasta \\" >> mpileup.sh
echo -e "--bam-list ${DIR}/bam_list_\${1}.txt  > ${DIR}/VCF/\${1}.mpileup" >> mpileup.sh
chmod +x mpileup.sh
time \
parallel ./mpileup.sh {} ::: $(cat $PGRP)
rm mpileup.sh

#######################
### (3) SYNCHRONIZE ###
#######################
echo -e "################################################"
echo -e "(3) SYNCHRONIZE MPILEUP"
echo -e "################################################"
# (3a) set up popoolation2
if [ $(ls ${POPOOLATION_DIR} | wc -l) == 0 ]
then
  cd ${POPOOLATION_DIR%popoolation2_1201*}
  wget https://sourceforge.net/projects/popoolation2/files/popoolation2_1201.zip
  unzip popoolation2_1201.zip -d popoolation2_1201
  rm popoolation2_1201.zip
  cd -
fi
# (3b) synchronize mpileup/s
USABLE_RAM=$(echo "($(grep MemTotal /proc/meminfo | awk '{print $2}') / 1000000) - 2" | bc)
NCORES=$(grep -c ^processor /proc/cpuinfo)
NPROCS=$(cat $PGRP | wc -l)
if [ $NCORES -gt $NPROCS ]
then
  MEMORY=$(echo $USABLE_RAM / $NPROCS | bc)
  NTHREADS=$(echo $NCORES / $NPROCS | bc)
else
  MEMORY=$(echo $USABLE_RAM / $NCORES | bc)
  NTHREADS=1
fi
time \
parallel \
java -ea -Xmx${MEMORY}g -jar ${POPOOLATION_DIR}/mpileup2sync.jar \
  --input ${DIR}/VCF/{}.mpileup \
  --output ${DIR}/VCF/{}_MAPQ${MAPQ}_BASQ${BASQ}.sync \
  --fastq-type sanger \
  --min-qual ${BASQ} \
  --threads ${NTHREADS} ::: $(cat $PGRP)

###################################################################
### (4) FILTER BY MAF AND DEPTH AND CALCULTE ALLELE FREQUENCIES ###
###################################################################
# (4a) set up genomic_prediction repo clone
if [ $(ls ${SRC_DIR} | wc -l) == 0 ]
then
  cd ${SRC_DIR%popoolation2_1201*}
  git clone https://gitlab.com/jeffersonfparil/genomic_prediction.git
  cd -
fi
# (4b) filter sync by MAF and DEPTH and calculate the alee frequencies
echo -e 'SRC_DIR = ARGS[1]' > SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e 'FNAME_SYNC = ARGS[2]' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e 'MAF = parse(Float64, ARGS[3])' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e 'DEPTH = parse(Int64, ARGS[4])' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e '# SRC_DIR = "/data/Lolium/Softwares/genomic_prediction/src"' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e '# FNAME_SYNC = "/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/test.sync"' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e '# MAF = 0.001' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e '# DEPTH = 10' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e 'push!(LOAD_PATH, SRC_DIR)' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e 'using filter_sync_module' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e 'using sync_parsing_module' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e 'filter_sync_module.filter_sync(filename_sync=FNAME_SYNC, MAF=MAF, DEPTH=DEPTH)' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e 'FILTERED_SYNC = string(join(split(FNAME_SYNC, ".")[1:(end-1)], "."), "_MAF", MAF, "_DEPTH", DEPTH, ".sync")' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e 'sync_parsing_module.sync_parse(FILTERED_SYNC)' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
echo -e '# OUTPUT: *((_MAF${MAF}_DEPTH${DEPTH}.sync))' >> SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl
USABLE_RAM=$(echo "($(grep MemTotal /proc/meminfo | awk '{print $2}') / 1000000) - 2" | bc)
SPLIT_SIZE=1000 # 1000 megabytes which translates to ~10 GB of memory usage in julia
MEMORY_PER_SPLIT=$(echo "($SPLIT_SIZE / 1000) * 10" | bc) # gigabytes: 1GB split size is to 10 GB RAM required by julia
NTHREADS=$(echo "$USABLE_RAM / $MEMORY_PER_SPLIT" | bc)
time \
for i in $(cat $PGRP)
do
  # i=IF
  split --numeric-suffixes=0000 -C ${SPLIT_SIZE}m ${DIR}/VCF/${i}_MAPQ${MAPQ}_BASQ${BASQ}.sync ${DIR}/VCF/${i}_MAPQ${MAPQ}_BASQ${BASQ}_SPLIT
  for j in $(ls ${DIR}/VCF/${i}_MAPQ${MAPQ}_BASQ${BASQ}_SPLIT*); do mv $j ${j}.sync; done ### suffix the splits with .sync
  time \
  parallel --jobs ${NTHREADS} julia SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl \
  ${SRC_DIR} \
  {1} \
  ${MAF} \
  ${DEPTH} ::: $(ls ${DIR}/VCF/${i}_MAPQ${MAPQ}_BASQ${BASQ}_SPLIT*.sync)
  cat ${DIR}/VCF/${i}_MAPQ${MAPQ}_BASQ${BASQ}_SPLIT*_MAF${MAF}_DEPTH${DEPTH}.sync > ${DIR}/VCF/${i}_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}.sync
  cat ${DIR}/VCF/${i}_MAPQ${MAPQ}_BASQ${BASQ}_SPLIT*_MAF${MAF}_DEPTH${DEPTH}_ALLELEFREQ.csv > ${DIR}/VCF/${i}_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_ALLELEFREQ.csv
  rm ${DIR}/VCF/${i}_MAPQ${MAPQ}_BASQ${BASQ}_SPLIT*
done
rm SYNC_2_SYNC_2_ALLELEFREQ_MAF_DEPTH_FILTERED.jl

####################################################################
### (5) POOL-GWAS AND POOL-GP [FIXED_{LS, RR, GLMNET and LASSO}] ###
####################################################################
# (5a) Generate covariates: PC and Fst
echo -e 'sync_fname = ARGS[1]' > FST_PAIRWISE.jl
echo -e 'window_size = parse(Int64, ARGS[2])' >> FST_PAIRWISE.jl
echo -e 'pool_sizes = parse.(Int64, split(ARGS[3], ",")) #poolsizes delimited by comma e.g. 10,20,30,20,10' >> FST_PAIRWISE.jl
echo -e 'METHOD = ARGS[4]' >> FST_PAIRWISE.jl
echo -e 'SRC_DIR = ARGS[5]' >> FST_PAIRWISE.jl
echo -e '# sync_fname = "/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/test.sync"' >> FST_PAIRWISE.jl
echo -e '# window_size = 10000' >> FST_PAIRWISE.jl
echo -e '# pool_sizes = repeat([42], inner=60)' >> FST_PAIRWISE.jl
echo -e '# METHOD = "Hivert" ### WeirCock or Hivert' >> FST_PAIRWISE.jl
echo -e '# SRC_DIR = "/data/Lolium/Softwares/genomic_prediction/src"' >> FST_PAIRWISE.jl
echo -e 'using CSV' >> FST_PAIRWISE.jl
echo -e 'using DataFrames' >> FST_PAIRWISE.jl
echo -e 'using Statistics' >> FST_PAIRWISE.jl
echo -e 'using DelimitedFiles' >> FST_PAIRWISE.jl
echo -e 'push!(LOAD_PATH, SRC_DIR)' >> FST_PAIRWISE.jl
echo -e 'using poolFST_module' >> FST_PAIRWISE.jl
echo -e 'poolFST_module.Fst_pairwise(sync_fname=sync_fname, window_size=window_size, pool_sizes=pool_sizes, METHOD=METHOD)' >> FST_PAIRWISE.jl
grep -f $PGRP $POOLSIZES | cut -f2 > GREPPED_POOLSIZES.temp
time \
parallel --link julia FST_PAIRWISE.jl \
${DIR}/VCF/{1}_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}.sync \
10000 \
{2} \
Hivert \
${SRC_DIR} ::: $(cat $PGRP) ::: $(cat GREPPED_POOLSIZES.temp)
rm FST_PAIRWISE.jl GREPPED_POOLSIZES.temp
# (5b) Pool-GWAS and Pool-GP
echo -e "### INPUTS" > POOL_GWAS_GP.jl
echo -e "SRC_DIR = ARGS[1]" >> POOL_GWAS_GP.jl
echo -e "filename_sync = ARGS[2]" >> POOL_GWAS_GP.jl
echo -e "filename_phen_py = ARGS[3]" >> POOL_GWAS_GP.jl
echo -e "filename_phen_csv = ARGS[4]" >> POOL_GWAS_GP.jl
echo -e "filename_covariate_csv = ARGS[5]" >> POOL_GWAS_GP.jl
echo -e "MAF = parse(Float64, ARGS[6])" >> POOL_GWAS_GP.jl
echo -e "DEPTH = parse(Int64, ARGS[7])" >> POOL_GWAS_GP.jl
echo -e "MODEL = ARGS[8]" >> POOL_GWAS_GP.jl
echo -e "# ### TESTS" >> POOL_GWAS_GP.jl
echo -e "# cd(\"/data/Lolium/Quantitative_Genetics/GWAS_GP_2018_Inverleigh_Urana\")" >> POOL_GWAS_GP.jl
echo -e "# SRC_DIR = \"/data/Lolium/Softwares/genomic_prediction/src\"" >> POOL_GWAS_GP.jl
echo -e "# filename_sync = \"VCF/UG_MAPQ20_BASQ20_MAF0.001_DEPTH10.sync\"" >> POOL_GWAS_GP.jl
echo -e "# filename_phen_py = \"UG_pheno.py\"" >> POOL_GWAS_GP.jl
echo -e "# filename_phen_csv = \"UG_pheno.csv\"" >> POOL_GWAS_GP.jl
echo -e "# filename_covariate_csv = 'VCF/UG_MAPQ20_BASQ20_MAF0.001_DEPTH10_COVARIATE_FST.csv'" >> POOL_GWAS_GP.jl
echo -e "# MAF = 0.001" >> POOL_GWAS_GP.jl
echo -e "# DEPTH = 10" >> POOL_GWAS_GP.jl
echo -e "# MODEL = \"FIXED_GWAlpha\" ### FIXED_GWAlpha, FIXED_LS, FIXED_RR, FIXED_GLMNET, FIXED_LASSO, MIXED_RR, MIXED_GLMNET, and MIXED_LASSO" >> POOL_GWAS_GP.jl
echo -e "### LOAD LIBRARIES" >> POOL_GWAS_GP.jl
echo -e "using DelimitedFiles" >> POOL_GWAS_GP.jl
echo -e "using DataFrames" >> POOL_GWAS_GP.jl
echo -e "using CSV" >> POOL_GWAS_GP.jl
echo -e "push!(LOAD_PATH, SRC_DIR)" >> POOL_GWAS_GP.jl
echo -e "using PoolGPAS_module" >> POOL_GWAS_GP.jl
echo -e "### LOAD COVARIATE" >> POOL_GWAS_GP.jl
echo -e "if MODEL == \"FIXED_GWAlpha\"" >> POOL_GWAS_GP.jl
echo -e "COVARIATE = nothing" >> POOL_GWAS_GP.jl
echo -e "else" >> POOL_GWAS_GP.jl >> POOL_GWAS_GP.jl
echo -e "COVARIATE = DelimitedFiles.readdlm(filename_covariate_csv, ',')" >> POOL_GWAS_GP.jl
echo -e "  end" >> POOL_GWAS_GP.jl >> POOL_GWAS_GP.jl
echo -e "### EXECUTE" >> POOL_GWAS_GP.jl >> POOL_GWAS_GP.jl
echo -e "println(\"#######################################################\")" >> POOL_GWAS_GP.jl >> POOL_GWAS_GP.jl
echo -e "println(MODEL)" >> POOL_GWAS_GP.jl >> POOL_GWAS_GP.jl
echo -e "println(\"#######################################################\")" >> POOL_GWAS_GP.jl >> POOL_GWAS_GP.jl
echo -e "if MODEL == \"FIXED_GWAlpha\"" >> POOL_GWAS_GP.jl
echo -e "  OUT, COVAR_EFF = PoolGPAS_module.PoolGPAS(filename_sync, filename_phen_py, MAF, DEPTH, MODEL=MODEL, COVARIATE=COVARIATE)" >> POOL_GWAS_GP.jl
echo -e "else" >> POOL_GWAS_GP.jl
echo -e "  OUT, COVAR_EFF = PoolGPAS_module.PoolGPAS(filename_sync, filename_phen_csv, MAF, DEPTH, MODEL=MODEL, COVARIATE=COVARIATE)" >> POOL_GWAS_GP.jl
echo -e "end" >> POOL_GWAS_GP.jl
echo -e "### OUTPUTS:" >> POOL_GWAS_GP.jl
echo -e "### (1) Allelic effects file (*-${MODEL}_Alphas.csv)" >> POOL_GWAS_GP.jl
echo -e "### (2) Manhattan plot (*-${MODEL}_Manhattan.png)" >> POOL_GWAS_GP.jl
time \
parallel julia POOL_GWAS_GP.jl \
${SRC_DIR} \
${DIR}/VCF/{1}_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}.sync \
${DIR}/{1}_pheno.py \
${DIR}/{1}_pheno.csv \
${DIR}/VCF/{1}_MAPQ${MAPQ}_BASQ${BASQ}_MAF${MAF}_DEPTH${DEPTH}_COVARIATE_FST.csv \
${MAF} \
${DEPTH} \
{2} ::: \
$(cat $PGRP) ::: \
FIXED_GWAlpha FIXED_LS FIXED_RR FIXED_GLMNET FIXED_LASSO MIXED_RR MIXED_GLMNET MIXED_LASSO
rm POOL_GWAS_GP.jl
# (5c) Clean up
mv ${DIR}/*.png ${DIR}/GPAS/
mv ${DIR}/*_Alphas.csv ${DIR}/GPAS/

#########################
### (7) PEAK ANALYSIS ###
#########################
### EXPLORATIONS in R:
ARGS = commandArgs(trailing=TRUE)
dat = read.csv(ARGS[1])
# dat = read.csv("IS_pheno-FIXED_GWAlpha_Alphas.csv")



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

# #################################
# ### SETUP WORKING DIRECTORIES ###
# #################################NOTE: awk filering of pileup file by depth does not translate to sync file - some reads seem to be lost
# DEPTH=10
# MAF=$(echo 0$(echo "scale=2; 1/(2*42)" | bc))
# DIR=/data/Lolium/Quantitative_Genetics/GPWAS_60_SEAUS
# SRC_DIR=/data/Lolium/Softwares/genomic_prediction/misc
# SYNC_FNAME=/data/Lolium/Population_Genetics/2018November_Ryegrass_Collection/VCF/Lolium2019_${DEPTH}X_30PHRED_filtered.sync
# NITERATIONS=100
# # BAMLIST_FNAME=/data/Lolium/Population_Genetics/2018November_Ryegrass_Collection/bam.list
# mkdir ${DIR}
# cd ${DIR}
# cut -f1-63 ${SYNC_FNAME} > input.sync ### only the 60 accessions (exclude the 4 resistant-individuals bulks)
#
# #########################
# ### PREPARE GENOTYPES ###
# #########################
# julia ${SRC_DIR}/parse_sync2csv_filter_MAF_DEPTH.jl \
# input.sync \
# ${MAF} \
# ${DEPTH}
# ### Outputs:
# # (1) input_${MAF}MAF_filtered.sync
# # (2) input_${MAF}MAF_filtered_ALLELEFREQ.csv
# mv input_${MAF}MAF_filtered_ALLELEFREQ.csv geno.csv
#
# ##########################
# ### PREPARE PHENOTYPES ###
# ##########################
# echo -e "POPULATION,CLETHODIM_RESISTANCE,GLYPHOSATE_RESISTANCE,SULFOMETURON_RESISTANCE,TERBUTHYLAZINE_RESISTANCE" > pheno.csv
# echo -e "ACC01,0,0,0,NA" >> pheno.csv
# echo -e "ACC02,5.26315789473684,8.786741713571,NA,NA" >> pheno.csv
# echo -e "ACC03,7.31707317073171,38.0952380952381,0,NA" >> pheno.csv
# echo -e "ACC04,0,50,36.3636363636364,48.1458113526864" >> pheno.csv
# echo -e "ACC05,0,0,NA,NA" >> pheno.csv
# echo -e "ACC06,0,0,8.10810810810811,65" >> pheno.csv
# echo -e "ACC07,48.3870967741936,40.9188034188034,10,2.38095238095238" >> pheno.csv
# echo -e "ACC08,97.2222222222222,2.38095238095238,42.1052631578947,NA" >> pheno.csv
# echo -e "ACC09,58.6206896551724,95.1219512195122,15.7894736842105,50" >> pheno.csv
# echo -e "ACC10,25.7142857142857,5.12820512820513,8.10810810810811,NA" >> pheno.csv
# echo -e "ACC11,100,14.2857142857143,46.1538461538462,24.4158113526864" >> pheno.csv
# echo -e "ACC12,13.1578947368421,23.0769230769231,81.0810810810811,33.7558113526864" >> pheno.csv
# echo -e "ACC13,54.5454545454546,7.69230769230769,22.8571428571429,45.8458113526864" >> pheno.csv
# echo -e "ACC14,0,10.8108108108108,NA,NA" >> pheno.csv
# echo -e "ACC15,0,0,13.8888888888889,21.4058113526864" >> pheno.csv
# echo -e "ACC16,0,0,NA,NA" >> pheno.csv
# echo -e "ACC17,0,13.6363636363636,0,NA" >> pheno.csv
# echo -e "ACC18,0,97.5,NA,NA" >> pheno.csv
# echo -e "ACC19,0,2.38095238095238,0,NA" >> pheno.csv
# echo -e "ACC20,0,0,0,NA" >> pheno.csv
# echo -e "ACC21,0,0,0,0" >> pheno.csv
# echo -e "ACC22,0,0,48.7804878048781,15.7894736842105" >> pheno.csv
# echo -e "ACC23,65.7894736842105,10,NA,NA" >> pheno.csv
# echo -e "ACC24,0,100,13.3333333333333,NA" >> pheno.csv
# echo -e "ACC26,0,16.6666666666667,6.66666666666667,2.5" >> pheno.csv
# echo -e "ACC27,65.3846153846154,21.7391304347826,23.8095238095238,NA" >> pheno.csv
# echo -e "ACC28,0,32.4324324324324,NA,NA" >> pheno.csv
# echo -e "ACC29,2.63157894736842,7.89473684210526,9.09090909090909,NA" >> pheno.csv
# echo -e "ACC30,2.4390243902439,0,7.31707317073171,NA" >> pheno.csv
# echo -e "ACC31,5.12820512820513,26.6341256366723,12.1951219512195,18.9189189189189" >> pheno.csv
# echo -e "ACC32,0,0,0,NA" >> pheno.csv
# echo -e "ACC33,50,5.4054054054054,12.5,33.3333333333333" >> pheno.csv
# echo -e "ACC34,41.3793103448276,0,NA,NA" >> pheno.csv
# echo -e "ACC35,0,16.6666666666667,25.6410256410256,12.5058113526864" >> pheno.csv
# echo -e "ACC36,72.4137931034483,18.6206896551724,NA,NA" >> pheno.csv
# echo -e "ACC37,17.2413793103448,16.2162162162162,2.56410256410256,44.6558113526864" >> pheno.csv
# echo -e "ACC38,0,0,0,NA" >> pheno.csv
# echo -e "ACC39,0,0,43.3333333333333,NA" >> pheno.csv
# echo -e "ACC40,0,9.52380952380952,0,NA" >> pheno.csv
# echo -e "ACC41,0,100,NA,NA" >> pheno.csv
# echo -e "ACC42,0,100,NA,NA" >> pheno.csv
# echo -e "ACC43,0,66.6666666666667,NA,NA" >> pheno.csv
# echo -e "ACC44,0,0,NA,15.4858113526864" >> pheno.csv
# echo -e "ACC45,0,4.16666666666667,25,NA" >> pheno.csv
# echo -e "ACC46,0,3.84615384615385,9.09090909090909,NA" >> pheno.csv
# echo -e "ACC47,0,90,10.5263157894737,NA" >> pheno.csv
# echo -e "ACC48,0,31.7073170731707,26.8292682926829,NA" >> pheno.csv
# echo -e "ACC49,10,16.6666666666667,NA,36.3158113526864" >> pheno.csv
# echo -e "ACC50,0,85.7142857142857,28.2051282051282,29.7297297297297" >> pheno.csv
# echo -e "ACC51,0,11.1111111111111,NA,NA" >> pheno.csv
# echo -e "ACC52,12.5,2.38095238095238,5,NA" >> pheno.csv
# echo -e "ACC53,13.1578947368421,5.12820512820513,43.5897435897436,NA" >> pheno.csv
# echo -e "ACC54,0,0,0,31.5789473684211" >> pheno.csv
# echo -e "ACC55,0,0,6.9069069069069,85.7142857142857" >> pheno.csv
# echo -e "ACC57,0,7.14285714285714,0,NA" >> pheno.csv
# echo -e "ACC58,0,0,NA,NA" >> pheno.csv
# echo -e "ACC59,15.3846153846154,19.0476190476191,41.1764705882353,NA" >> pheno.csv
# echo -e "ACC60,85,9.52380952380952,5.12820512820513,91.8918918918919" >> pheno.csv
# echo -e "ACC61,0,10.7142857142857,35,47.3684210526316" >> pheno.csv
# echo -e "ACC62,0,74.7967479674797,20.5882352941176,58.5365853658537" >> pheno.csv
# ### separate phenotypes per herbicide resistance trait
# cut -d, -f1 pheno.csv > col1.temp
# for i in $(seq 2 $(awk -F"," '{print NF}' pheno.csv | head -n1))
# do
# 	cut -d, -f${i} pheno.csv > col2.temp
# 	fname=$(head -n1 col2.temp | cut -d\" -f2)
# 	paste -d, col1.temp col2.temp > pheno_${fname}.csv
# 	echo $fname
# done
# rm *.temp
#
# ###############################
# ###   GENOMIC PREDICTION		###
# ### 5-FOLD CROSS-VALIDATION ###
# ###############################
# ########################
# ### BUILD COVARIATES ###
# cd ${DIR}
# julia ${SRC_DIR}/build_PC_covariates.jl geno.csv
# ##########################
# ### WITHOUT COVARIATES ###
# cd ${DIR}
# mkdir NO_COVARIATE
# cd NO_COVARIATE/
# time \
# parallel julia ${SRC_DIR}/genomic_prediction_5fold_CV.jl \
# ${SRC_DIR%\/misc*}/src \
# ${DIR}/geno.csv \
# {1} \
# ${MAF} \
# ${NITERATIONS} \
# not_an_existing_file.txt ::: $(ls ${DIR}/pheno_*.csv)
# ls CROSS_VALIDATION_5_FOLD_OUTPUT_*.csv > OUTPUT_LIST.temp
# head -n1 $(head -n1 OUTPUT_LIST.temp) > CROSS_VALIDATION_5_FOLD_OUTPUT.csv
# for i in $(cat OUTPUT_LIST.temp)
# do
# tail -n+2 $i >> CROSS_VALIDATION_5_FOLD_OUTPUT.csv
# done
# Rscript ${SRC_DIR}/genomic_prediction_5fold_CV_plotting.r
# ##########################
# ### WITH PC1 COVARIATE ###
# cd ${DIR}
# mkdir PC1_X_COVARIATE
# cd PC1_X_COVARIATE/
# time \
# parallel julia ${SRC_DIR}/genomic_prediction_5fold_CV.jl \
# ${SRC_DIR%\/misc*}/src \
# ${DIR}/geno.csv \
# {1} \
# ${MAF} \
# ${NITERATIONS} \
# ${DIR}/PC1_X_cova.csv ::: $(ls ${DIR}/pheno_*.csv)
# ls CROSS_VALIDATION_5_FOLD_OUTPUT_*.csv > OUTPUT_LIST.temp
# head -n1 $(head -n1 OUTPUT_LIST.temp) > CROSS_VALIDATION_5_FOLD_OUTPUT.csv
# for i in $(cat OUTPUT_LIST.temp)
# do
# tail -n+2 $i >> CROSS_VALIDATION_5_FOLD_OUTPUT.csv
# done
# Rscript ${SRC_DIR}/genomic_prediction_5fold_CV_plotting.r
# ############################
# ### WITH PC1-5 COVARIATE ###
# cd ${DIR}
# mkdir PC1to5_X_COVARIATE
# cd PC1to5_X_COVARIATE/
# time \
# parallel julia ${SRC_DIR}/genomic_prediction_5fold_CV.jl \
# ${SRC_DIR%\/misc*}/src \
# ${DIR}/geno.csv \
# {1} \
# ${MAF} \
# ${NITERATIONS} \
# ${DIR}/PC1to5_X_cova.csv ::: $(ls ${DIR}/pheno_*.csv)
# ls CROSS_VALIDATION_5_FOLD_OUTPUT_*.csv > OUTPUT_LIST.temp
# head -n1 $(head -n1 OUTPUT_LIST.temp) > CROSS_VALIDATION_5_FOLD_OUTPUT.csv
# for i in $(cat OUTPUT_LIST.temp)
# do
# tail -n+2 $i >> CROSS_VALIDATION_5_FOLD_OUTPUT.csv
# done
# Rscript ${SRC_DIR}/genomic_prediction_5fold_CV_plotting.r
# #############################
# ### WITH PC1-10 COVARIATE ###
# cd ${DIR}
# mkdir PC1to10_X_COVARIATE
# cd PC1to10_X_COVARIATE/
# time \
# parallel julia ${SRC_DIR}/genomic_prediction_5fold_CV.jl \
# ${SRC_DIR%\/misc*}/src \
# ${DIR}/geno.csv \
# {1} \
# ${MAF} \
# ${NITERATIONS} \
# ${DIR}/PC1to10_X_cova.csv ::: $(ls ${DIR}/pheno_*.csv)
# ls CROSS_VALIDATION_5_FOLD_OUTPUT_*.csv > OUTPUT_LIST.temp
# head -n1 $(head -n1 OUTPUT_LIST.temp) > CROSS_VALIDATION_5_FOLD_OUTPUT.csv
# for i in $(cat OUTPUT_LIST.temp)
# do
# tail -n+2 $i >> CROSS_VALIDATION_5_FOLD_OUTPUT.csv
# done
# Rscript ${SRC_DIR}/genomic_prediction_5fold_CV_plotting.r
#
# ###########################################
# ### BLASTING PUTATIVE QTL FROM CV TESTS ###
# ###########################################
# DIR=/data/Lolium/Quantitative_Genetics/GPWAS_60_SEAUS
# SRC_DIR=/data/Lolium/Softwares/genomic_prediction/misc
# CORR=0.95 ### filter betas by model performance: i.e. >= to rediction accuracy defined by this correlation coefficient
# GENOME=/data/Lolium/Genomics/SEQUENCES/DNA/REFERENCE_GENOMES/lope_V1.0.fasta
# DB=/data/BlastDB/NCBI_NT0060
# cd ${DIR}
#
# #####################################################################################
# ### Extract putative QTL from cross-validation genomic prediction allelic effects ###
# #####################################################################################
# time \
# parallel -j 10 ${SRC_DIR}/extract_putative_QTL_parallel.sh ${DIR} ${SRC_DIR} ${CORR} {1} {2} {3} ::: \
# NO_COVARIATE PC1_X_COVARIATE PC1to5_X_COVARIATE PC1to10_X_COVARIATE ::: \
# CLETHODIM_RESISTANCE GLYPHOSATE_RESISTANCE SULFOMETURON_RESISTANCE TERBUTHYLAZINE_RESISTANCE ::: \
# FIXED_LS FIXED_RR FIXED_GLMNET FIXED_LASSO
#
# ###########################################################
# ### Reformat the genome assembly (no sequence wrapping) ###
# ###########################################################
#   grep "^>" ${GENOME} > SCAFFOLDS.temp
#   csplit ${GENOME} /\>/ '{*}'
#   ls xx* > SPLIT_FNAMES.temp
#   rm $(head -n1 SPLIT_FNAMES.temp) ### remove the first file which contains nada!
#   ls xx* > SPLIT_FNAMES.temp
#
#   rm REFORMATED_GENOME.fasta
#   touch REFORMATED_GENOME.fasta
#   # for i in $(seq 1 $(cat SCAFFOLDS.temp | wc -l))
#   time for i in $(seq 1 $(cat SPLIT_FNAMES.temp | wc -l))
#   do
#     echo $i
#     seq=$(head -n${i} SPLIT_FNAMES.temp | tail -n1)
#     head -n1 ${seq} > temp_name
#     tail -n+2 ${seq} > temp_seq
#     sed -zi 's/\n//g' temp_seq
#     cat temp_name >> REFORMATED_GENOME.fasta
#     cat temp_seq >> REFORMATED_GENOME.fasta
#     echo -e "" >> REFORMATED_GENOME.fasta ### insert newline charcter that was removed by sed
#   done
#   rm *temp* xx*
#
# ################
# ### Blasting ###
# ################
# cd ${DIR}
# query_halflen=3000 #3kb upstream and downstream of the SNP (putative QTL)
# parallel ${SRC_DIR}/blast_putative_QTL_parallel.sh \
# {1} \
# {2} \
# {3} \
# ${query_halflen} \
# ${DB} \
# REFORMATED_GENOME.fasta \
# ::: NO_COVARIATE PC1_X_COVARIATE PC1to5_X_COVARIATE PC1to10_X_COVARIATE \
# ::: CLETHODIM_RESISTANCE GLYPHOSATE_RESISTANCE SULFOMETURON_RESISTANCE TERBUTHYLAZINE_RESISTANCE \
# ::: FIXED_LS FIXED_RR FIXED_GLMNET FIXED_LASSO
#
# ### consolidate
# echo -e "COVARIATE\tTRAIT\tMODEL\tPUTATIVE_QTL\tSUBJECT_ID\tPERCENT_IDENTITY\tE_VALUE\tQUERY_COVERAGE\tBITSCORE\tHIT_INFO" > BLAST_OUT_CONSOLIDATED.txt
# for i in $(ls BLASTOUT-*)
# do
#   # i=$(ls BLASTOUT-* | head -n1)
#   echo $i
#   covariate=$(cut -d"-" -f2 <<<$i)
#   trait=$(cut -d"-" -f3 <<<$i)
#   model=$(cut -d"-" -f4 <<<$i)
#   nLines=$(cat $i | wc -l)
#   touch temp
#   for j in $(seq 1 $nLines); do echo -e "${covariate}\t${trait}\t${model}" >> temp; done
#   paste temp $i >> BLAST_OUT_CONSOLIDATED.txt
#   rm temp
# done
#
# ################
# ### Analysis ###
# ################
# ### Generate a word cloud of blastn top hits
# Rscript ${SRC_DIR}/generate_word_cloud_of_blastn_putative_QTL_hits.r BLAST_OUT_CONSOLIDATED.txt 10
#
# ################
# ### Clean-up ###
# ################
# mkdir PUTATIVE_QTL_OUTPUT; mv PUTATIVE_QTL-* PUTATIVE_QTL_OUTPUT/
# mkdir BLASTN_OUTPUT;  mv BLASTOUT-* BLASTN_OUTPUT/
# mkdir BLASTN_FASTA;  mv *.query.fasta BLASTN_FASTA/
# mkdir PLOTS;  mv *.svg PLOTS/
# rm *temp*
#
# # ###############################################################
# # ### IDENTIFYING THE VARIANTS FOR B-1,6 GLUCANASE-LIKE LOCUS ###
# # ###############################################################
# # ### setwd
# # DIR=/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/GP
# # cd ${DIR}
# # mkdir GLUCANASE_VARIANTS
# # cd ${DIR}/GLUCANASE_VARIANTS
# # ### Extract glucanase loci info
# # grep "glucanase" ../BLAST_OUT_CONSOLIDATED.txt > GLUCANASE_BLASTN_HITS.txt
# # cut -f4 GLUCANASE_BLASTN_HITS.txt | cut -d: -f1-2 | uniq > GLUCANASE_LOCI.txt
# # sed -i 's/:/\t/g' GLUCANASE_LOCI.txt
# # cut -f4 ../BLAST_OUT_CONSOLIDATED.txt | cut -d: -f1-2 | tail -n+2 | uniq > HITS_LOCI.txt
# # sed -i 's/:/\t/g' HITS_LOCI.txt
# # # ### copy alignments
# # # DIR_SAM=/data/Lolium/Population_Genetics/2018November_Ryegrass_Collection/SAM
# # # cp ${DIR_SAM}/*.bam .
# # ### copy sync file
# # DIR_SYNC=/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia
# # cp ${DIR_SYNC}/Lolium_2019_60pop.sync .
# # ### TESTS
# # grep -f GLUCANASE_LOCI.txt Lolium_2019_60pop.sync > GLUCANASE_SYNC.sync
# # grep -f HITS_LOCI.txt Lolium_2019_60pop.sync > HITS_SYNC.sync
# # sed 's/\t/,/g' GLUCANASE_LOCI.txt > GLUCANASE_LOCI.csv
# # grep -f GLUCANASE_LOCI.csv ../geno_ALLELEFREQ.csv > GLUCANASE_FREQ.csv
#
# ##################
# ##################
# ##################
# ################## NOPE NOPE NOPE!!!! SOMETHING WENT WRONG NUMBER OF READS IN THE SYNC FILE IS TOO SMALLLLL!!!!
# ################## RESTART RESTART!!!! 2019 10 18
# ##################
# ##################
# ##################
#
# #########################################
# ### NULL HYPOTHESIS: TEST RANDOM SNPS ###
# #########################################
# DIR=/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/GP/NULL
# mkdir $DIR
# cd $DIR
# cp ../AWK_10X/REFORMATED_GENOME.fasta .
# cp ../AWK_10X/LOCI_info.csv .
#
# NSAMPLES=1000
# query_halflen=3000 #3kb up- and down-stream
# NLOCI=$(cat LOCI_info.csv | wc -l)
# shuf -i1-${NLOCI} -n${NSAMPLES} | sort -g > IDX.temp
# touch null_loci.csv
# for i in $(cat IDX.temp)
# do
# 	echo $i
# 	head -n${i} LOCI_info.csv | tail -n1 >> null_loci.csv
# done
#
# rm null_query.fasta
# touch null_query.fasta
# for locus in $(cat null_loci.csv)
# do
#   # locus=$(cat null_loci.csv | head -n1)
# 	echo $locus
#   chrom=$(cut -d, -f1 <<<$locus)
#   # echo ${covariate}-${herbi}-${model}-${chrom}
#   pos=$(cut -d, -f2 <<<$locus)
#   grep -A 1 ${chrom} REFORMATED_GENOME.fasta | tail -n1 > null_seq.temp
#   bp=$(cat null_seq.temp | wc -c)
#   start=$(echo "${pos} - ${query_halflen}" | bc)
#   end=$(echo "${pos} + ${query_halflen}" | bc)
#   if [ $start -lt 1 ]
#   then
#     start=1
#   fi
#   if [ $end -gt $bp ]
#   then
#     end=${bp}
#   fi
#   echo -e "> ${chrom}:${pos}:${start}-${end}" >> null_query.fasta
#   cut -c${start}-${end} null_seq.temp >> null_query.fasta ||
#   cut -c1-${end} null_seq.temp >> null_query.fasta ||
#   cut -c${start}-${bp} null_seq.temp >> null_query.fasta ||
#   cut -c1-${bp} null_seq.temp >> null_query.fasta
#   rm null_seq.temp
# done
#
# DB=/data/BlastDB/NCBI_NT0060
# blastn -db ${DB} \
#   -query null_query.fasta \
#   -perc_identity 90 \
#   -outfmt '6 qseqid staxids pident evalue qcovhsp bitscore stitle' \
#   -out BLASTOUT-null.temp
# echo -e "COVARIATE\tTRAIT\tMODEL\tPUTATIVE_QTL\tSUBJECT_ID\tPERCENT_IDENTITY\tE_VALUE\tQUERY_COVERAGE\tBITSCORE\tHIT_INFO" > BLASTOUT-null.txt
# COVARIATE="NONE"
# TRAIT="NONE"
# MODEL="NONE"
# rm COL123.temp; touch COL123.temp
# for i in $(seq 1 $(cat BLASTOUT-null.temp | wc -l)); do echo -e "${COVARIATE}\t${TRAIT}\t${MODEL}" >> COL123.temp; done
# paste COL123.temp BLASTOUT-null.temp > COL123_BLASTOUT-null.temp
# cat COL123_BLASTOUT-null.temp >> BLASTOUT-null.txt
# rm *.temp
#
# SRC_DIR=/data/Lolium/Softwares/genomic_prediction/misc
# Rscript ${SRC_DIR}/generate_word_cloud_of_blastn_putative_QTL_hits.r BLASTOUT-null.txt 10
#
# ### RUBBISHSHSHSHSSHSHHSHS!!!!!!!!!
#
# ##############################################################
# ### MODELLING USING THE FULL DATA WITHOUT CROSS-VALIDATION ###
# ##############################################################
# DIR=/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/GP/FULL_DATA_NO_CV
# mkdir $DIR
# cd $DIR
# cp ../AWK_10X/geno_ALLELEFREQ.csv .
# cp ../AWK_10X/LOCI_info.csv .
# cp ../AWK_10X/pheno_* .
# cp ../AWK_10X/PC1*_X_cova.csv .
# cp ../AWK_10X/REFORMATED_GENOME.fasta .
#
# SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src
# MISC_DIR=/data/Lolium/Softwares/genomic_prediction/misc
# for i in $(ls pheno_*)
# do
# echo $i
# julia ${MISC_DIR}/genomic_prediction_full_data.jl \
# ${SRC_DIR} \
# geno_ALLELEFREQ.csv \
# ${i} \
# 0.001 \
# NULL \
# FIXED_LS,FIXED_RR,FIXED_GLMNET,FIXED_LASSO,FIXED_ITERATIVE
# # FIXED_ITERATIVE
# done
#
# ####################
# ### GFF3 THING-O ###
# ####################
# ### extracting loci info part og SNP_distribution_stats_plots.sh script too:
# DIR=/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/GP/FULL_DATA_NO_CV
# cd $DIR
# SNP_SYNC=/data/Lolium/Population_Genetics/2018November_Ryegrass_Collection/VCF/Lolium2019_10X_30PHRED_filtered.sync
# cut -f1-3 ${SNP_SYNC} > SYNC_LOCI.info
# cut -f1 SYNC_LOCI.info | uniq > SYNC_SCAFF.info
# ### download gff3
# wget http://185.45.23.197:5080/ryegrassdata/GENE_ANNOTATIONS/Lp_Annotation_V1.1.mp.gff3
# ### extract gff3 slices corresponding to the loci captured
# grep -f SYNC_SCAFF.info Lp_Annotation_V1.1.mp.gff3 > ANNOTATION_CAPTURED.gff3
