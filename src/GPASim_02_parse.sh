#!/bin/bash

####################################
### PARSE THE QUANTINEMO2 OUTPUT ###
####################################

### INPUT
OUTDIR=$1             #quantinemo output directory containing the initiation file, spec files and the output subdirectory with the phenotype and genotype files
fname_dat=$2          #filename of the genotype (*.dat) file within the subdirectory of $OUTDIR
nPools=$3             #number of pools to subdivided each population
GEN_PRED_SRC_DIR=$4   #full path to the genomic_prediction/src directory
INI_FNAME=$5          #full path and filename of the quantinemo initiation (*.ini) file
nCores=$6             #number of cores to use for parallelisation

### SET WORKING DIRECTORY AND EXTRACT BASIC INFO
MOST_RECENT_OUTPUT_DIR=$(ls -dt ${OUTDIR}/*/ | grep ${INI_FNAME%.ini*} | tail -n1)
cd ${MOST_RECENT_OUTPUT_DIR}
fname_phe=${fname_dat%.*}.phe
fname_array=(${fname_dat//_g/ })
nGen=${fname_array[-1]%.*}
nPop=$(head -n1 $fname_dat | cut -d' ' -f1)
nLoci=$(head -n1 $fname_dat | cut -d' ' -f2)
nAlleles=$(head -n1 $fname_dat | cut -d' ' -f3)

####################
### SAMPLE EXECUTION
# DIR=/data/Lolium/Quantitative_Genetics/LOLSIM_2019
# rep=1
# nGen=500
# nQTL=10
# migration=0.001
# selection=0.25
# bg_selection=-0.25
# GRADIENT=0
# OUTPREFIX=LOLIUM_${rep}rep_${nQTL}QTL_${migration}mr_${selection}fgs_${bg_selection}bgs_${GRADIENT}grad
# OUTDIR=${DIR}/Output_${rep}rep_${nQTL}QTL_${migration}mr_${selection}fgs_${bg_selection}bgs_${GRADIENT}grad/
# fname_dat=${OUTPREFIX}_g${nGen}.dat
# nPools=5
# GEN_PRED_SRC_DIR=${DIR}/Softwares/genomic_prediction/src
# INI_FNAME=${OUTPREFIX}.ini
# nCores=$(echo $(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l) / 2 | bc) #half the number of parallel threads for RAM-intensive parsing
#
# time \
# ${GEN_PRED_SRC_DIR}/GPASim_02_parse.sh \
# ${OUTDIR} \
# ${fname_dat} \
# ${nPools} \
# ${GEN_PRED_SRC_DIR} \
# ${INI_FNAME} \
# ${nCores}
####################

echo "#########################################################################################"
echo "PARSING QUANTINEMO2 OUTPUT, POOLING, AND BUILDING GENOMIC PREDCTION INPUT FILES"

### separate genotype data of each patch (population) in bash
NHEADERS=$(echo $nLoci + 1 | bc)
NTAIL=$(echo "$NHEADERS + 1" | bc)
head -n  ${NHEADERS} ${fname_dat} > DAT_HEADER.temp
tail -n +${NTAIL} ${fname_dat} > DAT_GENOTYPES.temp
cut -d' ' -f1 DAT_GENOTYPES.temp > POP_ID.txt ### strings of population ID of each individual across populations

for pop in $(uniq POP_ID.txt)
do
  grep "^${pop} " DAT_GENOTYPES.temp > POP_GREP.temp
  cat DAT_HEADER.temp POP_GREP.temp > ${fname_dat%.dat*}_p${pop}.dat
done

### separate phenotype data of each patch (population) in bash
### NOTE: bug in *.phe file - no grouping based on patches in the first column!
NHEADERS=2
NHEADERS_PLUS_ONE=$(echo $NHEADERS + 2 | bc) # line1: nPop nTraits; line2: phenotype of interest (QTL); line3: background selection phenotype
head -n  ${NHEADERS} ${fname_phe} > PHE_HEADER.temp
tail -n +${NHEADERS_PLUS_ONE} ${fname_phe} | cut -f2 > PHE_PHENOTYPES_NOID.temp ### isolating the QTL-derived phenotypic values (the background-selection phenotypes are in the 3rd column)
paste POP_ID.txt PHE_PHENOTYPES_NOID.temp > PHE_PHENOTYPES.temp

for pop in $(uniq POP_ID.txt)
do
  grep ^${pop} PHE_PHENOTYPES.temp > POP_GREP.temp
  cat PHE_HEADER.temp POP_GREP.temp > ${fname_phe%.phe*}_p${pop}.phe
done

### cleanup
mkdir ORIGINAL_OUTPUT/
mv ${fname_dat} ORIGINAL_OUTPUT/
mv ${fname_phe} ORIGINAL_OUTPUT/
mv PHE_PHENOTYPES.temp PHENOTYPES_WITH_POPID.txt ### save phonotype file with population grouping id column 1
rm *.temp

### julia parsing script for each population
cd $OUTDIR
uniq ${MOST_RECENT_OUTPUT_DIR}/POP_ID.txt > POP_UNIQ.temp ### list the population ID names to be used for parallel execution of the julia parsing script
parallel -j $nCores julia ${GEN_PRED_SRC_DIR}/GPASim_02_parse.jl \
                                            ${MOST_RECENT_OUTPUT_DIR}/${fname_dat%.dat*}_p{1}.dat \
                                            GENOME_SPEC.csv \
                                            QTL_min_max_GEBV.spec \
                                            ${nPools} \
                                            100 ::: $(cat POP_UNIQ.temp)
rm POP_UNIQ.temp
        ##############
        ### inputs ###
        ##############
        ### (1) dat_fname = quantinemo2 genotype FSTAT file outpu
        ### (2) genome_spec_fname = genome specificications: chromosome length in cM and bp
        ### (3) QTL_min_max_GEBV = minimum and maximum possible genotypic value specification based on the QTL effects
        ### (4) Number of equally sized pools to group individual genotypes into
        ### (5) Read depth to simulate for building pileup files from the FSTAT genotype output file
        ###############
        ### outputs ###
        ###############
        ### (1) Genotype matrix (*_GENO_.csv) (nind x nloci*5) (NO HEADER)
        ### (2) Loci information (*_LOCI_SPEC.csv) (nloci x 3) (HEADER: chromosome, position in cM, position in bp)
        ### (3) Phenotype data (*_PHENO.csv) (nind x 3) (HEADER: individual ID, untransformed phenotypic value, transformed phenotyp value to range from 0-1 based on possible extremes based on allele effects)
        ### (4) Pooling percentiles (*_POOLS_PERCENTILES.csv) (npools+1 x 1) (NO HEADER)
        ### (5) Pool sizes and the corresponding mean phenotypic values per pool (*_POOLS_PHENO.csv) (npools x 2) (NO HEADER)
        ### (6) Allele frequency matrix (*_POOLS_GENO.csv) (npools x nloci*5) (NO HEADER)
        ### (7) Synchronized pileup file (*_POOLS_GENO.sync) (nloci x npools+3) (NO HEADER)
        ### (8) Pooled phenotype file in tandem with the sync file (*_POOLS_PHENO.py) (6 x 1) (NO HEADER: pheno_name, sig, min, max, perc, & q)
        ### (9) QTL information (*_QTL_SPEC.csv) (nQTL*5 x 4) (HEADER: chromosome,position in bp, allele, effect)
        ### (10) Pileup file format combining individual genotype into a single population of simulated sequencing reads (*_POPULATION.pileup) (nloci x 6) (NO HEADER)
### generate synchronized pileup file where each pool is one whole population --> to simulate Pool-seq per population
ls ${MOST_RECENT_OUTPUT_DIR}/${fname_dat%.dat*}_p*_GENO.csv | grep -v POOL > fname_list.txt
julia ${GEN_PRED_SRC_DIR}/GPASim_02_pool.jl fname_list.txt $nLoci $nAlleles ${MOST_RECENT_OUTPUT_DIR}/PHENOTYPES_WITH_POPID.txt QTL_min_max_GEBV.spec
        ################
        #### inputs ####
        ################
        #### fname_list_txt = text file listing the filenames of the genotype data output of quantinemo2_geno_sim_01a_OUTPUTPARSE.jl
        #### nLoci = number of loci in the genotype files
        #### nAlleles = number of alleles per loci in the genotype files
        #### PHENOTYPES_WITH_POPID.txt = tab-delimited file of phenotype values of each individual per population
        #### QTL_min_max_GEBV.spec = the minimum ad maximum possible phenotype values
        #################
        #### outputs ####
        #################
        #### (1) synchronized pileup file (*_ALLPOP_GENO.sync) (nLoci x nPop+3) (NO HEADER)
        #### (2) phenotype means per population (*_ALLPOP_PHENO.pool) (HEADER: pop_id{string}, y_untr{untransformed phenotype}, y{tranformed phenotype ranges from 0 to 1})
echo "#########################################################################################"
