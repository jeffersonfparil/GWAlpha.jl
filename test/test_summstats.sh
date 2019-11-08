#!/bin/bash

################################
###                          ###
### SUMMARY STATISTICS TESTS ###
###                          ###
################################

DIR=/data/Lolium/Quantitative_Genetics/LOLSIM_2019_TEST
OUTDIR=${DIR}/TEST_SUMMSTATS                              # existing output directory to dump all the output files and folder
QUANTINEMO_DIR=/data/Lolium/Softwares/quantinemo_linux          # location of the quantinemo executable
GEN_PRED_SRC_DIR=${DIR}/Softwares/genomic_prediction/src  # full path to the genomic_prediction/src directory
nIndividuals=1000                                         # number of individuals to simulate
nLoci=1000                                                # number of loci to simulate (neutral + QTL)
nQTL=10                                                   # number of QTL among the all the loci simulated
nBGS=100                                                  # number of background selection loci (should be highly polygenic >= 100)
nAlleles=5                                                # number of alleles per loci, e.g. 5 for A,T,C,G, and DEL
allele_eff_model=CHISQ                                    # probability distribution to sample the allele effects from (i.e. CHISQ, NORM, & UNIF)
nGen=500                                                  # number of generations to simulate
nPop=16                                                    # number of populations/subpopulations to simulate (NOTE: must have a natural number square-root)
migration=0.001                                           # migration rate across the populations (surrently using the 1D stepping stone model see line 117)
selection=0.25                                             # selection intensity (trait of interest) defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
bg_selection=0.0                                          # background selection intensity defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
OUTPREFIX=LOLIUM_${nQTL}QTL_${migration}mr_${selection}fgs_${bg_selection}bgs   # prefix for the quantiNemo2 initiation (*.ini) file and the output folder
nPools=5                                                  # number of pools to subdivided each population
GEN_PRED_SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src  # full path to the location of the julia 1.0.5 binary
INI_FNAME=${OUTPREFIX}.ini                                # full path and filename of the quantinemo initiation (*.ini) file
nCores=5                                                  # number of cores to use for parallelisation
OUT_SUBDIR=${OUTDIR}/${OUTPREFIX}*/                       # output directory containing the parsed genotype and phenotype data
POP_ID_txt=${OUT_SUBDIR}/POP_ID.txt                       # full path and file name of the text file containing the population ID (strings: e.g. 1 or 01 or 16) of each individual across populations - one ID per line
GENOME_SPEC_FNAME=${OUTDIR}/Lperenne_genome.spec          # full path and file name of the genome specification file
QTL_SPEC_FNAME=${OUTDIR}/QTL_SPEC.csv                     # full path and file name of the QTL specification file
NPSTAT_DIR=${DIR}/Softwares/npstat                        # full path to the npstat directory
N=500                                                     # maximum number of simulated libraries sequenced (e.g. N=500)
GRADIENT=0

##############################
### GPWASim_01_simulate.sh ###
##############################
echo "#########################################################################################"
echo "SIMULATING GENOTYPE, QTL, AND PHENOTYPE DATA USING QUANTINENO2"

### prepare the quantinemo2 initialization file
### setwd
cd $OUTDIR

### Setup mating system
echo -e "mating_system 0" > ${OUTPREFIX}.ini                    #hermaphrodite random mating

### Demography
echo -e "patch_number $nPop" >> ${OUTPREFIX}.ini                # Number of populations or subpopulations to simulate: must have a natural number square-root
echo -e "generations $nGen" >> ${OUTPREFIX}.ini                 #Lolium has been introduced to Australia around 1880
echo -e "patch_capacity $nIndividuals" >> ${OUTPREFIX}.ini      #carrying capacity set to a billion
echo -e "regulation_model_offspring 1" >> ${OUTPREFIX}.ini      #regulation of population size to carrying capacity (patch_capacity) via random culling for offsprings
echo -e "regulation_model_adults 1" >> ${OUTPREFIX}.ini         #regulation of population size to carrying capacity (patch_capacity) via random culling for adults
echo -e "mating_nb_offspring_model 0" >> ${OUTPREFIX}.ini       #total number of offspring is set to carrying capacity
echo -e "mean_fecundity 1" >> ${OUTPREFIX}.ini                  #... 1 which means a constant population size which is just an approximation though I don't know how bad of an approximation it is
echo -e "dispersal_rate $migration" >> ${OUTPREFIX}.ini         #migration rate across adjacent patches [@Wang2016]
echo -e "dispersal_model 3" >> ${OUTPREFIX}.ini                 #2D stepping-stone
echo -e "dispersal_lattice_range 1" >> ${OUTPREFIX}.ini         #2D stepping-stone: dispersal range of 8 adjacent patches (vertical, horizontal and dagonal; m/8)
echo -e "dispersal_border_model 2" >> ${OUTPREFIX}.ini          #absorbing boundaries: migration beyond the border is lost
### dispersal_lattice_dims by defualt sets a square patches structure sqrt(nPop) x sqrt(nPop)
###### Setting the genetic map
Rscript ${GEN_PRED_SRC_DIR}/GPWASim_01_build_genome.r \
  7 \
  97.7,151.5,63.3,119.2,89.1,115.2,113.7  \
  $nLoci
  ### Input:
    ### 7 chromosomes in Lolium
    ### Chromosome length in cM from Pfeifer, 2015
    ### number of loci
  ### Output: genome_loci_*.temp
echo -e "quanti_genome { {1: $(cat genome_loci_1.temp)}" >> ${OUTPREFIX}.ini
echo -e "\t\t{2: $(cat genome_loci_2.temp)}" >> ${OUTPREFIX}.ini
echo -e "\t\t{3: $(cat genome_loci_3.temp)}" >> ${OUTPREFIX}.ini
echo -e "\t\t{4: $(cat genome_loci_4.temp)}" >> ${OUTPREFIX}.ini
echo -e "\t\t{5: $(cat genome_loci_5.temp)}" >> ${OUTPREFIX}.ini
echo -e "\t\t{6: $(cat genome_loci_6.temp)}" >> ${OUTPREFIX}.ini
echo -e "\t\t{7: $(cat genome_loci_7.temp)} }" >> ${OUTPREFIX}.ini
### Genotype configuration
echo -e "quanti_loci_1 $(echo $nLoci - $nBGS | bc)" >> ${OUTPREFIX}.ini                #[trait of interest: QTL-based] total number of loci: neutral and QTL
echo -e "quanti_loci_2 $nBGS" >> ${OUTPREFIX}.ini                #[background selection trait]total number of loci: neutral and QTL
echo -e "quanti_all $nAlleles" >> ${OUTPREFIX}.ini              #number of alleles per locus; and since we're simulating SNPs we want 5: A,T,C,G and DEL excluding N
echo -e "quanti_nb_trait 2" >> ${OUTPREFIX}.ini                 #number of traits: 2 - one for the quantitative trait of interest and the other as background selection simulator
    ###### Build the quanti_allelic_file for the QTL and background selection (BGS) loci specs
    echo -e "# Quantitative Alleles Specifications File" > QTL.spec ### loci specifications for the trait of interest
    echo -e "[FILE_INFO]{" >> QTL.spec
    echo -e "  col_locus 1" >> QTL.spec
    echo -e "  col_allele 2" >> QTL.spec
    echo -e "  col_allelic_value 3" >> QTL.spec
    echo -e "  col_mut_freq 4" >> QTL.spec
    echo -e "  col_ini_freq 5" >> QTL.spec
    echo -e "}" >> QTL.spec
    echo -e "#locus\tallele\tvalue\tmut_freq\tini_freq" >> QTL.spec
    echo -e "# Quantitative Alleles Specifications File" > BGS.spec ### loci specifications for the background selection trait
    echo -e "[FILE_INFO]{" >> BGS.spec
    echo -e "  col_locus 1" >> BGS.spec
    echo -e "  col_allele 2" >> BGS.spec
    echo -e "  col_allelic_value 3" >> BGS.spec
    echo -e "  col_mut_freq 4" >> BGS.spec
    echo -e "  col_ini_freq 5" >> BGS.spec
    echo -e "}" >> BGS.spec
    echo -e "#locus\tallele\tvalue\tmut_freq\tini_freq" >> BGS.spec
    ###### Empirical allele count distribution in Lolium (Lolium2018 Illumina sequencing): script found in "misc/quantinemo2_geno_sim_00_QTL_build_testing.r" commented-out
    echo -e "ALLELE_COUNT_1,0.0" > Allele_counts_dist.spec
    echo -e "ALLELE_COUNT_2,0.673749892530331" >> Allele_counts_dist.spec
    echo -e "ALLELE_COUNT_3,0.257344297077753" >> Allele_counts_dist.spec
    echo -e "ALLELE_COUNT_4,0.0660887235376553" >> Allele_counts_dist.spec
    echo -e "ALLELE_COUNT_5,0.00281708685426013" >> Allele_counts_dist.spec
    echo -e "ALLELE_COUNT_6,0.0" >> Allele_counts_dist.spec ### input for the GPWASim_01_build_loci_spec.r script
    ###### Simulate the loci specifications for trait of interest (QTL + zero-effect) and background selection trait
    Rscript ${GEN_PRED_SRC_DIR}/GPWASim_01_build_loci_spec.r \
      $nLoci \
      $nQTL \
      $nBGS \
      Allele_counts_dist.spec \
      $allele_eff_model \
      $nIndividuals
      ### Input:
        ### total number of loci in the genome
        ### number of QTL for the trait of interest
        ### number of QTL for the background selection trait
        ### frequency table of the number of alleles per locus
        ### probability distribution model to sample the allele effects from (the same for both trait of interest and the background selection trait)
      ### Output:
        ### (1) "BGS.spec"
        ### (2) "BGS_idx_map_to_genome.spec"
        ### (3) "QTL.spec"
        ### (4) "QTL_idx_map_to_genome.spec"
        ### (5) "QTL_min_max_GEBV.spec"
echo -e "quanti_allelic_file_1 QTL.spec" >> ${OUTPREFIX}.ini      #allelic file for the QTL (trait of interest) + zero-effect loci
echo -e "quanti_allelic_file_2 BGS.spec" >> ${OUTPREFIX}.ini      #alellic file for the background selection loci
echo -e "quanti_locus_index_1 {$(uniq QTL_idx_map_to_genome.spec)}" >> ${OUTPREFIX}.ini
echo -e "quanti_locus_index_2 {$(uniq BGS_idx_map_to_genome.spec)}" >> ${OUTPREFIX}.ini
echo -e "quanti_ini_allele_model 0" >> ${OUTPREFIX}.ini         #genotypes are set to be maximally polymorph while set to 1 for monomorph or fixed alleles per genotype
echo -e "quanti_mutation_rate 0" >> ${OUTPREFIX}.ini            #mutation rate (set to no mutation for simplicity now)
echo -e "quanti_mutation_model 0" >> ${OUTPREFIX}.ini           #random mutation model (see pages 56-60)
### Phenotype settings
echo -e "quanti_environmental_model 0" >> ${OUTPREFIX}.ini      #environmental variance constant set as quanti_heritability=50%
echo -e "quanti_heritability 0.50" >> ${OUTPREFIX}.ini          #narrow-sense heritability set to 50%
echo -e "quanti_selection_model 2" >> ${OUTPREFIX}.ini          #directional selection (Richards, 1959)
echo -e "quanti_dir_sel_min 0.0" >> ${OUTPREFIX}.ini          #minimum fitness
echo -e "quanti_dir_sel_max 1.0" >> ${OUTPREFIX}.ini          #maximum fitness
MIN_GEBV=$(tail -n1 QTL_min_max_GEBV.spec | cut -f1)                      #minimum genomic breeding value
MAX_GEBV=$(tail -n1 QTL_min_max_GEBV.spec | cut -f2)                      #maximum genomic breeding value
MID_GEBV=$(echo "(($MAX_GEBV - $MIN_GEBV)) / 2" | bc)                     #half-way between minimum and maximum GEBV
echo -e "quanti_dir_sel_max_growth $MID_GEBV" >> ${OUTPREFIX}.ini       #maximal slope half-way between minimum and maximu possible GEBV
# echo -e "quanti_dir_sel_max_growth_1 $MID_GEBV" >> ${OUTPREFIX}.ini       #maximal slope half-way between minimum and maximu possible GEBV
# MIN_GEBV=$(tail -n1 BGS_min_max_GEBV.spec | cut -f1)                      #same thing as above but for the background selection (BGS) trait and QTL
# MAX_GEBV=$(tail -n1 BGS_min_max_GEBV.spec | cut -f2)
# MID_GEBV=$(echo "(($MAX_GEBV - $MIN_GEBV)) / 2" | bc)
# echo -e "quanti_dir_sel_max_growth_2 $MID_GEBV" >> ${OUTPREFIX}.ini
# echo -e "quanti_dir_sel_symmetry 1" >> ${OUTPREFIX}.ini                   #symmetric logistic selection curve
echo -e "quanti_dir_sel_growth_rate_1 $selection" >> ${OUTPREFIX}.ini     #selection intensity as the slope of the directional selection logistic curve [-inf, +inf]--> but we'll try 0.5, 1, 2, 10??? [QTL]
echo -e "quanti_dir_sel_growth_rate_2 $bg_selection" >> ${OUTPREFIX}.ini  #selection intensity as the slope of the directional selection logistic curve [-inf, +inf]--> but we'll try 0.5, 1, 2, 10??? [background selection]
### Summary Statistics Output
### Save genotype and phenotype data
gen_interval=$(echo "$nGen / 10" | bc)
if [ $gen_interval -lt 1 ]; then gen_interval=$(echo $nGen); fi
echo -e "quanti_genot_logtime $gen_interval" >> ${OUTPREFIX}.ini
echo -e "quanti_save_genotype 1" >> ${OUTPREFIX}.ini            #output in FSTAT format: main genotype block after the loci ID: n x l --> number of individuals x number of loci: [1:5][1:5] --> genotype format, e.g. 25 for a heterozygote containing the 2nd allele and the 5th allele
echo -e "quanti_genot_filename $OUTPREFIX" >> ${OUTPREFIX}.ini
echo -e "quanti_phenot_logtime $gen_interval" >> ${OUTPREFIX}.ini
echo -e "quanti_save_phenotype 1" >> ${OUTPREFIX}.ini
echo -e "quanti_phenot_filename $OUTPREFIX" >> ${OUTPREFIX}.ini
echo -e "Chrom\tMbp\tcM" > Lperenne_genome.spec
echo -e "1\t470.30\t97.7" >> Lperenne_genome.spec
echo -e "2\t429.90\t151.5" >> Lperenne_genome.spec
echo -e "3\t406.56\t63.3" >> Lperenne_genome.spec
echo -e "4\t369.03\t119.2" >> Lperenne_genome.spec
echo -e "5\t339.41\t89.1" >> Lperenne_genome.spec
echo -e "6\t322.62\t115.2" >> Lperenne_genome.spec
echo -e "7\t284.07\t113.7" >> Lperenne_genome.spec
### generate loci specifications fle: CHROM, r (cM), POS (bp) and the QTL spec file for the trait of interest (excluding the allele effect of the BGS QTL)
ls genome_loci_*.temp | sed s/.temp//g | rev | cut -d'_' -f1 > chrom_id.temp
cat genome_loci_*.temp > r_pos.temp
paste -d' ' chrom_id.temp r_pos.temp > chrom_r.temp ### locus positions in cM
Rscript ${GEN_PRED_SRC_DIR}/GPWASim_01_genome_QTL_spec_parsing.r \
  chrom_r.temp \
  Lperenne_genome.spec \
  QTL.spec \
  QTL_idx_map_to_genome.spec
### Input:
### (1) locus positions per chromosome (NO HEADER;space-delimited; 1 chrom/line)
### (2) chromosome specifications (HEADER: Chrom, Mbp, cM; tab-delimited)
### (3) QTL specifications file (HEADER: metadata-see above; tab-delimited; excludes BGS QTL)
### (4) QTL indices that maps to the genome (NO HEADER; tab-delimited; 1 column)
### Output:
### (1) "GENOME_SPEC.csv" (HEADER: CHROM, r, POS)
### (2) "QTL_SPEC.csv" (HEADER: CHROM, POS, ALLELE, EFFECT)
### clean-up
rm *.temp
### EXECUTION
${QUANTINEMO_DIR}/quantinemo ${OUTPREFIX}.ini

###########################
### GPWASim_02_parse.sh ###
###########################
cd ${OUTDIR}/${OUTPREFIX}*/
ls *.dat > DAT_LIST
for i in $(cat DAT_LIST)
do
  # fname_dat=${OUTPREFIX}_g${i}.dat
  echo ${i}
  ${GEN_PRED_SRC_DIR}/GPASim_02_parse.sh \
        ${OUTDIR} \
        ${i} \
        ${nPools} \
        ${GEN_PRED_SRC_DIR} \
        ${INI_FNAME} \
        ${nCores}
  mv ORIGINAL_OUTPUT/*.* .
done
for i in $(cat DAT_LIST); do mv $i ORIGINAL_OUTPUT/; done

##############
### NPSTAT ###
##############
cd ${OUTDIR}/${OUTPREFIX}*/
nChrom=$(echo $(cat ${GENOME_SPEC_FNAME} | wc -l) - 1 | bc)
# window_size=10000000 #10 Mb
# window_size=1000000 #1 Mb
window_size=100000 #100 kb
for i in $(ls *_POPULATION.pileup)
do
  pop_size=10 #fixed to 10 for quick estimations though highly over-estimated summary statistics
  parallel \
  ${GEN_PRED_SRC_DIR}/GPASim_03_NPSTAT_parallel.sh \
          {1} \
          ${i} \
          ${pop_size} \
          ${window_size} \
          ${NPSTAT_DIR} ::: $(seq 1 $nChrom)
  echo -e "chrom\twindow\tcoverage\toutgroup_coverage\tdepth\tS\tWatterson_estimator\tTajimas_Pi\tTajimas_D\tFay_Wu_H_unorm\tFay_Wu_H_norm\tS_var\tV_theta\toutgroup_divergence\tnonsynonimous_polymorph\tsynonymous_polymorph\tnonsynonimous_divergence\tsynonymous_divergence\talpha" > ${i%.pileup*}_NPSTAT.stats
  cat ${i%.pileup*}_CHR*_NPSTAT.stats >> ${i%.pileup*}_NPSTAT.stats
  rm ${i%.pileup*}_CHR*_NPSTAT.stats
done
### MERGE NPSTAT OUTPUT ACROSS GENERATIONS
echo -e "gen\tpop\tchrom\twindow\tcoverage\toutgroup_coverage\tdepth\tS\tWatterson_estimator\tTajimas_Pi\tTajimas_D\tFay_Wu_H_unorm\tFay_Wu_H_norm\tS_var\tV_theta\toutgroup_divergence\tnonsynonimous_polymorph\tsynonymous_polymorph\tnonsynonimous_divergence\tsynonymous_divergence\talpha" > ALL_POP_NPSTAT.stats
for i in $(ls *_POPULATION_NPSTAT.stats)
do
  l=$(cat $i | wc -l)
  gen=$(cut -d_ -f6 <<<$i)
  pop=$(cut -d_ -f7 <<<$i)
  touch COL1.temp
  for j in $(seq 1 $l); do echo -e "${gen}\t${pop}" >> COL1.temp; done
  paste COL1.temp ${i} > npstat_with_genID.temp
  tail -n+2 npstat_with_genID.temp >> ALL_POP_NPSTAT.stats
  rm *.temp
done

###############
### CLEANUP ###
###############
cd ${OUTDIR}/${OUTPREFIX}*/
mkdir MISC/
mv *.* MISC/
mv DAT_LIST MISC/
mv MISC/ALL_POP_NPSTAT.stats .

################
### ANALYSIS ###
################
cd ${OUTDIR}/${OUTPREFIX}*/
Rscript ${GEN_PRED_SRC_DIR%src*}test/test_summstats.r ALL_POP_NPSTAT.stats g${nGen} ${window_size}
