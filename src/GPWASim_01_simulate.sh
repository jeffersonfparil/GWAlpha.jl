#!/bin/bash

#####################################################################
### SIMULATE GENOTYPE, QTL, AND PHENOTYPE DATA USIONG QUANTINEMO2 ###
#####################################################################

### INPUT
QUANTINEMO_DIR=${1}         # location of the quantinemo executable
GEN_PRED_SRC_DIR=${2}       #full path to the genomic_prediction/src directory
OUTDIR=${3}                 # existing output directory to dump all the output files and folder
OUTPREFIX=${4}              # prefix for the quantiNemo2 initiation (*.ini) file and the output folder
nIndividuals=${5}           # number of individuals to simulate
nLoci=${6}                  # number of loci to simulate (neutral + QTL)
nQTL=${7}                   # number of QTL among the all the loci simulated
nBGS=${8}                   # number of background selection loci (should be highly polygenic >= 100)
nAlleles=${9}               # number of alleles per loci, e.g. 5 for A,T,C,G, and DEL
allele_eff_model=${10}      # probability distribution to sample the allele effects from (i.e. CHISQ, NORM, & UNIF)
nGen=${11}                  # number of generations to simulate
nPop=${12}                  # number of populations/subpopulations to simulate (NOTE: must have a natural number square-root)
migration=${13}             # migration rate across the populations (surrently using the 2D stepping stone model see line 117) [@Wang2017]
selection=${14}             # selection intensity (nQTL) defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
bg_selection=${15}          # background selection intensity (nBGS) defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
GRADIENT=${16}              # gradient of QTL (foregrund and background selection) across the landscape: 0 for uniform, 1 for first row of populations, and 2 for first and last rows of populations are the sources of the the non-zero (non-wild type) alleles

####################
### SAMPLE EXECUTION
# OUTDIR=/data/Lolium/Quantitative_Genetics/TEST_LANDSCAPING      # existing output directory to dump all the output files and folder
# DIR=/data/Lolium
# QUANTINEMO_DIR=${DIR}/Softwares/quantinemo_linux                # location of the quantinemo executable
# GEN_PRED_SRC_DIR=${DIR}/Softwares/genomic_prediction/src        # full path to the genomic_prediction/src directory
# nIndividuals=1000                                               # number of individuals to simulate
# nLoci=1000                                                      # number of loci to simulate (neutral + QTL)
# nQTL=10                                                         # number of QTL among the all the loci simulated
# nBGS=100                                                        # number of background selection loci (should be highly polygenic >= 100)
# nAlleles=5                                                      # number of alleles per loci, e.g. 5 for A,T,C,G, and DEL
# allele_eff_model=CHISQ                                          # probability distribution to sample the allele effects from (i.e. CHISQ, NORM, & UNIF)
# nGen=500                                                        # number of generations to simulate
# nPop=16                                                         # number of populations/subpopulations to simulate (NOTE: must have a natural number square-root)
# migration=0.00                                                  # migration rate across the populations (surrently using the 1D stepping stone model see line 117)
# selection=0.25                                                  # selection intensity (trait of interest) defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
# bg_selection=0.00                                               # background selection intensity defined as the slope of the directional selection logistic curve (Richards, 1959): ranges from -Inf (select for low phen) to +Inf (select for high phen)
# OUTPREFIX=LOLSIM_${nQTL}QTL_${migration}mr_${selection}fgs_${bg_selection}bgs   # prefix for the quantiNemo2 initiation (*.ini) file and the output folder
# GRADIENT=0                                                      # uniformly distributed non-wildtype alleles
#
# time \
# ${GEN_PRED_SRC_DIR}/GPWASim_01_simulate.sh \
#       $QUANTINEMO_DIR \
#       $GEN_PRED_SRC_DIR \
#       $OUTDIR \
#       $OUTPREFIX \
#       $nIndividuals \
#       $nLoci \
#       $nQTL \
#       $nBGS \
#       $nAlleles \
#       $allele_eff_model \
#       $nGen \
#       $nPop \
#       $migration \
#       $selection \
#       $bg_selection \
#       $GRADIENT
####################

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
# echo -e "mating_nb_offspring_model 3" >> ${OUTPREFIX}.ini       #simple fecundity by rouding the number of offsrpings) with fecundity rate of...
echo -e "mating_nb_offspring_model 0" >> ${OUTPREFIX}.ini       #total number of offspring is set to carrying capacity
echo -e "mean_fecundity 1" >> ${OUTPREFIX}.ini                  #... 1 which means a constant population size which is just an approximation though I don't know how bad of an approximation it is
echo -e "dispersal_rate $migration" >> ${OUTPREFIX}.ini         #migration rate across adjacent patches [@Wang2016]
# echo -e "dispersal_model 2" >> ${OUTPREFIX}.ini                 #1-dimensional stepping-stone model
# echo -e "dispersal_border_model 2" >> ${OUTPREFIX}.ini          #migrants from the border gets lost for the 1D stepping-stone model
echo -e "dispersal_model 3" >> ${OUTPREFIX}.ini                 #2D stepping-stone
echo -e "dispersal_lattice_range 1" >> ${OUTPREFIX}.ini         #2D stepping-stone: dispersal range of 8 adjacent patches (vertical, horizontal and dagonal; m/8)
echo -e "dispersal_border_model 2" >> ${OUTPREFIX}.ini          #absorbing boundaries: migration beyond the border is lost
### dispersal_lattice_dims by defualt sets a square patches structure sqrt(nPop) x sqrt(nPop)
# n_rows_cols=$(echo "sqrt(($nPop))" | bc)
# echo -e "dispersal_lattice_dims ($n_rows_cols, $n_rows_cols)" >> ${OUTPREFIX}.ini        #structure of population in a sqrt(n) x sqrt(n) square grid

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
    cp QTL.spec QTL_WT.spec
    cp BGS.spec BGS_WT.spec
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
        ### effective population size or in this case the number of individuals
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
echo -e "quanti_dir_sel_max_growth $MID_GEBV" >> ${OUTPREFIX}.ini         #maximal slope half-way between minimum and maximu possible GEBV
# echo -e "quanti_dir_sel_symmetry 1" >> ${OUTPREFIX}.ini                   #symmetric logistic selection curve
echo -e "quanti_dir_sel_growth_rate_1 $selection" >> ${OUTPREFIX}.ini     #selection intensity as the slope of the directional selection logistic curve [-inf, +inf]--> but we'll try 0.5, 1, 2, 10??? [QTL]
echo -e "quanti_dir_sel_growth_rate_2 $bg_selection" >> ${OUTPREFIX}.ini  #selection intensity as the slope of the directional selection logistic curve [-inf, +inf]--> but we'll try 0.5, 1, 2, 10??? [background selection]

### Summary Statistics Output
# echo -e "stat_log_time $nGen" >> ${OUTPREFIX}.ini               #output summary statistics for the last generation only
# echo -e "stat {q.adlt.fst" >> ${OUTPREFIX}.ini
# echo -e "q.adlt.fst.wc_pair" >> ${OUTPREFIX}.ini
# echo -e "q.adlt.R2}" >> ${OUTPREFIX}.ini
### Save genotype and phenotype data
echo -e "quanti_genot_logtime $nGen" >> ${OUTPREFIX}.ini        #write genotype output for the final generation only
echo -e "quanti_save_genotype 1" >> ${OUTPREFIX}.ini            #output in FSTAT format: main genotype block after the loci ID: n x l --> number of individuals x number of loci: [1:5][1:5] --> genotype format, e.g. 25 for a heterozygote containing the 2nd allele and the 5th allele
echo -e "quanti_genot_filename $OUTPREFIX" >> ${OUTPREFIX}.ini
echo -e "quanti_phenot_logtime $nGen" >> ${OUTPREFIX}.ini       #write phenotype output for the final generation only
# echo -e "quanti_phenot_logtime 10" >> ${OUTPREFIX}.ini        #test
echo -e "quanti_save_phenotype 1" >> ${OUTPREFIX}.ini
echo -e "quanti_phenot_filename $OUTPREFIX" >> ${OUTPREFIX}.ini
### write out chromosome lengths in bp (@Ansari2016) and cM (@Pfeifer2013)
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

### Generate intital genotypes
if [ $GRADIENT == 0 ]
then
  echo -e "Using NULL gradient of allele effects across the landscape."
elif [ $GRADIENT == 1 ]
then
  echo -e "Using a gradient of allele effects across the landscape."
  echo -e "Alleles diffusing from the first row of populations."
  ### WILDTYPE
  sed s/BGS.spec/BGS_WT.spec/g ${OUTPREFIX}.ini > ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/QTL.spec/QTL_WT.spec/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"generations $nGen"/"generations 1"/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"quanti_genot_logtime $nGen"/"quanti_genot_logtime 1"/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"quanti_phenot_logtime $nGen"/"quanti_phenot_logtime 1"/g ${OUTPREFIX}_WILDTYPE.ini
  ${QUANTINEMO_DIR}/quantinemo ${OUTPREFIX}_WILDTYPE.ini
  cp ${OUTPREFIX}_*/*.dat WILDTYPE.dat ### P(alleles with effects) = 0.00
  rm -R ${OUTPREFIX}_*/
  ### UNIFORM
  sed s/"generations $nGen"/"generations 1"/g ${OUTPREFIX}.ini > ${OUTPREFIX}_UNIFORM.ini
  sed -i s/"quanti_genot_logtime $nGen"/"quanti_genot_logtime 1"/g ${OUTPREFIX}_UNIFORM.ini
  sed -i s/"quanti_phenot_logtime $nGen"/"quanti_phenot_logtime 1"/g ${OUTPREFIX}_UNIFORM.ini
  ${QUANTINEMO_DIR}/quantinemo ${OUTPREFIX}_UNIFORM.ini
  cp ${OUTPREFIX}_*/*.dat UNIFORM.dat ### P(allele_max) = 1/2N
  rm -R ${OUTPREFIX}_*/
  ### MERGE
  NCOL=$(echo "sqrt ($nPop)" | bc)
  ROW1_NLINES=$(echo "(1 + ${nLoci}) + (${NCOL} * ${nIndividuals})" | bc)
  TAILROW1_NLINES=$(echo "${ROW1_NLINES} + 1" | bc)
  head -n${ROW1_NLINES} UNIFORM.dat > GENOTYPES_INI.dat           ### header +row 1 diffusers
  tail -n+${TAILROW1_NLINES} WILDTYPE.dat >> GENOTYPES_INI.dat   ### wildtype populations
  ### APPEND TO THE INI FILE
  echo -e "quanti_ini_genotypes GENOTYPES_INI.dat" >> ${OUTPREFIX}.ini
elif [ $GRADIENT == 2 ]
then
  echo -e "Using a gradient of allele effects across the landscape."
  echo -e "Alleles diffusing from the first and last rows of populations."
  ### WILDTYPE
  sed s/BGS.spec/BGS_WT.spec/g ${OUTPREFIX}.ini > ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/QTL.spec/QTL_WT.spec/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"generations $nGen"/"generations 1"/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"quanti_genot_logtime $nGen"/"quanti_genot_logtime 1"/g ${OUTPREFIX}_WILDTYPE.ini
  sed -i s/"quanti_phenot_logtime $nGen"/"quanti_phenot_logtime 1"/g ${OUTPREFIX}_WILDTYPE.ini
  ${QUANTINEMO_DIR}/quantinemo ${OUTPREFIX}_WILDTYPE.ini
  cp ${OUTPREFIX}_*/*.dat WILDTYPE.dat ### P(alleles with effects) = 0.00
  rm -R ${OUTPREFIX}_*/
  ### UNIFORM
  sed s/"generations $nGen"/"generations 1"/g ${OUTPREFIX}.ini > ${OUTPREFIX}_UNIFORM.ini
  sed -i s/"quanti_genot_logtime $nGen"/"quanti_genot_logtime 1"/g ${OUTPREFIX}_UNIFORM.ini
  sed -i s/"quanti_phenot_logtime $nGen"/"quanti_phenot_logtime 1"/g ${OUTPREFIX}_UNIFORM.ini
  ${QUANTINEMO_DIR}/quantinemo ${OUTPREFIX}_UNIFORM.ini
  cp ${OUTPREFIX}_*/*.dat UNIFORM.dat ### P(allele_max) = 1/2N
  rm -R ${OUTPREFIX}_*/
  ### MERGE
  NCOL=$(echo "sqrt ($nPop)" | bc)
  ROW1_NLINES=$(echo "(1 + ${nLoci}) + (${NCOL} * ${nIndividuals})" | bc)
  TAILROW1_NLINES=$(echo "${ROW1_NLINES} + 1" | bc)
  TOTAL_NLINES=$(echo "(1 + ${nLoci}) + (${nPop} * ${nIndividuals})" | bc)
  SANDWICH_NLINES=$(echo "$TOTAL_NLINES - $ROW1_NLINES - (${NCOL} * ${nIndividuals})" | bc)
  ROW2_NLINES=$(echo "${NCOL} * ${nIndividuals}" | bc)
  head -n${ROW1_NLINES} UNIFORM.dat > GENOTYPES_INI.dat                                       ### header +row1 diffusers
  tail -n+${TAILROW1_NLINES} WILDTYPE.dat | head -n${SANDWICH_NLINES}  >> GENOTYPES_INI.dat  ### wild type populations
  tail -n${ROW2_NLINES} UNIFORM.dat >> GENOTYPES_INI.dat                                     ### last row diffusers
  ### APPEND TO THE INI FILE
  echo -e "quanti_ini_genotypes GENOTYPES_INI.dat" >> ${OUTPREFIX}.ini
else
  echo -e "Sorry ${GRADIENT} is not a valid gradient input."
  echo -e "   - use 0 for none (i.e. alleles with effects randomly scattered across the landscape)"
  echo -e "   - use 1 for alleles with effects diffusing from the first row of populations"
fi

### clean-up
rm *.temp

#################
### EXECUTION ###
#################
${QUANTINEMO_DIR}/quantinemo ${OUTPREFIX}.ini
echo "#########################################################################################"
