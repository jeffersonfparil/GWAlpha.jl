####################################
###                              ###
### SIMULATE LOCI SPECIFICATIONS ###
###                              ###
####################################

#############
### INPUT ###
#############
args = commandArgs(trailing=TRUE)
nLoci = as.numeric(args[1])                   ### total number of loci to simulate which will be partitioned into QTL of the trait of interest and QTL of the background-selection trait
nQTL = as.numeric(args[2])                    ### number of QTL for the trait of interest
BGS_nQTL = as.numeric(args[3])                ### number of QTL for the background-selection trait --->>> nLoci = nQTL + BGS_nQTL + zero_effect_loci
allele_dist = read.csv(args[4], header=FALSE) ### frequency table of the number of alleles in a locus
allele_eff_model = args[5]                    ### probability distribution to sample the allele effects from (i.e. CHISQ, NORM, & UNIF)
Ne = as.numeric(args[6])                      ### effective population size
### NOTE: For simplicity the same probability distribution is used for both the trait of interest and the background selection trait
### test
# nLoci = 1000
# nQTL = 10
# BGS_nQTL = 100
# allele_dist = read.csv("Allele_counts_dist.spec", header=FALSE)
# allele_eff_model = "CHISQ"
# Ne = 1000

########################
### SAMPLE EXECUTION ###
########################
# Rscript GPWASim_01_build_loci_spec.r \
#   1000 \
#   10 \
#   100 \
#   Allele_counts_dist.spec \
#   CHISQ
### Output:
  ### (1) "BGS.spec"
  ### (2) "BGS_idx_map_to_genome.spec"
  ### (3) "QTL.spec"
  ### (4) "QTL_idx_map_to_genome.spec"
  ### (5) "QTL_min_max_GEBV.spec"

########################################
### INITIALIZE OUTPUT VECTORS' BASIS ###
########################################
nAlleles = max((1:nrow(allele_dist))[allele_dist[,2] > 1.0e-10])  ### maximum number of alleles per locus
col_locus = rep(1:nLoci, each=nAlleles)                           ### loci indices that maps to genome ---> but should be consecutive numbers for each QTL spec file (trait of interest and BGS trait)
col_allele = rep(1:nAlleles, times=nLoci)
col_mut_freq = rep(c(0), times=nLoci*nAlleles)

#########################################################################################
### GENERATE A SPARSITY MATRIX REFLECTING THE ALLELE COUNT DISTRIBUTION (allele_dist) ###_
#########################################################################################
n_possible_sparsity_designs = 100
sparsity_list = c()
for (i in 2:nAlleles) { #start with 2 alleles - do not simulate monomorphic sites
  for (j in 1:round(allele_dist[i,2]*n_possible_sparsity_designs)) {
    sparse_vec = c(rep(0, time=nAlleles-i), rep(1, times=i))
    sparse_vec = sparse_vec[order(runif(nAlleles))]
    sparsity_list = c(sparsity_list, sparse_vec)
  }
}
sparsity_df = matrix(sparsity_list, ncol=nAlleles, byrow=TRUE) ### nrows of possible allele configurations based on allele counts X ncols of nAlleles
ALLELE_COUNTS_MAT = sparsity_df[sample(1:nrow(sparsity_df), size=nLoci, replace=TRUE), ] ### randomly sample nLoci locus configurations
ALLELE_COUNTS_VEC = matrix(t(ALLELE_COUNTS_MAT), nrow=1) ### reshape into a vector grouped into nAlleles per locus

####################################################
### SET INITIAL ALLELE FREQUENCIES AS 1/NALLELES ###
####################################################
col_ini_freq = rep(1/rowSums(ALLELE_COUNTS_MAT), each=nAlleles) * ALLELE_COUNTS_VEC[1,] ### initial allele frequence as the reciprocal of the number of alleles based on the sampled sparsity matrix

###################################################################
### SET THE INDEX LOCATION OF THE QTL FOR THE TRAIT OF INTEREST ###
###################################################################
idxQTL_raw = sample(1:nLoci, nQTL, replace=FALSE) * nAlleles  ### sample the nAllele-th QTL indices
idxQTL = rep(idxQTL_raw, each=nAlleles) ### intermediate step to...
idxQTL = idxQTL - rep((nAlleles-1):0, times=nQTL) ### trying to represent the index of each allele across all QTL
idxQTL_nonZero_iniFreq = col_ini_freq[idxQTL] != 0 ### indexing the QTL allele indices with no intial allele frequency of alleles that are simply absent in the QTlocus
nQTL_nonZero_iniFreq = sum(idxQTL_nonZero_iniFreq) ### counting the number of alleles across QTL with non-zero allele frequencies (number of alleles persent across QTL)

#######################################################
### SET QTL ALLELE EFFECTS OF THE TRAIT OF INTEREST ###
#######################################################
col_allelic_value = rep(c(0), times=nLoci*nAlleles) ### initial allele effects to zero
qtl_effects = c()
for (i in 1:nQTL){
  if (allele_eff_model == "CHISQ"){
    qtl_effects = c(qtl_effects, round(rchisq(n=nAlleles, df=1),2))
  } else if (allele_eff_model == "NORM"){
    qtl_effects = c(qtl_effects, round(rnorm(n=nAlleles, mean=0, sd=1),2))
  } else if (allele_eff_model == "UNIF"){
    qtl_effects = c(qtl_effects, round(runif(n=nAlleles, min=0, max=5),2))
  } else {
    print(paste0("#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"))
    print(paste0("Invalid allelic effect model: ", allele_eff_model, "!"))
    print(paste0("Please choose from: CHISQ, NORM, and UNIF."))
    print(paste0("Exiting now."))
    print(paste0("#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"))
    quit()
  }
}
qtl_effects = qtl_effects[idxQTL_nonZero_iniFreq]
### set the effects of the minimum-squared-effect alleles to zero to simulate the wildtype allele
QTL_ID = rep(1:nQTL, each=nAlleles)[idxQTL_nonZero_iniFreq]
for (i in unique(QTL_ID)){
  # i = 1
  sub = qtl_effects[QTL_ID == i]              ### subset
  sub[abs(sub) == min(abs(sub))] = 0.00       ### set the allele with the lowest effect to zero
  qtl_effects[QTL_ID == i] = sub              ### replace the effects in the original vector
}
### insert the effects into the output column
# qtl_effects = qtl_effects / sum(qtl_effects) #scale by the sum of QTL effects
col_allelic_value[idxQTL[idxQTL_nonZero_iniFreq]] = qtl_effects

##################################################################################
### SET INITIAL ALLELE FREQUENCY OF LARGEST EFFECT ALLELE PER QTLOCUS TO 1/2Ne ###
##################################################################################
# Ne = 1000
freq_min = 1/(2*Ne)
QTL_ALL_EFF = rep(c(0), times=nQTL*nAlleles)
QTL_ALL_EFF[idxQTL_nonZero_iniFreq] = qtl_effects
QTL_DF = data.frame(LOCUS=rep(1:nQTL, each=nAlleles), ALLELE=rep(1:nAlleles, times=nQTL), EFF=QTL_ALL_EFF)[idxQTL_nonZero_iniFreq, ]
QTL_FREQS = c() ### QTL frequencies setting frequency of the allele with maximum effect to 1/2N
WT_FREQS = c()  ### setting the frequency of the wildtype allele to 1.00
for (i in unique(QTL_DF$LOCUS)) {
  # print(i)
  sub = subset(QTL_DF, LOCUS==i)
  ### set P(allele_max) = 1/2N
  eff_max = max(abs(sub$EFF))
  freqs =  rep(0, times=nrow(sub))
  freqs[eff_max == abs(sub$EFF)] = freq_min
  freqs[eff_max != abs(sub$EFF)] = (1 - freq_min) / (nrow(sub)-1)
  QTL_FREQS = c(QTL_FREQS, freqs)
  ### set P(wild type allele) = 1.00
  freqs =  rep(0, times=nrow(sub))
  freqs[abs(sub$EFF) == 0.0][1] = 1.00
  WT_FREQS = c(WT_FREQS, freqs)
}
col_ini_freq[idxQTL[idxQTL_nonZero_iniFreq]] = QTL_FREQS ### column output for P(allele_max) = 1/Ne, i.e. low freq. of favourable allele
col_ini_freq_WT = col_ini_freq; col_ini_freq_WT[idxQTL[idxQTL_nonZero_iniFreq]] = WT_FREQS ### column output for P(allele == 0) = 1.00, i.e. wildtype population

############################################################################
### SET THE INDEX LOCATION OF THE QTL FOR THE BACKGROUND SELECTION TRAIT ###
############################################################################
loci_list = 1:nLoci * nAlleles ### the vector of indices of all possible loci
BGS_idxQTL_raw = sample(loci_list[!(loci_list %in% idxQTL_raw)], BGS_nQTL, replace=FALSE) ### sample the nAllele-th QTL indices BBBUT only from the indices not already taken by the the QTL of the trait of interest
BGS_idxQTL = rep(BGS_idxQTL_raw, each=nAlleles) ### intermediate step to...
BGS_idxQTL = BGS_idxQTL - rep((nAlleles-1):0, times=BGS_nQTL) ### trying to represent the index of each allele across all QTL
BGS_idxQTL_nonZero_iniFreq = col_ini_freq[BGS_idxQTL] != 0 ### indexing the QTL allele indices with no intial allele frequency of alleles that are simply absent in the QTlocus
BGS_nQTL_nonZero_iniFreq = sum(BGS_idxQTL_nonZero_iniFreq) ### counting the number of alleles across QTL with non-zero allele frequencies (number of alleles persent across QTL)

#########################################################
### SET QTL EFFECTS OF THE BACKGROUND SELECTION TRAIT ###
#########################################################
BGS_col_allelic_value = rep(c(0), times=nLoci*nAlleles) ### initial allele effects to zero
# if        (allele_eff_model == "CHISQ"){
#   ### sample QTL allele effects from a Chi-square distribution at df=1
#   BGS_col_allelic_value[BGS_idxQTL[BGS_idxQTL_nonZero_iniFreq]] = round(rchisq(n=BGS_nQTL_nonZero_iniFreq, df=1),2)
# } else if (allele_eff_model == "NORM"){
#   ### sample QTL allele effects from a Gaussian distribution at mean=0 & sd=1
#   BGS_col_allelic_value[BGS_idxQTL[BGS_idxQTL_nonZero_iniFreq]] = round(rnorm(n=BGS_nQTL_nonZero_iniFreq, mean=0, sd=1),2)
# } else if (allele_eff_model == "UNIF"){
#   ### sample QTL allele effects from a Uniform distribution with range 0 to 5
#   BGS_col_allelic_value[BGS_idxQTL[BGS_idxQTL_nonZero_iniFreq]] = round(runif(n=BGS_nQTL_nonZero_iniFreq, min=0, max=5),2)
# } else {
#   print(paste0("#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"))
#   print(paste0("Invalid allelic effect model: ", allele_eff_model, "!"))
#   print(paste0("Please choose from: CHISQ, NORM, and UNIF."))
#   print(paste0("Exiting now."))
#   print(paste0("#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"))
#   quit()
# }
bgs_effects = c()
for (i in 1:BGS_nQTL){
  if (allele_eff_model == "CHISQ"){
    bgs_effects = c(bgs_effects, round(rchisq(n=nAlleles, df=1),2))
  } else if (allele_eff_model == "NORM"){
    bgs_effects = c(bgs_effects, round(rnorm(n=nAlleles, mean=0, sd=1),2))
  } else if (allele_eff_model == "UNIF"){
    bgs_effects = c(bgs_effects, round(runif(n=nAlleles, min=0, max=5),2))
  } else {
    print(paste0("#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"))
    print(paste0("Invalid allelic effect model: ", allele_eff_model, "!"))
    print(paste0("Please choose from: CHISQ, NORM, and UNIF."))
    print(paste0("Exiting now."))
    print(paste0("#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"))
    quit()
  }
}
bgs_effects = bgs_effects[BGS_idxQTL_nonZero_iniFreq]
### set the effects of the minimum-effect alleles to zero to simulate the wildtype allele
BGS_ID = rep(1:BGS_nQTL, each=nAlleles)[BGS_idxQTL_nonZero_iniFreq]
for (i in unique(BGS_ID)){
  # i = 1
  sub = bgs_effects[BGS_ID == i]              ### subset
  sub[abs(sub) == min(abs(sub))] = 0.00       ### set the allele with the lowest effect to zero
  bgs_effects[BGS_ID == i] = sub              ### replace the effects in the original vector
}
### insert the effects into the output column
BGS_col_allelic_value[BGS_idxQTL[BGS_idxQTL_nonZero_iniFreq]] = bgs_effects

####################################################################################
### SET INITIAL ALLELE FREQUENCY OF LARGEST EFFECT ALLELE PER BGS LOCUS TO 1/2Ne ###
####################################################################################
# Ne = 1000
freq_min = 1/(2*Ne)
BGS_ALL_EFF = rep(c(0), times=BGS_nQTL*nAlleles)
BGS_ALL_EFF[BGS_idxQTL_nonZero_iniFreq] = bgs_effects
BGS_DF = data.frame(LOCUS=rep(1:BGS_nQTL, each=nAlleles), ALLELE=rep(1:nAlleles, times=BGS_nQTL), EFF=BGS_ALL_EFF)[BGS_idxQTL_nonZero_iniFreq, ]
BGS_FREQS = c()
WT_FREQS = c()  ### setting the frequency of the wildtype allele to 1.00
for (i in unique(BGS_DF$LOCUS)) {
  # print(i)
  sub = subset(BGS_DF, LOCUS==i)
  ### set P(allele_max) = 1/2N
  eff_max = max(abs(sub$EFF))
  freqs =  rep(0, times=nrow(sub))
  freqs[eff_max == abs(sub$EFF)] = freq_min
  freqs[eff_max != abs(sub$EFF)] = (1 - freq_min) / (nrow(sub)-1)
  BGS_FREQS = c(BGS_FREQS, freqs)
  ### set P(wild type allele) = 1.00
  freqs =  rep(0, times=nrow(sub))
  freqs[abs(sub$EFF) == 0.0][1] = 1.00
  WT_FREQS = c(WT_FREQS, freqs)
}
col_ini_freq[BGS_idxQTL[BGS_idxQTL_nonZero_iniFreq]] = BGS_FREQS
col_ini_freq_WT[BGS_idxQTL[BGS_idxQTL_nonZero_iniFreq]] = WT_FREQS

####################################################################
### SET THE FINAL QTL INDICES FOR THE BACKGROUND SELECTION TRAIT ###
####################################################################
BGS_loci = col_locus[BGS_idxQTL] ### the indices of QTL for the background selection trait that maps to the genome
COL1 = rep(1:BGS_nQTL, each=nAlleles) ### column1: consecutive numbers of QTL indices
COL2 = col_allele[BGS_idxQTL][order(BGS_loci, decreasing=FALSE)] ### column2 of BGS spec file: allele ID
COL3 = BGS_col_allelic_value[BGS_idxQTL][order(BGS_loci, decreasing=FALSE)] ### column3 of BGS spec file: allele effect
COL4 = col_mut_freq[BGS_idxQTL][order(BGS_loci, decreasing=FALSE)] ### column4 of BGS spec file: mutation rate
COL5 = col_ini_freq[BGS_idxQTL][order(BGS_loci, decreasing=FALSE)] ### column5 of BGS spec file: initial allele frequency
BGS_loci_map_to_genome = BGS_loci[order(BGS_loci, decreasing=FALSE)] ### sorting the QTL indices mapping to the genome
if (length(COL2) == 0){
  COL1 = c()
}
BGS_out = data.frame(COL1, COL2, COL3, COL4, COL5) ### the output BGS spec matrix
write.table(BGS_out, file='BGS.spec', col.names=FALSE, row.names=FALSE, append=TRUE, sep='  ')
write.table(BGS_loci_map_to_genome, file='BGS_idx_map_to_genome.spec', col.names=FALSE, row.names=FALSE, sep='  ')
### wild type population specs
BGS_out_WT = data.frame(COL1, COL2, COL3, COL4, COL5=col_ini_freq_WT[BGS_idxQTL][order(BGS_loci, decreasing=FALSE)])
write.table(BGS_out_WT, file='BGS_WT.spec', col.names=FALSE, row.names=FALSE, append=TRUE, sep='  ')

###########################################################
### SET THE FINAL QTL INDICES FOR THE TRAIT OF INTEREST ###
###########################################################
COL1 = rep(1:(nLoci-BGS_nQTL), each=nAlleles) ### column1: consecutive numbers of locus indices for (nLoci - nBGS) = (nQTL + zero_effect_loci)
COL2 = col_allele[!(col_locus %in% BGS_loci_map_to_genome)] ### column2: allele ID excluding the BGS loci
COL3 = col_allelic_value[!(col_locus %in% BGS_loci_map_to_genome)] ### column3: allele effects excluding the BGS loci
COL4 = col_mut_freq[!(col_locus %in% BGS_loci_map_to_genome)] ### column4: mutation rate excluding BGS loci
COL5 = col_ini_freq[!(col_locus %in% BGS_loci_map_to_genome)] ### column5: initial allele frequencies excluding the BGS loci
if (length(COL2) == 0){
  COL1 = c()
}
out = data.frame(COL1, COL2, COL3, COL4, COL5) # the output QTL spec file including the QTL of the trait of interest and zero-effect loci and EXCLUDING the BGS loci
loci_map_to_genome = col_locus[!(col_locus %in% BGS_loci_map_to_genome)] ### the QTL indices that map to the genome
write.table(out, file='QTL.spec', col.names=FALSE, row.names=FALSE, append=TRUE, sep='	')
write.table(loci_map_to_genome, file='QTL_idx_map_to_genome.spec', col.names=FALSE, row.names=FALSE, sep='	')
### wild type population specs
out_WT = data.frame(COL1, COL2, COL3, COL4, COL5=col_ini_freq_WT[!(col_locus %in% BGS_loci_map_to_genome)])
write.table(out_WT, file='QTL_WT.spec', col.names=FALSE, row.names=FALSE, append=TRUE, sep='  ')

#######################################################################################################
### CALCULATE THE MINIMUM AND MAXIMUM POSSIBLE GENOTYPE BREEDING VALUES BASED ON THE ALLELE EFFECTS ###
#######################################################################################################
sub = out[(out[,3] != 0), ]
MIN_GEBV = 2*sum(aggregate(sub[,3] ~ sub[,1], FUN=min)[,2])
MAX_GEBV = 2*sum(aggregate(sub[,3] ~ sub[,1], FUN=max)[,2])
write.table(data.frame(MIN_GEBV, MAX_GEBV), file='QTL_min_max_GEBV.spec', row.names=FALSE, sep='	')

sub = BGS_out[(BGS_out[,3] != 0), ]
MIN_GEBV = 2*sum(aggregate(sub[,3] ~ sub[,1], FUN=min)[,2])
MAX_GEBV = 2*sum(aggregate(sub[,3] ~ sub[,1], FUN=max)[,2])
write.table(data.frame(MIN_GEBV, MAX_GEBV), file='BGS_min_max_GEBV.spec', row.names=FALSE, sep='	')
