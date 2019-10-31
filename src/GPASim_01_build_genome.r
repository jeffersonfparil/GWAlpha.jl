###########################
###                     ###
### SIMULATE THE GENOME ###
###                     ###
###########################

#############
### INPUT ###
#############
args = commandArgs(trailing=TRUE)
nChrom = as.numeric(args[1])                                    #INT
chrom_len_cM = as.numeric(unlist(strsplit(args[2], ",")))       #STR: FLOAT1,FLOAT2,FLOAT3, ...FLOATnChrom
nLoci = as.numeric(args[3])                                     #INT
###test
# nChrom = 7
# chrom_len_cM = as.numeric(unlist(strsplit("97.7,151.5,63.3,119.2,89.1,115.2,113.7", ",")))
# nLoci = 1000

########################
### SAMPLE EXECUTION ###
########################
# Rscript ${GEN_PRED_SRC_DIR}/GPWASim_01_build_genome.r \
#   7 \
#   97.7,151.5,63.3,119.2,89.1,115.2,113.7 \
#   $nLoci
#   ### Output: genome_loci_*.temp

##############################################################
### RANDOM SAMPLING OF LOCUS POSITIONS PER CHROMOSOME (cM) ###
##############################################################
nAllelesPerChrom = round(nLoci / nChrom) #may cause the number of positions sampled not equal to nLoci -> fixed by the if-statement below
counter_chrom=1
for (i in chrom_len_cM){
  if (i != chrom_len_cM[nChrom]){ ### testing if we are at the last chromosome
    out = matrix(sort(sample(seq(from=0.00001, to=i, by=0.00001), size=nAllelesPerChrom)), nrow=1)
  } else { ### if we are at the last chromosome sample the correct number of positions so that we have the correct number of loci: nLoci
    out = matrix(sort(sample(seq(from=0.00001, to=i, by=0.00001), size=(nAllelesPerChrom + (nLoci - (nChrom*nAllelesPerChrom))))), nrow=1)
  }
  write.table(out, file=paste0('genome_loci_', counter_chrom, '.temp'), col.names=FALSE, row.names=FALSE, sep=' ') ### write-out the locus positions per chromosome
  counter_chrom = counter_chrom + 1
}
