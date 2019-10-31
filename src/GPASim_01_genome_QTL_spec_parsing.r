#################################################
###                                           ###
### PARSING THE GENOME AND QTL SPECIFICATIONS ###
###                                           ###
#################################################

#############
### INPUT ###
#############
args = commandArgs(trailing=TRUE)
chrom_r = read.delim(args[1], header=FALSE, sep=" ")
genome_spec = read.delim(args[2], header=TRUE, sep="\t")
QTL_spec = read.delim(args[3], skip=8, header=TRUE, sep="\t")
QTL_idx_map_to_genome_spec  = read.delim(args[4], header=FALSE, sep="\t")
# ### test
# chrom_r = read.delim("chrom_r.temp", header=FALSE, sep=" ")
# genome_spec = read.delim("Lperenne_genome.spec", header=TRUE, sep="\t")
# QTL_spec = read.delim("QTL.spec", skip=8, header=TRUE, sep="\t")
# QTL_idx_map_to_genome_spec  = read.delim("QTL_idx_map_to_genome.spec", header=FALSE, sep="\t")

##############
### OUTPUT ###
##############
### (1) "GENOME_SPEC.csv" (HEADER: CHROM, r, POS)
### (2) "QTL_SPEC.csv" (HEADER: CHROM, POS, ALLELE, EFFECT)

########################
### SAMPLE EXECUTION ###
########################
# Rscript GPASim_01_genome_QTL_spec_parsing.r \
#   chrom_r.temp \
#   Lperenne_genome.spec \
#   QTL.spec \
#   QTL_idx_map_to_genome.spec

################################
### calculate position in bp ###
################################
CHROM = rep(chrom_r[,1], each=ncol(chrom_r)-1)
CF = matrix((genome_spec$Mbp * 1.0e6) / genome_spec$cM, ncol=1)
r = matrix(t(chrom_r[, 2:ncol(chrom_r)]), ncol=1)
POS = matrix(t(round(chrom_r[, 2:ncol(chrom_r)] * CF, 0)), ncol=1) #each group of nChrom r is 1 chromosome
OUT_LOCI_SPEC = data.frame(CHROM=CHROM, r=r, POS=POS)
OUT_LOCI_SPEC = OUT_LOCI_SPEC[complete.cases(OUT_LOCI_SPEC),] # when nLoci is not a multiple of nChrom - remove missing r (cM) and POS
write.table(OUT_LOCI_SPEC, "GENOME_SPEC.csv", row.names=FALSE, sep=",")

###############################################################################################
### replace the consecutive locus ID from QTL_spec with QTL_idx_map_to_genome_spec locus ID ###
###############################################################################################
# QTL_spec = QTL_spec[, 1:3] ### retaining only locus ID, allele_ID and allele_effect
# colnames(QTL_spec) = c("ID", "allele", "eff")
# QTL_spec$ID = QTL_idx_map_to_genome_spec[,1]
# idx_ID_non_zero = unique(QTL_spec$ID[QTL_spec$eff != 0.0])
# QTL_spec = QTL_spec[(QTL_spec$ID %in% idx_ID_non_zero), ]
# CHROM = OUT_LOCI_SPEC$CHROM[QTL_spec$ID]
# POS = OUT_LOCI_SPEC$POS[QTL_spec$ID]
# ALLELE = c("A", "T", "C", "G", "DEL")[QTL_spec$allele]
# EFFECT = QTL_spec$eff
# OUT_QTL_SPEC = data.frame(CHROM, POS, ALLELE, EFFECT)
# write.table(OUT_QTL_SPEC, "QTL_SPEC.csv", row.names=FALSE, sep=",")

### OR OR:
QTL_spec = QTL_spec[, 1:3] ### retaining only locus ID, allele_ID and allele_effect
colnames(QTL_spec) = c("LOCI_ID", "ALLELE", "EFFECT")
# QTL_spec$ID = QTL_idx_map_to_genome_spec[,1]

nAlleles = max(QTL_spec$ALLELE)
LOCI_ID = rep(1:nrow(OUT_LOCI_SPEC), each=nAlleles)
CHROM = rep(OUT_LOCI_SPEC$CHROM, each=nAlleles)
POS = rep(OUT_LOCI_SPEC$POS, each=nAlleles)
ALLELE = rep(1:nAlleles, times=nrow(OUT_LOCI_SPEC))
FOR_MERGING = data.frame(LOCI_ID, CHROM, POS, ALLELE)

QTL_SPEC = merge(FOR_MERGING, QTL_spec, by=c("LOCI_ID", "ALLELE"))
QTL_SPEC = QTL_SPEC[order(QTL_SPEC$LOCI_ID), ]
idx = unique(QTL_SPEC$LOCI_ID[QTL_SPEC$EFFECT != 0.0])
OUT_QTL_SPEC = QTL_SPEC[(QTL_SPEC$LOCI_ID %in% idx), c(3,4,2,5)]
OUT_QTL_SPEC$ALLELE = c("A", "T", "C", "G", "DEL", "N")[OUT_QTL_SPEC$ALLELE]
write.table(OUT_QTL_SPEC, "QTL_SPEC.csv", row.names=FALSE, sep=",")
