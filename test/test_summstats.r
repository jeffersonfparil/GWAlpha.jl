args = commandArgs(trailing=TRUE)
fname=args[1]
generation=args[2] ### e.g. g01000 and g10000
window_size=as.numeric(args[3])

### FUNCTION DEFINITION: finding the best dimensions of the allele frequency change through generations plots
FACTORIZE = function(n){
  # n = nlevels(QTL_FREQ_DF$LOCI)
  factor_1 = c()
  factor_2 = c()
  for (i in seq(from=1, to=n)){
    if ((n %% i) == 0){
      factor_1 = c(factor_1, i)
      factor_2 = c(factor_2, n/i)
    }
  }
  diff = abs(factor_1-factor_2)
  nrow=factor_1[c(1:length(diff))[diff==min(diff)][1]]
  ncol=factor_2[c(1:length(diff))[diff==min(diff)][1]]
  return(c(nrow, ncol))
}

dat = read.table(fname, header=TRUE)
# dat = read.table("ALL_POP_NPSTAT.stats", header=TRUE)
# generation = "g500"
# window_size=100000 #100 kb

#### Extract window id of the QTL
qtl = read.csv("../QTL_SPEC.csv")
genome = read.csv("../GENOME_SPEC.csv")
QTL_WINDOW_LIST=c()
for (i in unique(genome$CHROM)){
  print(i)
  sub_qtl = subset(qtl, CHROM==i)
  if (nrow(sub_qtl) != 0) {
      sub_genome = subset(genome, CHROM==i)
      n_windows = ceiling(max(sub_genome$POS) / window_size)
      POS=1:max(sub_genome$POS)
      WINDOW=rep(c(1:n_windows), each=window_size)[POS]
      windows_genome = WINDOW[POS %in% sub_genome$POS]
      windows_qtl = rep(windows_genome[sub_genome$POS %in% sub_qtl$POS], each=5)
      QTL_WINDOW_LIST = c(QTL_WINDOW_LIST, windows_qtl)
      sub_dat = max(subset(dat, chrom==i)$window)
    }
}
qtl$WINDOW = QTL_WINDOW_LIST

### Plot Tajima's D across the genome and per chromosome
svg("Tajimas_D_across_populations.svg", width=5, height=7)
# AGGREGATE = aggregate(Tajimas_D ~ window + chrom, FUN=mean, na.action=na.omit, data=subset(subset(dat, pop==population), gen==generation))
AGGREGATE_ALLPOP = aggregate(Tajimas_D ~ window + chrom, FUN=mean, na.rm=TRUE, data=dat)
nrow = FACTORIZE(n=length(unique(AGGREGATE_ALLPOP$chrom))+1)[2]
ncol = FACTORIZE(n=length(unique(AGGREGATE_ALLPOP$chrom))+1)[1]
par(mfrow=c(nrow, ncol))
plot(AGGREGATE_ALLPOP$Tajimas_D, type="p", pch=19, col=rgb(0.1, 0.9, 0.2, alpha=0.5), main=paste0("WHOLE GENOME"), xlab="window", ylab="Tajima's D")
y_lim=c(min(AGGREGATE_ALLPOP$Tajimas_D), max(AGGREGATE_ALLPOP$Tajimas_D))
for (i in unique(AGGREGATE_ALLPOP$chrom)) {
  sub = subset(AGGREGATE_ALLPOP, chrom==i)
  plot(x=c(1, max(sub$window)), y=y_lim, type="n", ylim=y_lim, main=paste0("CHROM_", i), xlab="window", ylab="Tajima's D")
  lines(x=sub$window, y=sub$Tajimas_D, lty=1, col="black")
  abline(v=unique(subset(qtl, CHROM==i)$WINDOW), col="red")
}
dev.off()

### Plot Tajima's D across the genome and per chromosome PER POPULATION
for (population in levels(dat$pop)){
  # population = levels(dat$pop)[1]
  svg(paste0("Tajimas_D_", population,".svg"), width=5, height=7)
  AGGREGATE = aggregate(Tajimas_D ~ window + chrom, FUN=mean, na.action=na.omit, data=subset(subset(dat, pop==population), gen==generation))
  # AGGREGATE = aggregate(Tajimas_D ~ window + chrom, FUN=mean, na.action=na.omit, data=dat)
  AGGREGATE_ALLPOP = aggregate(Tajimas_D ~ window + chrom, FUN=mean, na.rm=TRUE, data=dat)
  nrow = FACTORIZE(n=length(unique(AGGREGATE_ALLPOP$chrom))+1)[2]
  ncol = FACTORIZE(n=length(unique(AGGREGATE_ALLPOP$chrom))+1)[1]
  par(mfrow=c(nrow, ncol))
  plot(AGGREGATE$Tajimas_D, type="p", pch=19, col=rgb(0.1, 0.9, 0.2, alpha=0.5), main=paste0("WHOLE GENOME"), xlab="window", ylab="Tajima's D")
  y_lim=c(min(AGGREGATE$Tajimas_D), max(AGGREGATE$Tajimas_D))
  for (i in unique(AGGREGATE$chrom)) {
    sub = subset(AGGREGATE, chrom==i)
    plot(x=c(1, max(sub$window)), y=y_lim, type="n", ylim=y_lim, main=paste0("CHROM_", i), xlab="window", ylab="Tajima's D")
    lines(x=sub$window, y=sub$Tajimas_D, lty=1, col="black")
    abline(v=unique(subset(qtl, CHROM==i)$WINDOW), col="red")
  }
  dev.off()
}

### Manhattan plot based no the distribution of Tajima's D across the genome

### set x-axis by window number idx
svg("Tajimas_D_LOD_Manhattan_plot.svg", width=15, height=4)
pval_vec = pnorm(q=AGGREGATE_ALLPOP$Tajimas_D, mean=0, sd=sd(AGGREGATE_ALLPOP$Tajimas_D))
lod_vec = -log(pval_vec, base=10)
STATS_DF = data.frame(CHROM=AGGREGATE_ALLPOP$chrom, WINDOW=AGGREGATE_ALLPOP$window, PVAL=pval_vec, LOD=lod_vec)
COLORS = rainbow(length(unique(AGGREGATE_ALLPOP$chrom)), alpha=0.5)
nrow = FACTORIZE(n=length(unique(AGGREGATE_ALLPOP$chrom)))[1]
ncol = FACTORIZE(n=length(unique(AGGREGATE_ALLPOP$chrom)))[2]
par(mfrow=c(nrow, ncol))
for (i in 1:length(unique(STATS_DF$CHROM))){
  chromosome = unique(STATS_DF$CHROM)[i]
  sub = subset(STATS_DF, CHROM==chromosome)
  plot(x=sub$WINDOW, y=sub$LOD, type="p", pch=19, col=COLORS[i], ylim=c(min(STATS_DF$LOD), max(STATS_DF$LOD)), xlab="Locus ID", ylab="-log10(p value)", main=paste0("Chromosome: ", chromosome))
  abline(v=unique(subset(qtl, CHROM==chromosome)$WINDOW), col="red")
}
dev.off()


### extract QTL allele frequencies per population across generations
QTL_FREQS = c()
GEN = c()
POP = c()
CHROM = c()
POS = c()
ALLELE = c()
EFF = c()
FNAMES_GENO_LIST = system("ls MISC/*_GENO.csv | grep -v POOL", intern=TRUE)
pb = txtProgressBar(min=0, max=length(FNAMES_GENO_LIST), initial=0, style=3); counter=1
for (i in 1:length(FNAMES_GENO_LIST)){
  # i = 1
  fname=FNAMES_GENO_LIST[i]
  fname_split = unlist(strsplit(fname, "_"))
  pop = fname_split[length(fname_split)-1]
  gen = fname_split[length(fname_split)-2]
  X = read.csv(fname, header=FALSE)
  freqs = rowMeans(X[,4:ncol(X)])/2
  qtl_freqs = freqs[(X[,1] %in% qtl$CHROM) & (X[,2] %in% qtl$POS)]
  QTL_FREQS = c(QTL_FREQS, qtl_freqs)
  GEN = c(GEN, rep(gen, times=length(qtl_freqs)))
  POP = c(POP, rep(pop, times=length(qtl_freqs)))
  CHROM = c(CHROM, qtl$CHROM)
  POS = c(POS, qtl$POS)
  ALLELE= c(ALLELE, as.character(qtl$ALLELE))
  EFF = c(EFF, qtl$EFFECT)
  setTxtProgressBar(pb, i)
}
close(pb)
GEN = as.numeric(gsub("^g", "", GEN))
LOCI=rep(rep(paste0(qtl$CHROM, "_", qtl$POS), times=4), times=10)
QTL_FREQ_DF = data.frame(GEN=GEN, POP=POP, CHROM=CHROM, POS=POS, LOCI=LOCI, ALLELE=ALLELE, EFF=EFF, FREQ=QTL_FREQS)
LOCI_TO_TRACK = unique(QTL_FREQ_DF$LOCI)
### plot
COLORS_ALLELES = rainbow(nlevels(QTL_FREQ_DF$ALLELE), alpha=1.0)
for (pop in levels(QTL_FREQ_DF$POP)){
  # pop = levels(QTL_FREQ_DF$POP)[1]
  print(pop)
  sub = subset(QTL_FREQ_DF, POP==pop)
  nrow = FACTORIZE(n=nlevels(QTL_FREQ_DF$LOCI))[1]
  ncol = FACTORIZE(n=nlevels(QTL_FREQ_DF$LOCI))[2]
  svg(paste0("Allele_Frequency_Across_Generation_of_Non_Zero_Effect_Alleles_", pop, ".svg"), width=10, height=7)
  par(mfrow=c(nrow,ncol))
  for (j in levels(sub$LOCI)){
    # j = levels(sub$LOCI)[1]
    sub_locus = subset(sub, LOCI==j)
    plot(x=c(min(sub$GEN), max(sub$GEN)), y=c(1, 0), type="n", main=paste0("LOCUS ", j), xlab="Generation", ylab="Allele Frequency")
    for (k in 1:nlevels(sub$ALLELE)) {
      # k = 1
      allele = levels(sub$ALLELE)[k]
      color = COLORS_ALLELES[k]
      sub_allele = subset(sub_locus, ALLELE==allele)
      if (mean(sub_allele$EFF) != 0.0){
        lines(x=sub_allele$GEN, y=sub_allele$FREQ, col=color)
        text(x=mean(sub_allele$GEN), y=mean(sub_allele$FREQ), labels=paste0(allele, "=", round(sub_allele$EFF[1], 2)))
      }
    }
  }
  dev.off()
}
