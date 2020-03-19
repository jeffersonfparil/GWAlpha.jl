##############################
### ANALYSIS PER LANDSCAPE ### MODULE 3 of 3:
############################## ABC posterior distribution estimation via aceptance-rejection algorithm

####################
### user inputs ###
args = commandArgs(trailingOnly=TRUE)
DIR = args[1]   ### GPASIM output directory
ABC_OUT_FNAME = args[2]   ### ABC sampling output in csv format
nLIB = as.numeric(args[3]) ### same as the number of librries used in "assessing_CROSS_VALIDATION_OUTPUT_MERGED_2.jl"
PLOT = args[4]  ### save plots?
# ### TEST:
# DIR = getwd()
# ABC_OUT_FNAME = "ABC_SAMPLING_OUTPUT_100REPS.csv"
# nLIB = 10
# PLOT = TRUE

#############################
### set working directory ###
setwd(DIR)

###########################################################
### extract landscape parameters from th directory name ###
REP = as.numeric(sub("rep", "", unlist(strsplit(basename(DIR), "_"))[2]))   ### replication number
QTL = as.numeric(sub("qtl", "", unlist(strsplit(basename(DIR), "_"))[3]))   ### number of QTL
MGR = as.numeric(sub("mr", "", unlist(strsplit(basename(DIR), "_"))[4]))    ### migration rate
FGS = as.numeric(sub("fgs", "", unlist(strsplit(basename(DIR), "_"))[5]))   ### foreground selection intensity
BGS = as.numeric(sub("bgs", "", unlist(strsplit(basename(DIR), "_"))[6]))   ### backegroun selection intensity
GRAD = as.numeric(sub("grad", "", unlist(strsplit(basename(DIR), "_"))[7])) ### resistance alleles landscape gradient

######################
### load libraries ###
library(abc)
library(RColorBrewer)

##########################################
### define observed summary statistics ###
### i.e. metrics we want to optimize: correlation, log10(RMSD +1), true positive rate and false positive rate
TARGET_VECTOR = c(CORRELATION=1.0,
                  LOG10_RMSD=0.0,
                  TRUE_POSITIVE_RATE=1.0,
                  FALSE_DISCOVERY_RATE=0.0)

##########################################
### load the simulated data from Julia ###
### header: REP, ACROSS_POOL, ACROSS_INDI, WITHIN_POOL, WITHIN_INDI, CORRELATION, LOG10_RMSD, TRUE_POSITIVE_RATE, FALSE_DISCOVERY_RATE
# DATA = read.csv("ABC_SAMPLING_OUTPUT_1000REPS.csv")
DATA = read.csv(ABC_OUT_FNAME)
# DATA = read.csv("ABC_SAMPLING_OUTPUT_100REPS.csv")

### bi-plots across the metrics we will be using as summary statistics
if (PLOT==TRUE){png(paste0(basename(DIR), "-ABC_summstats_biplots.png"), width=2000, height=1333)}
METRICS = colnames(DATA)[6:ncol(DATA)]
TARGET_QUADRANT = c("top", "bottom", "top", "bottom")
par(mfrow=c(2, 3), cex=2)
for (i in 1:(length(METRICS)-1)) {
  # i = 1
  x = eval(parse(text=paste0("DATA$", METRICS[i])))
  for (j in (i+1):length(METRICS)) {
    # j = 2
    y = eval(parse(text=paste0("DATA$", METRICS[j])))
    plot(x, y, pch=20, col=rgb(0.25, 0.8, 0.15, alpha=0.1), xlab=METRICS[i], ylab=METRICS[j], main=paste0(METRICS[i], "\nx\n", METRICS[j]))
    if (TARGET_QUADRANT[i] == "top") {
      pos_x = "right"
    } else {
      pos_x = "left"
    }
    pos_y = TARGET_QUADRANT[j]
    legend(paste0(pos_y, pos_x), legend=c("TARGET:", paste0("Quadrant ", pos_y, pos_x)), text.col=rgb(0.7, 0.1, 0.1, alpha=1.0))
  }
}
if (PLOT==TRUE){dev.off()}

### heuristically determine an appropriate tolerance value via odd-one-out cross-validation
start = Sys.time()
cv = cv4abc(param=DATA[,2:5], sumstat=DATA[,6:9], nval=nLIB, tols=c(0.001, 0.01, 0.1), statistic="mean", method="rejection")
end = Sys.time()
print(end - start)
tol = as.numeric(names(which.min(rowMeans(scale(summary(cv))))))
### perform ABC rejection algorithm and estimate the posterior distribution of the parameters
out = abc(target=TARGET_VECTOR, param=DATA[,2:5], sumstat=DATA[,6:9], tol=tol, method="rejection")
DATA_SELECTED = DATA[out$region,]

#####################################################################################################
### plot summary statistics distribution and the posterior distribution of the summary statistics ###
if (PLOT==TRUE){svg(paste0(basename(DIR), "-ABC_resource_optim.svg"), width=6, height=10)}
par(mfrow=c(4,2))
LABELS = c("Across Populations & Pool-Level",
           "Across Populations & Individual-Level",
           "Within Populations & Pool-Level",
           "Within Populations & Individual-Level",
           "Correlation",
           "log10(RMSD + 1)",
           "True Positive Rate",
           "False Positive Rate")
ncolors = length(LABELS)
colors_fill = rep(brewer.pal(min(9, max(3, ncolors)), "Pastel1"), times=ceiling(ncolors/9))
colors_bord = rep(brewer.pal(min(9, max(3, ncolors)), "Set1"), times=ceiling(ncolors/9))
for (i in 1:length(LABELS)){
  hist(DATA_SELECTED[, (i+1)], main=LABELS[i], xlab=colnames(DATA_SELECTED)[i+1], col=colors_fill[i], border=colors_bord[i])
  legend("topright", legend=round(mean(DATA_SELECTED[, (i+1)]),2))
}
if (PLOT==TRUE){dev.off()}

############################
### plot NULL parameters ### equal partitioning
nREP = max(DATA$REP)
NULL_SUBSET = droplevels(DATA[(DATA$ACROSS_POOL==0.25) & (DATA$ACROSS_INDI==0.25) & (DATA$WITHIN_POOL==0.25) & (DATA$WITHIN_INDI==0.25), ])
NULL_CORRELATION = c()
NULL_LOG10_RMSD = c()
NULL_TRUE_POSITIVE_RATE = c()
NULL_FALSE_DISCOVERY_RATE = c()
for (i in 1:nREP){
  idx = sample(1:nrow(NULL_SUBSET), nLIB, replace=FALSE)
  NULL_CORRELATION = c(NULL_CORRELATION, mean(NULL_SUBSET$CORRELATION[idx], na.rm=TRUE))
  NULL_LOG10_RMSD = c(NULL_LOG10_RMSD, mean(NULL_SUBSET$LOG10_RMSD[idx], na.rm=TRUE))
  NULL_TRUE_POSITIVE_RATE = c(NULL_TRUE_POSITIVE_RATE, mean(NULL_SUBSET$TRUE_POSITIVE_RATE[idx], na.rm=TRUE))
  NULL_FALSE_DISCOVERY_RATE = c(NULL_FALSE_DISCOVERY_RATE, mean(NULL_SUBSET$FALSE_DISCOVERY_RATE[idx], na.rm=TRUE))
}
NULL_PARTITION = rep(0.25, times=nREP)
DATA_NULL = data.frame(REP=1:nREP,
                       ACROSS_POOL=NULL_PARTITION,
                       ACROSS_INDI=NULL_PARTITION,
                       WITHIN_POOL=NULL_PARTITION,
                       WITHIN_INDI=NULL_PARTITION,
                       CORRELATION = NULL_CORRELATION,
                       LOG10_RMSD = NULL_LOG10_RMSD,
                       TRUE_POSITIVE_RATE = NULL_TRUE_POSITIVE_RATE,
                       FALSE_DISCOVERY_RATE = NULL_FALSE_DISCOVERY_RATE)
if (PLOT==TRUE){svg(paste0(basename(DIR), "-ABC_NULL.svg"), width=6, height=10)}
par(mfrow=c(4,2))
LABELS = c("Across Populations & Pool-Level",
           "Across Populations & Individual-Level",
           "Within Populations & Pool-Level",
           "Within Populations & Individual-Level",
           "Correlation",
           "log10(RMSD + 1)",
           "True Positive Rate",
           "False Positive Rate")
ncolors = length(LABELS)
colors_fill = rep(brewer.pal(min(9, max(3, ncolors)), "Pastel1"), times=ceiling(ncolors/9))
colors_bord = rep(brewer.pal(min(9, max(3, ncolors)), "Set1"), times=ceiling(ncolors/9))
for (i in 1:length(LABELS)){
  hist(DATA_NULL[, (i+1)], main=LABELS[i], xlab=colnames(DATA_NULL)[i+1], col=colors_fill[i], border=colors_bord[i])
  legend("topright", legend=round(mean(DATA_NULL[, (i+1)]),2))
}
if (PLOT==TRUE){dev.off()}

###########################################################################################################################
### Separate GP and GWAS ABC (using the tolerance value derived from the ABC using both GWAS and GP summary statistics) ###
out_GP   = abc(target=c(1.0, 0.0), param=DATA[,2:5], sumstat=DATA[,6:7], tol=tol, method="rejection")
out_GWAS = abc(target=c(1.0, 0.0), param=DATA[,2:5], sumstat=DATA[,8:9], tol=tol, method="rejection")
summary(out_GP)
summary(out_GWAS)
DATA_SELECTED_GP = DATA[out_GP$region,]
DATA_SELECTED_GWAS = DATA[out_GWAS$region,]

### plot summary statistics distribution and the posterior distribution of the summary statistics ###
if (PLOT==TRUE){svg(paste0(basename(DIR), "-ABC_PER_GP_GWAS.svg"), width=8, height=10)}
par(mfrow=c(6,2))
LABELS_GP = c("Across Populations & Pool-Level",
           "Across Populations & Individual-Level",
           "Within Populations & Pool-Level",
           "Within Populations & Individual-Level",
           "Correlation",
           "log10(RMSD + 1)")
LABELS_GWAS = c("Across Populations & Pool-Level",
          "Across Populations & Individual-Level",
          "Within Populations & Pool-Level",
          "Within Populations & Individual-Level",
          "True Positive Rate",
          "False Positive Rate")
ncolors = length(LABELS_GP) + length(LABELS_GWAS)
colors_fill = rep(brewer.pal(min(9, max(3, ncolors)), "Pastel1"), times=ceiling(ncolors/9))
colors_bord = rep(brewer.pal(min(9, max(3, ncolors)), "Set1"), times=ceiling(ncolors/9))
for (i in 1:length(LABELS_GP)){
  hist(DATA_SELECTED_GP[, (i+1)], main=paste0("GP Optim - ", LABELS_GP[i]), xlab=colnames(DATA_SELECTED_GP)[i+1], col=colors_fill[i], border=colors_bord[i])
  legend("topright", legend=round(mean(DATA_SELECTED_GP[, (i+1)]),2))
}
for (i in 1:length(LABELS_GWAS)){
  hist(DATA_SELECTED_GWAS[, (i+3)], main=paste0("GWAS Optim - ", LABELS_GWAS[i]), xlab=colnames(DATA_SELECTED_GWAS)[i+3], col=colors_fill[i], border=colors_bord[i])
  legend("topright", legend=round(mean(DATA_SELECTED_GWAS[, (i+3)]),2))
}
if (PLOT==TRUE){dev.off()}

######################################################################################################################################################################################
###################
### MAIN OUTPUT ###
###################
### ABC partitioning and metrics summary statistics for each ABC_CLASS testing the:
### (1) NULL hyposthesis,
### (2) optmizing for GP (corr and RMSD) only,
### (3) optimizing for GWAS (tpr and fpr) only, and
### (4) optimizing for both GP and GWAS simultaneously.
OUT = data.frame(REP = rep(REP, each=4),
                 QTL = rep(QTL, each=4),
                 MGR = rep(MGR, each=4),
                 FGS = rep(FGS, each=4),
                 BGS = rep(BGS, each=4),
                 GRAD = rep(GRAD, each=4),
                 ABC_CLASS = c("NULL", "GP", "GWAS", "GPAS"),
                 ACROSS_POOL = c(mean(DATA_NULL$ACROSS_POOL, na.rm=TRUE), mean(DATA_SELECTED_GP$ACROSS_POOL, na.rm=TRUE), mean(DATA_SELECTED_GWAS$ACROSS_POOL, na.rm=TRUE), mean(DATA_SELECTED$ACROSS_POOL, na.rm=TRUE)),
                 ACROSS_INDI = c(mean(DATA_NULL$ACROSS_INDI, na.rm=TRUE), mean(DATA_SELECTED_GP$ACROSS_INDI, na.rm=TRUE), mean(DATA_SELECTED_GWAS$ACROSS_INDI, na.rm=TRUE), mean(DATA_SELECTED$ACROSS_INDI, na.rm=TRUE)),
                 WITHIN_POOL = c(mean(DATA_NULL$WITHIN_POOL, na.rm=TRUE), mean(DATA_SELECTED_GP$WITHIN_POOL, na.rm=TRUE), mean(DATA_SELECTED_GWAS$WITHIN_POOL, na.rm=TRUE), mean(DATA_SELECTED$WITHIN_POOL, na.rm=TRUE)),
                 WITHIN_INDI = c(mean(DATA_NULL$WITHIN_INDI, na.rm=TRUE), mean(DATA_SELECTED_GP$WITHIN_INDI, na.rm=TRUE), mean(DATA_SELECTED_GWAS$WITHIN_INDI, na.rm=TRUE), mean(DATA_SELECTED$WITHIN_INDI, na.rm=TRUE)),
                 CORRELATION = c(mean(DATA_NULL$CORRELATION, na.rm=TRUE), mean(DATA_SELECTED_GP$CORRELATION, na.rm=TRUE), mean(DATA_SELECTED_GWAS$CORRELATION, na.rm=TRUE), mean(DATA_SELECTED$CORRELATION, na.rm=TRUE)),
                 LOG10_RMSD = c(mean(DATA_NULL$LOG10_RMSD, na.rm=TRUE), mean(DATA_SELECTED_GP$LOG10_RMSD, na.rm=TRUE), mean(DATA_SELECTED_GWAS$LOG10_RMSD, na.rm=TRUE), mean(DATA_SELECTED$LOG10_RMSD, na.rm=TRUE)),
                 TRUE_POSITIVE_RATE = c(mean(DATA_NULL$TRUE_POSITIVE_RATE, na.rm=TRUE), mean(DATA_SELECTED_GP$TRUE_POSITIVE_RATE, na.rm=TRUE), mean(DATA_SELECTED_GWAS$TRUE_POSITIVE_RATE, na.rm=TRUE), mean(DATA_SELECTED$TRUE_POSITIVE_RATE, na.rm=TRUE)),
                 FALSE_DISCOVERY_RATE = c(mean(DATA_NULL$FALSE_DISCOVERY_RATE, na.rm=TRUE), mean(DATA_SELECTED_GP$FALSE_DISCOVERY_RATE, na.rm=TRUE), mean(DATA_SELECTED_GWAS$FALSE_DISCOVERY_RATE, na.rm=TRUE), mean(DATA_SELECTED$FALSE_DISCOVERY_RATE, na.rm=TRUE)))
write.table(OUT, file=paste0("ABC_RESOURCE_OPTIM-", REP, "rep-", QTL, "qtl-", MGR, "mgr-", FGS, "fgs-", BGS, "bgs-", GRAD, "grad.csv"), sep=",", quote=FALSE, row.names=FALSE)
######################################################################################################################################################################################

# #####################
# ### MISCELLANEOUS ###
# #####################
# ####################################
# ### plot slected parameter space ###
# par(mfrow=c(2,3))
# for (x in 1:3){
#   for (y in (x+1):4){
#     plot(x=DATA_SELECTED[, x+1], y=DATA_SELECTED[, y+1], xlab=LABELS[x], ylab=LABELS[y], main="", type="p", pch=20, col=c(colors_fill[x], colors_fill[y]))
#   }
# }
#
# library(rgl)
# plot3d(
#   x=DATA_SELECTED$WITHIN_POOL,
#   y=DATA_SELECTED$WITHIN_INDI,
#   z=DATA_SELECTED$TRUE_POSITIVE_RATE,
#   col=rgb(0.1, 0.9, 0.5, alpha=0.5)
# )
# ##################################
# ### lanscape-based model plots ###

### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# ##########################################################################
# ### USING ONLY LOG10_RMSD and TRUE_POSITIVE_RATE AS SUMMARY STATISTICS ### 20200104
# ##########################################################################

#####################################################################
### USING ONLY LOG10_RMSD and TPR/(TPR+FDR) AS SUMMARY STATISTICS ### 20200311
#####################################################################
### define the new summary statistic TPR/(TPR+FDR)
DATA$TPR_OVER_TPR_PLUS_FDR = DATA$TRUE_POSITIVE_RATE / (DATA$TRUE_POSITIVE_RATE + DATA$FALSE_DISCOVERY_RATE)
### heuristically determine an appropriate tolerance value via odd-one-out cross-validation
start = Sys.time()
### columns 7 and 8 are LOG10_RMSD and TRUE_POSITIVE_RATE, respectively
# cv = cv4abc(param=DATA[,2:5], sumstat=DATA[,7:8], nval=nLIB, tols=c(0.001, 0.01, 0.1), statistic="mean", method="rejection")
cv = cv4abc(param=DATA[,2:5], sumstat=DATA[,c(7,10)], nval=nLIB, tols=c(0.001, 0.01, 0.1), statistic="mean", method="rejection")
end = Sys.time()
print(end - start)
tol = as.numeric(names(which.min(rowMeans(scale(summary(cv))))))
### perform ABC rejection algorithm and estimate the posterior distribution of the parameters
# out = abc(target=c(0.0, 1.0), param=DATA[,2:5], sumstat=DATA[,7:8], tol=tol, method="rejection")
out = abc(target=c(0.0, 1.0), param=DATA[,2:5], sumstat=DATA[,c(7,10)], tol=tol, method="rejection")
# DATA_SELECTED = DATA[out$region, c(1,2,3,4,5,7,8)]
DATA_SELECTED = DATA[out$region, c(1,2,3,4,5,7,10)]
#####################################################################################################
### plot summary statistics distribution and the posterior distribution of the summary statistics ###
if (PLOT==TRUE){svg(paste0(basename(DIR), "-ABC_resource_optim_JUST2SUMMSTATS.svg"), width=6, height=10)}
par(mfrow=c(3,2))
LABELS = c("Across Populations & Pool-Level",
           "Across Populations & Individual-Level",
           "Within Populations & Pool-Level",
           "Within Populations & Individual-Level",
           "log10(RMSD + 1)",
           "True Positive Rate")
ncolors = length(LABELS)
colors_fill = rep(brewer.pal(min(9, max(3, ncolors)), "Pastel1"), times=ceiling(ncolors/9))
colors_bord = rep(brewer.pal(min(9, max(3, ncolors)), "Set1"), times=ceiling(ncolors/9))
for (i in 1:length(LABELS)){
  hist(DATA_SELECTED[, (i+1)], main=LABELS[i], xlab=colnames(DATA_SELECTED)[i+1], col=colors_fill[i], border=colors_bord[i])
  legend("topright", legend=round(mean(DATA_SELECTED[, (i+1)]),2))
}
if (PLOT==TRUE){dev.off()}
############################
### plot NULL parameters ### equal partitioning
nREP = max(DATA$REP)
NULL_SUBSET = droplevels(DATA[(DATA$ACROSS_POOL==0.25) & (DATA$ACROSS_INDI==0.25) & (DATA$WITHIN_POOL==0.25) & (DATA$WITHIN_INDI==0.25), ])
NULL_CORRELATION = c()
NULL_LOG10_RMSD = c()
NULL_TRUE_POSITIVE_RATE = c()
NULL_FALSE_DISCOVERY_RATE = c()
NULL_TPR_OVER_TPR_PLUS_FDR = c()
for (i in 1:nREP){
  idx = sample(1:nrow(NULL_SUBSET), nLIB, replace=FALSE)
  NULL_LOG10_RMSD = c(NULL_LOG10_RMSD, mean(NULL_SUBSET$LOG10_RMSD[idx], na.rm=TRUE))
  # NULL_TRUE_POSITIVE_RATE = c(NULL_TRUE_POSITIVE_RATE, mean(NULL_SUBSET$TRUE_POSITIVE_RATE[idx], na.rm=TRUE))
  NULL_TPR_OVER_TPR_PLUS_FDR = c(NULL_TPR_OVER_TPR_PLUS_FDR, mean(NULL_SUBSET$TPR_OVER_TPR_PLUS_FDR[idx], na.rm=TRUE))
}
NULL_PARTITION = rep(0.25, times=nREP)
DATA_NULL = data.frame(REP=1:nREP,
                       ACROSS_POOL=NULL_PARTITION,
                       ACROSS_INDI=NULL_PARTITION,
                       WITHIN_POOL=NULL_PARTITION,
                       WITHIN_INDI=NULL_PARTITION,
                       LOG10_RMSD = NULL_LOG10_RMSD,
                       TPR_OVER_TPR_PLUS_FDR = NULL_TPR_OVER_TPR_PLUS_FDR)
                      #  TRUE_POSITIVE_RATE = NULL_TRUE_POSITIVE_RATE)
if (PLOT==TRUE){svg(paste0(basename(DIR), "-ABC_NULL_JUST2SUMMSTATS.svg"), width=6, height=10)}
par(mfrow=c(3,2))
LABELS = c("Across Populations & Pool-Level",
           "Across Populations & Individual-Level",
           "Within Populations & Pool-Level",
           "Within Populations & Individual-Level",
           "log10(RMSD + 1)",
           "True Positive Rate\n/(True Positive Rate + False Discovery Rate)")
          #  "True Positive Rate")
ncolors = length(LABELS)
colors_fill = rep(brewer.pal(min(9, max(3, ncolors)), "Pastel1"), times=ceiling(ncolors/9))
colors_bord = rep(brewer.pal(min(9, max(3, ncolors)), "Set1"), times=ceiling(ncolors/9))
for (i in 1:length(LABELS)){
  hist(DATA_NULL[, (i+1)], main=LABELS[i], xlab=colnames(DATA_NULL)[i+1], col=colors_fill[i], border=colors_bord[i])
  legend("topright", legend=round(mean(DATA_NULL[, (i+1)]),2))
}
if (PLOT==TRUE){dev.off()}
###########################################################################################################################
### Separate GP and GWAS ABC (using the tolerance value derived from the ABC using both GWAS and GP summary statistics) ###
out_GP   = abc(target=c(0.0), param=DATA[,2:5], sumstat=DATA[,7], tol=tol, method="rejection")
out_GWAS = abc(target=c(1.0), param=DATA[,2:5], sumstat=DATA[,8], tol=tol, method="rejection")
summary(out_GP)
summary(out_GWAS)
# DATA_SELECTED_GP = DATA[out_GP$region, c(1,2,3,4,5,7,8)]
# DATA_SELECTED_GWAS = DATA[out_GWAS$region, c(1,2,3,4,5,7,8)]
DATA_SELECTED_GP = DATA[out_GP$region, c(1,2,3,4,5,7,10)]
DATA_SELECTED_GWAS = DATA[out_GWAS$region, c(1,2,3,4,5,7,10)]

### plot summary statistics distribution and the posterior distribution of the summary statistics ###
if (PLOT==TRUE){svg(paste0(basename(DIR), "-ABC_PER_GP_GWAS_JUST2SUMMSTATS.svg"), width=8, height=10)}
par(mfrow=c(6,2))
LABELS_GP = c("Across Populations & Pool-Level",
           "Across Populations & Individual-Level",
           "Within Populations & Pool-Level",
           "Within Populations & Individual-Level",
           "log10(RMSD + 1)",
           "True Positive Rate")
LABELS_GWAS = c("Across Populations & Pool-Level",
          "Across Populations & Individual-Level",
          "Within Populations & Pool-Level",
          "Within Populations & Individual-Level",
          "log10(RMSD + 1)",
          "True Positive Rate")
ncolors = length(LABELS_GP) + length(LABELS_GWAS)
colors_fill = rep(brewer.pal(min(9, max(3, ncolors)), "Pastel1"), times=ceiling(ncolors/9))
colors_bord = rep(brewer.pal(min(9, max(3, ncolors)), "Set1"), times=ceiling(ncolors/9))
for (i in 1:length(LABELS_GP)){
  hist(DATA_SELECTED_GP[, (i+1)], main=paste0("GP Optim - ", LABELS_GP[i]), xlab=colnames(DATA_SELECTED_GP)[i+1], col=colors_fill[i], border=colors_bord[i])
  legend("topright", legend=round(mean(DATA_SELECTED_GP[, (i+1)]),2))
}
for (i in 1:length(LABELS_GWAS)){
  hist(DATA_SELECTED_GWAS[, (i+1)], main=paste0("GWAS Optim - ", LABELS_GWAS[i]), xlab=colnames(DATA_SELECTED_GWAS)[i+1], col=colors_fill[i], border=colors_bord[i])
  legend("topright", legend=round(mean(DATA_SELECTED_GWAS[, (i+1)]),2))
}
if (PLOT==TRUE){dev.off()}
######################################################################################################################################################################################
###################
### MAIN OUTPUT ###
###################
### ABC partitioning and metrics summary statistics for each ABC_CLASS testing the:
### (1) NULL hyposthesis,
### (2) optmizing for GP (corr and RMSD) only,
### (3) optimizing for GWAS (tpr and fpr) only, and
### (4) optimizing for both GP and GWAS simultaneously.
OUT = data.frame(REP = rep(REP, each=4),
                 QTL = rep(QTL, each=4),
                 MGR = rep(MGR, each=4),
                 FGS = rep(FGS, each=4),
                 BGS = rep(BGS, each=4),
                 GRAD = rep(GRAD, each=4),
                 ABC_CLASS = c("NULL", "GP", "GWAS", "GPAS"),
                 ACROSS_POOL = c(mean(DATA_NULL$ACROSS_POOL, na.rm=TRUE), mean(DATA_SELECTED_GP$ACROSS_POOL, na.rm=TRUE), mean(DATA_SELECTED_GWAS$ACROSS_POOL, na.rm=TRUE), mean(DATA_SELECTED$ACROSS_POOL, na.rm=TRUE)),
                 ACROSS_INDI = c(mean(DATA_NULL$ACROSS_INDI, na.rm=TRUE), mean(DATA_SELECTED_GP$ACROSS_INDI, na.rm=TRUE), mean(DATA_SELECTED_GWAS$ACROSS_INDI, na.rm=TRUE), mean(DATA_SELECTED$ACROSS_INDI, na.rm=TRUE)),
                 WITHIN_POOL = c(mean(DATA_NULL$WITHIN_POOL, na.rm=TRUE), mean(DATA_SELECTED_GP$WITHIN_POOL, na.rm=TRUE), mean(DATA_SELECTED_GWAS$WITHIN_POOL, na.rm=TRUE), mean(DATA_SELECTED$WITHIN_POOL, na.rm=TRUE)),
                 WITHIN_INDI = c(mean(DATA_NULL$WITHIN_INDI, na.rm=TRUE), mean(DATA_SELECTED_GP$WITHIN_INDI, na.rm=TRUE), mean(DATA_SELECTED_GWAS$WITHIN_INDI, na.rm=TRUE), mean(DATA_SELECTED$WITHIN_INDI, na.rm=TRUE)),
                 LOG10_RMSD = c(mean(DATA_NULL$LOG10_RMSD, na.rm=TRUE), mean(DATA_SELECTED_GP$LOG10_RMSD, na.rm=TRUE), mean(DATA_SELECTED_GWAS$LOG10_RMSD, na.rm=TRUE), mean(DATA_SELECTED$LOG10_RMSD, na.rm=TRUE)),
                 TPR_OVER_TPR_PLUS_FDR = c(mean(DATA_NULL$TPR_OVER_TPR_PLUS_FDR, na.rm=TRUE), mean(DATA_SELECTED_GP$TPR_OVER_TPR_PLUS_FDR, na.rm=TRUE), mean(DATA_SELECTED_GWAS$TPR_OVER_TPR_PLUS_FDR, na.rm=TRUE), mean(DATA_SELECTED$TPR_OVER_TPR_PLUS_FDR, na.rm=TRUE)))
                #  TRUE_POSITIVE_RATE = c(mean(DATA_NULL$TRUE_POSITIVE_RATE, na.rm=TRUE), mean(DATA_SELECTED_GP$TRUE_POSITIVE_RATE, na.rm=TRUE), mean(DATA_SELECTED_GWAS$TRUE_POSITIVE_RATE, na.rm=TRUE), mean(DATA_SELECTED$TRUE_POSITIVE_RATE, na.rm=TRUE)))
write.table(OUT, file=paste0("ABC_RESOURCE_OPTIM-", REP, "rep-", QTL, "qtl-", MGR, "mgr-", FGS, "fgs-", BGS, "bgs-", GRAD, "grad_JUST2SUMMSTATS.csv"), sep=",", quote=FALSE, row.names=FALSE)
