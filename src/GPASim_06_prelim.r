##############################
### ANALYSIS PER LANDSCAPE ### MODULE 1 of 3:
############################## Raw data visualization, naive linear modeling, data streamlining, and preparation for ABC simulations

####################
### user inputs ###
args = commandArgs(trailingOnly=TRUE)
DIR = args[1]   ### GPASIM output directory
PLOT = args[2]  ### save plots?
# ### TEST:
# DIR = "/data/Lolium/Quantitative_Genetics/LOLSIM_2019_TEST/LOLSIM_1rep_5qtl_0.001mr_0.25fgs_0.00bgs_0grad"
# PLOT = TRUE

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
library(RColorBrewer)
library(agricolae)

#############################
### set working directory ###
setwd(DIR)

####################################################################
### load merged cross-validation  output per simulated landscape ###
dat = read.csv("CROSS_VALIDATION_OUTPUT_MERGED.csv")
### for old metric labels i.e. before 20191113
if (length(dat$TRUE_POSITIVE_RATE)==0) {
  dat$TRUE_POSITIVE_RATE = dat$QTL_DETECTED_PERC
  dat$TRUE_POSITIVE_ID = dat$QTL_DETECTED_ID
}

####################################################################
### set the metrics we will be assessing model performances with ###
dat$LOG10_RMSD = log10(dat$RMSD + 1) ### log10-transform RMSD + added 1 to avoid log10(0)
METRICS = list(CORRELATION=dat$CORRELATION, LOG10_RMSD=dat$LOG10_RMSD, TRUE_POSITIVE_RATE=dat$TRUE_POSITIVE_RATE, FALSE_POSITIVE_RATE=dat$FALSE_POSITIVE_RATE)
LABELS = c("Correlation", "log10(Root Mean Square Deviation)", "True Positive Rate", "False Positive Rate")

########################################
### prelimiaries: visualize the data ###
### metrics distributions
colors_fill = brewer.pal(length(METRICS), "Pastel1")
colors_bord = brewer.pal(length(METRICS), "Set1")
if (PLOT==TRUE){svg(paste0(basename(DIR), "-METRICS_DISTRIBUTIONS.svg"), width=10, height=6)}
par(mfrow=c(2,2))
for (i in 1:length(METRICS)){
  # hist(unlist(METRICS[i][!is.na(METRICS[i])]), col=colors_fill[i], border=colors_bord[i], main=LABELS[i], xlab=LABELS[i])
  hist(unlist(METRICS[i]), col=colors_fill[i], border=colors_bord[i], main=LABELS[i], xlab=LABELS[i])
}
if (PLOT==TRUE){dev.off()}
### boxplots
METRICS = c("CORRELATION", "LOG10_RMSD", "TRUE_POSITIVE_RATE", "FALSE_POSITIVE_RATE")
MODELS = unique(paste(dat$MODEL_ITERATION, dat$MODEL_MODEL, dat$MODEL_COVARIATE, sep="-"))
MODELS = MODELS[order(MODELS)]
GROUPING = levels(dat$GROUPING)
ALGORITHM = levels(dat$ALGORITHM)
for (metric in METRICS){
  # metric = "CORRELATION"
  if (PLOT==TRUE){svg(paste0(basename(DIR), "-BOXPLOTS_", metric, ".svg"), width=10, height=6)}
  par(mfrow=c(2,2), mar=c(20,5,3,1), cex=0.5)
  for (grouping in GROUPING){
    for (algorithm in ALGORITHM){
      # grouping = "ACROSS"
      # algorithm = "INDI"
      sub_dat = subset(subset(dat, GROUPING==grouping), ALGORITHM==algorithm,)
      sub_dat = droplevels(sub_dat)
      sub_dat$MODELS = paste(sub_dat$MODEL_ITERATION, sub_dat$MODEL_MODEL, sub_dat$MODEL_COVARIATE, sep="-")
      ncolors = length(unique(sub_dat$MODELS))
      colors_fill = rep(brewer.pal(min(9, ncolors), "Pastel1"), times=ceiling(ncolors/9))
      colors_bord = rep(brewer.pal(min(9, ncolors), "Set1"), times=ceiling(ncolors/9))
      boxplot(eval(parse(text=metric)) ~ MODELS, data=sub_dat, las=2, ylab=metric, main=paste0(paste(metric, grouping, algorithm, sep="-"), "\n(n=", nrow(sub_dat), ")"), col=colors_fill[1:ncolors], border=colors_bord[1:ncolors])
    }
  }
  if (PLOT==TRUE){dev.off()}
}

############################
### linear model fitting ###
GROUP = as.factor(dat$GROUPING)
LEVEL = as.factor(dat$ALGORITHM)
VARIABLE = as.factor(matrix(unlist(strsplit(as.character(dat$MODEL_MODEL), "_")), ncol=2, byrow=TRUE)[,1])
COVARIATE = as.factor(dat$MODEL_COVARIATE)
MODEL = as.factor(paste0(dat$MODEL_ITERATION, "_", matrix(unlist(strsplit(as.character(dat$MODEL_MODEL), "_")), ncol=2, byrow=TRUE)[,2]))
POP_TRAIN = as.factor(dat$POP_TRAIN)
POP_TEST = as.factor(dat$POP_TEST)
CORRELATION = dat$CORRELATION
LOG10_RMSD = dat$LOG10_RMSD
TRUE_POSITIVE_RATE = dat$TRUE_POSITIVE_RATE
FALSE_POSITIVE_RATE = dat$FALSE_POSITIVE_RATE
TRUE_POSITIVE_ID = dat$TRUE_POSITIVE_ID
FALSE_POSITIVE_ID = dat$FALSE_POSITIVE_ID
DATA = data.frame(GROUP, LEVEL, VARIABLE, COVARIATE, MODEL, POP_TRAIN, POP_TEST, CORRELATION, LOG10_RMSD, TRUE_POSITIVE_RATE, FALSE_POSITIVE_RATE, TRUE_POSITIVE_ID, FALSE_POSITIVE_ID)
# ### remove self-cross-validation
# DATA = DATA[(POP_TRAIN != POP_TEST), ]
# DATA = droplevels(DATA)
### build the linear models and perform mean comparisons
METRICS = c("CORRELATION", "LOG10_RMSD", "TRUE_POSITIVE_RATE", "FALSE_POSITIVE_RATE")
FACTORS = c("GROUP", "LEVEL", "VARIABLE", "COVARIATE", "MODEL")
for (i in 1:length(METRICS)){
  # i = 1
  metric = METRICS[i]
  print(metric)
  model = eval(parse(text=paste0("lm(", metric, " ~ ", paste(FACTORS, collapse="+"), ", data=DATA)")))
  if (PLOT==TRUE){svg(paste0(basename(DIR), "-HSD_", metric, ".svg"), width=7, height=10)}
  colors_fill = brewer.pal(length(FACTORS), "Pastel1")
  par(mfrow=c(length(FACTORS), 1), mar=c(5, 13, 2, 2))
  for (j in 1:length(FACTORS)){
    # j = 1
    factor = FACTORS[j]
    hsd_temp = eval(parse(text=paste0("HSD.test(model, trt='", factor, "')$groups")))
    hsd = data.frame(FACTOR=rep(factor, times=nrow(hsd_temp)), LEVEL=rownames(hsd_temp), MEAN=hsd_temp[,1], GROUP=hsd_temp[,2])
    if (!exists("HSD")){
      HSD =  hsd
    } else {
      HSD = rbind(HSD, hsd)
    }
    # y = eval(parse(text=metric))
    # x = eval(parse(text=factor))
    # dy = (max(y, na.rm=T) - min(y, na.rm=T))/(2^4)
    y = hsd$MEAN[order(hsd$MEAN, decreasing=FALSE)]
    # lim = c(max(c(0, min(eval(parse(text=metric)),na.rm=T))), max(eval(parse(text=metric)),na.rm=T))
    lim = c(max(c(0, min(eval(parse(text=metric)),na.rm=T))), max(y,na.rm=T))
    dy = diff(lim)/20
    lim = c(max(c(0, min(eval(parse(text=metric)),na.rm=T))), max(y,na.rm=T)+dy)
    x = hsd$LEVEL[order(y, decreasing=TRUE)]
    z = hsd$GROUP[order(y, decreasing=TRUE)]
    if (j != length(FACTORS)){
      # # plot_hsd = boxplot(y ~ x, horizontal=TRUE, las=2, ylim=c(min(y,na.rm=T), max(y,na.rm=T)+dy), xlab="", main=factor)
      # plot_bar = barplot(height=y, horiz=TRUE, names.arg=x, las=1, xlim=c(min(c(0, min(y,na.rm=T))), max(y,na.rm=T)+dy), main=factor)
      plot_bar = barplot(height=y, horiz=TRUE, names.arg=x, las=1, xlim=lim, main=factor, col=colors_fill[j])
    } else {
      # # plot_hsd = boxplot(y ~ x, horizontal=TRUE, las=2, ylim=c(min(y,na.rm=T), max(y,na.rm=T)+dy), xlab=metric, main=factor)
      # plot_bar = barplot(height=y, horiz=TRUE, names.arg=x, las=1, xlim=c(min(c(0, min(y,na.rm=T))), max(y,na.rm=T)+dy), main=factor, xlab=metric)
      plot_bar = barplot(height=y, horiz=TRUE, names.arg=x, las=1, xlim=lim, main=factor, xlab=metric, col=colors_fill[j])
    }
    # text(x=1.10, y=1:length(hsd$LEVEL), labels=hsd$GROUP[order(hsd$LEVEL)])
    text(x=y+(dy/2), y=plot_bar, labels=z)
  }
  if (PLOT==TRUE){dev.off()}
  write.table(HSD, file=paste0(basename(DIR), "-HSD_", metric, ".csv"), sep=",", row.names=FALSE, quote=FALSE)
  if(!exists("OUT")){
    OUT = HSD
    colnames(OUT) = c("FACTOR", "LEVEL", metric, paste0("GROUP_", metric))
  } else {
    colnames(HSD) = c("FACTOR", "LEVEL", metric, paste0("GROUP_", metric))
    OUT = merge(OUT, HSD, by=c("FACTOR", "LEVEL"), all=TRUE)
  }
  rm(HSD)
}
OUT = cbind(data.frame(REP = rep(REP,times=nrow(OUT)),
                       QTL = rep(QTL,times=nrow(OUT)),
                       MGR = rep(MGR,times=nrow(OUT)),
                       FGS = rep(FGS,times=nrow(OUT)),
                       BGS = rep(BGS,times=nrow(OUT)),
                       GRAD = rep(GRAD,times=nrow(OUT))),
                       OUT)
write.table(OUT, file=paste0("CV_OUTPUT_SUMMARY-", REP, "rep-", QTL, "qtl-", MGR, "mgr-", FGS, "fgs-", BGS, "bgs-", GRAD, "grad.csv"), sep=",", quote=FALSE, row.names=FALSE)
### May need to revise these models and account for the unbalancedness of the dataset e.g. more datapoints for within than across populations as a result of the square-quadrant sampling

########################################
### resource-efficiency optimisation ###
### remove datapoints with missing data (mostly from mixed models convergence failures) and same population cross-validation
filtered_DATA = droplevels(DATA[(complete.cases(DATA)) & (DATA$POP_TRAIN != DATA$POP_TEST), ])
### add fraction of training pops subset of the test populations and vice-versa; also the number of training and testing populations
train_pops = strsplit(as.character(filtered_DATA$POP_TRAIN), ";")
test_pops = strsplit(as.character(filtered_DATA$POP_TEST), ";")
TRAIN_SUB_TEST = c()
TEST_SUB_TRAIN = c()
TRAIN_SIZE = c()
TEST_SIZE = c()
for (i in 1:length(train_pops)){
  TRAIN_SUB_TEST = c(TRAIN_SUB_TEST, sum(train_pops[[i]] %in% test_pops[[i]]) / length(train_pops[[i]]))
  TEST_SUB_TRAIN = c(TEST_SUB_TRAIN, sum(test_pops[[i]] %in% train_pops[[i]]) / length(test_pops[[i]]))
  TRAIN_SIZE = c(TRAIN_SIZE, length(train_pops[[i]]))
  TEST_SIZE = c(TEST_SIZE, length(test_pops[[i]]))
}
### build the streamlined dataframe
STREAMLINED_DATA = data.frame(GROUP=filtered_DATA$GROUP,
                              LEVEL=filtered_DATA$LEVEL,
                              MODEL=paste(filtered_DATA$VARIABLE, filtered_DATA$COVARIATE, filtered_DATA$MODEL, sep="-"),
                              POP_TRAIN=filtered_DATA$POP_TRAIN,
                              POP_TEST=filtered_DATA$POP_TEST,
                              TRAIN_SIZE = TRAIN_SIZE,
                              TEST_SIZE = TEST_SIZE,
                              TRAIN_SUB_TEST=TRAIN_SUB_TEST,
                              TEST_SUB_TRAIN=TEST_SUB_TRAIN,
                              CORRELATION=filtered_DATA$CORRELATION,
                              LOG10_RMSD=filtered_DATA$LOG10_RMSD,
                              TRUE_POSITIVE_RATE=filtered_DATA$TRUE_POSITIVE_RATE,
                              FALSE_POSITIVE_RATE=filtered_DATA$FALSE_POSITIVE_RATE,
                              TRUE_POSITIVE_ID=filtered_DATA$TRUE_POSITIVE_ID,
                              FALSE_POSITIVE_ID=filtered_DATA$FALSE_POSITIVE_ID)
### subsets per path and do some mean comparisons for the models and effect sizes for pop sizes and the fractions of test/train pops subset of train/test pops
STREAMLINED_DATA_SUBSET_LIST = list()
METRICS = c("CORRELATION", "LOG10_RMSD", "TRUE_POSITIVE_RATE", "FALSE_POSITIVE_RATE")
GROUPS = c("ACROSS", "WITHIN")
LEVELS = c("POOL", "INDI")
EXPLANATORY_VARIABLES = c("MODEL", "TRAIN_SIZE", "TEST_SIZE", "TRAIN_SUB_TEST", "TEST_SUB_TRAIN")
for (i in 1:length(GROUPS)){
  for (j in 1:length(LEVELS)){
    # i = 1
    # j = 1
    group = GROUPS[i]
    level = LEVELS[j]
    sub = droplevels(subset(subset(STREAMLINED_DATA, GROUP==group), LEVEL==level))
    eval(parse(text=paste0("STREAMLINED_DATA_SUBSET_LIST$",  group, "_", level, "= sub")))
    for (k in 1:length(METRICS)){
      # k = 1
      metric = METRICS[k]
      print(paste(metric, group, level, sep="-"))
      mod  = lm(eval(parse(text=paste0(metric, "~", paste(EXPLANATORY_VARIABLES, collapse="+")))) , data=sub)
      anova_table = as.data.frame(cbind(rownames(anova(mod)), anova(mod)))
      rownames(anova_table) = NULL; colnames(anova_table) = c("SV", "DF", "SS", "MS", "F", "SIG")
      write.table(anova_table, file=paste0("MEAN_COMPARISON-", metric, "-", group, "-", level, ".csv"), sep=",", row.names=FALSE, quote=FALSE)
      for (explanatory_variable in EXPLANATORY_VARIABLES){
        tryCatch({
          HSD = HSD.test(mod, trt=explanatory_variable, alpha=0.01)$groups ### and then regress the continuous factors to gen the slopes! i.e. for the train/test sizes and fraction of train/test subsets
          HSD = as.data.frame(cbind(rownames(HSD), HSD))
          rownames(HSD) = NULL; colnames(HSD) = c(explanatory_variable, metric, "GROUPING")
          write.table("", file=paste0("MEAN_COMPARISON-", metric, "-", group, "-", level, ".csv"), sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
          write.table(HSD, file=paste0("MEAN_COMPARISON-", metric, "-", group, "-", level, ".csv"), sep=",", row.names=FALSE, quote=FALSE, append=TRUE)
        }, error = function(e){
          print(paste0("No variation for ", explanatory_variable, "!"))
        }
        )
      }
    }
  }
}
write.table(STREAMLINED_DATA_SUBSET_LIST$ACROSS_POOL, file="ACROSS_POOL.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(STREAMLINED_DATA_SUBSET_LIST$ACROSS_INDI, file="ACROSS_INDI.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(STREAMLINED_DATA_SUBSET_LIST$WITHIN_POOL, file="WITHIN_POOL.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(STREAMLINED_DATA_SUBSET_LIST$WITHIN_INDI, file="WITHIN_INDI.csv", sep=",", row.names=FALSE, quote=FALSE)
  # ### define the cost function --> maximizes: correlation & TPR; and minimizes log10(RMSD) & FPR
  # COST_FUNCTION = function(par, nLib=5, nRep=100, metric="TRUE_POSITIVE_RATE", data=STREAMLINED_DATA_SUBSET_LIST, best_models=FALSE){
  #   # par = c(1, 1, 1, 1)
  #   # nLib = 100
  #   # nRep = 10
  #   # metric = "TRUE_POSITIVE_RATE"
  #   # data=STREAMLINED_DATA_SUBSET_LIST
  #   # best_models = TRUE
  #   ### input check: the number of parameters should equal the number of streamlined data list elemenets
  #   # if (length(par) != length(data)){
  #   if ((length(par)+1) != length(data)){
  #     print("ERROR: The number of parameters is not equal to the number of streamlined data list elemenets")
  #     return(0)
  #   } else {
  #     ### scale the parameters to sum to 1
  #     params_scaled = par/sum(par)
  #     # ### testing odd-one-out contrained par optim
  #     # params_scaled = c(par, 1-sum(par))
  #     ### number of libraries to sample for each path
  #     params_counts = round(params_scaled * nLib)
  #     ### select the best model; then random sampling of populations or population-groups
  #     MEAN_METRIC_VALUES = c()
  #     for (j in 1:nRep){ ### replication
  #       metric_values_vec = c()
  #       for (i in 1:length(data)){
  #         # i = 1
  #         sub = data[[i]]
  #         if (best_models==TRUE){
  #           mod = eval(parse(text=paste0("lm(", metric, " ~ MODEL, data=sub)")))
  #           HSD = HSD.test(mod, trt="MODEL")$groups
  #           best_models = rownames(HSD)[HSD$groups == "a"]
  #           sub_best_mods = droplevels(sub[sub$MODEL %in% best_models, ])
  #         } else {
  #           sub_best_mods = sub
  #         }
  #         if (params_counts[i] > 0) {
  #           RAND_SAMP = droplevels(sub_best_mods[sample(x=c(1:nrow(sub_best_mods)), size=params_counts[i], replace=TRUE), ])
  #           if ((metric == "TRUE_POSITIVE_RATE") | (metric == "FALSE_POSITIVE_RATE")){
  #             metric_values_vec = unique(c(metric_values_vec, unlist(strsplit(as.character(eval(parse(text=paste0("RAND_SAMP$", sub("_RATE", "_ID", metric))))), ";"))))
  #           } else {
  #             metric_values_vec = c(metric_values_vec, mean(eval(parse(text=paste0("RAND_SAMP$", metric)))))
  #           }
  #         }
  #       }
  #       if ((metric == "TRUE_POSITIVE_RATE") | (metric == "FALSE_POSITIVE_RATE")){
  #         MEAN_METRIC_VALUES = c(MEAN_METRIC_VALUES, length(unique(metric_values_vec)))
  #       } else {
  #         MEAN_METRIC_VALUES = c(MEAN_METRIC_VALUES, mean(metric_values_vec))
  #       }
  #     }
  #     if ((metric == "TRUE_POSITIVE_RATE") | (metric == "CORRELATION")) {
  #       MINIMIZE_ME = 1 - mean(MEAN_METRIC_VALUES)
  #     } else if ((metric == "FALSE_POSITIVE_RATE") | (metric == "LOG10_RMSD")){
  #       MINIMIZE_ME = mean(MEAN_METRIC_VALUES)
  #     } else {
  #       print("Error metric: ", metric, " invalid!")
  #       print("Please use: 'CORRELATION', 'LOG10_RMSD', 'TRUE_POSITIVE_RATE', or 'FALSE_POSITIVE_RATE'")
  #     }
  #   }
  #   return(MINIMIZE_ME)
  # }
  #
  # ### optimize!
  # OPTIM_PAR_TEMP = c()
  # for (metric in METRICS){
  #   print("###############################################")
  #   print(metric)
  #   out_optim = optim(par=c(1,1,1,1), fn=COST_FUNCTION, method="Nelder-Mead", nLib=5, nRep=1000, metric=metric, data=STREAMLINED_DATA_SUBSET_LIST, best_models=FALSE)
  #
  #   ### test optimx
  #   library(optimx)
  #   # out_optim = opm(par=1:3/9, fn=COST_FUNCTION, method=c("Nelder-Mead", "BFGS"), gr="grcentral", control=list(trace=0), nLib=5, nRep=10, metric=metric, data=STREAMLINED_DATA_SUBSET_LIST, best_models=FALSE)
  #   out_optim = opm(par=rep(1/9, 3), fn=COST_FUNCTION, method=c("Nelder-Mead", "BFGS"), control=list(trace=0), nLib=5, nRep=10, metric=metric, data=STREAMLINED_DATA_SUBSET_LIST, best_models=FALSE)
  #
  #   print(out_optim$convergence)
  #   OPTIM_PAR_TEMP = c(OPTIM_PAR_TEMP, out_optim$par)
  # }
  # OPTIM_PAR = matrix(OPTIM_PAR_TEMP, ncol=4, byrow=TRUE)
  # OPTIM_PAR = OPTIM_PAR / rowSums(OPTIM_PAR)
  # OPTIM_PAR = data.frame(METRICS=METRICS,
  #                       AP=OPTIM_PAR[,1],
  #                       AI=OPTIM_PAR[,2],
  #                       WP=OPTIM_PAR[,3],
  #                       WI=OPTIM_PAR[,4])
  # colnames(OPTIM_PAR) = c("METRICS", names(STREAMLINED_DATA_SUBSET_LIST))
  # write.table(OPTIM_PAR, file=paste0(basename(DIR), "-RESOURCE_OPTIMISATION.csv"), sep=",", row.names=FALSE, quote=FALSE)
  #
  # ### show optimization results as a heatmap
  # ### the redder the higher the proportion (relative to the others) of the resources devoted for the level of experimental setup...
  # ### i.e. across-pool, across-indi, within-pool and within-indi
  # HEATMAP = as.matrix(OPTIM_PAR[,2:ncol(OPTIM_PAR)])
  # rownames(HEATMAP) = OPTIM_PAR[,1]
  # if (PLOT==TRUE){svg(paste0(basename(DIR), "-RESOURCE_OPTIMISATION.svg"), width=10, height=10)}
  # heatmap(HEATMAP, Colv = NA, Rowv = NA, scale="column", margins=c(20, 20))
  # if (PLOT==TRUE){dev.off()}
  #
  # ### visualizing the cost function across 4-dimensions i.e. ap, ai, wp and wi
  # n = 10
  # ARRAY = array(0, dim=c(n, n, n, n))
  # pb = txtProgressBar(min=0, max=(n^4), initial=1, style=3, width=30); counter=1
  # for (i in 1:n){
  #   for (j in 1:n){
  #     for (k in 1:n){
  #       for (l in 1:n){
  #         ARRAY[i,j,k,l] = COST_FUNCTION(par=c(i,j,k,l), nLib=5, nRep=10, metric="TRUE_POSITIVE_RATE", data=STREAMLINED_DATA_SUBSET_LIST, best_models=FALSE)
  #         setTxtProgressBar(pb, counter); counter = counter + 1
  #       }
  #     }
  #   }
  # }
  # close(pb)
  # heatmap(ARRAY[1:10,1:10,1,1], Colv = NA, Rowv = NA, xlab=colnames(OPTIM_PAR)[2], ylab=colnames(OPTIM_PAR)[3])
  # heatmap(ARRAY[1,1,1:10,1:10], Colv = NA, Rowv = NA, xlab=colnames(OPTIM_PAR)[4], ylab=colnames(OPTIM_PAR)[5])
