########################
### Linear modelling ###
########################

### (1) Semi-Exhaustive pairwise mixed models

# ### testing in spartan with beefy server:
# # sinteractive --x11 -c 1 --mem=500000 --partition=bigmem  --time=0-8:00:00
# # sinteractive --x11 -c 1 --mem=120000 --partition=physical --time=0-8:00:00
# sinteractive --x11 -c 32 --mem=120000 --partition=snowy --time=0-8:00:00
# sinteractive --x11 -c 12 --mem=100000 --partition=mig --time=7-0:00:00 --account=punim0594
# module load parallel/20181222-spartan_gcc-6.2.0.lua
# module load Julia/1.1.1-spartan_gcc-6.2.0.lua
# module load R/3.5.2-GCC-6.2.0
# module load GSL/2.5-intel-2018.u4
# module load tmux/2.9a-GCC-6.2.0
# module load gnuplot/5.0.0-GCC-6.2.0
# module load Python
# tmux new

### In R
#####################################################################################################
#################
### LOAD DATA ###
#################
args = commandArgs(trailingOnly=TRUE) ### csv file
# args = c("MERGED.csv")
time_0 = Sys.time()
dat=read.csv(args[1])
time_1 = Sys.time()
print(time_1 - time_0)
### Prepare variables
N_LIB = 500 ### total number of libraries simulated - for the computation of population sizes
str(dat)
dat$REP = as.factor(dat$REP)
dat$GRAD = as.factor(dat$GRAD) ### categorical variable!
dat$LOG10_RMSD = log10(dat$RMSD+1)
dat$TPR_OVER_TPR_PLUS_FDR = dat$TRUE_POSITIVE_RATE / (dat$TRUE_POSITIVE_RATE + dat$FALSE_DISCOVERY_RATE) ### ranges from 0: which means no QTL detected, to 1: which means all QTL were detected and no false positives
dat$NPOP_TRAIN = unlist(lapply(strsplit(as.character(dat$POP_TRAIN), ";"), FUN=length))
dat$NPOP_TEST = unlist(lapply(strsplit(as.character(dat$POP_TEST), ";"), FUN=length))
dat$POPSIZE_TRAIN = floor(N_LIB / dat$NPOP_TRAIN)
dat$POPSIZE_TEST = floor(N_LIB / dat$NPOP_TEST)
dat = dat[dat$TRUE_POSITIVE_RATE <= 1.0, ] ### for some reason I'm gettingTPR greater than 1.00 even after removing ERRORED datasets!
dat$MODEL = paste0(gsub("-", "", as.character(dat$MODEL_ITERATION)), "_", as.character(dat$MODEL_MODEL))
### clean-up and save
dat = droplevels(dat)
str(dat)
saveRDS(dat, "GPASim_data.rds")
# dat = readRDS("GPASim_data.rds")
######################
### LOAD LIBRARIES ###
######################
library(RColorBrewer)
library(lme4)
# remotes::install_github("jeffersonfparil/violinplotter")
library(violinplotter)
library(doParallel)
###########################################
### METRIC AND FACTORS NAMES AND LABELS ###
###########################################
### the response variables we measured (GPAS performance metrics):
METRICS = c("CORRELATION", "LOG10_RMSD", "TRUE_POSITIVE_RATE", "FALSE_DISCOVERY_RATE", "TPR_OVER_TPR_PLUS_FDR")
METRIC_LABELS = c("Correlation", "log10(RMSD+1)", "True Positive Rate", "False Discovery Rate", "TPR/(TPR + FDR)")
### the explanatory variables that we cannot control but interesting nonetheless (population genetics factors):
GENFACTORS = c("QTL", "MGR", "BGS", "GRAD")
GENFACTOR_LABELS = c("Number of QTL", "Migration Rate", "Background Selection Intensity", "QTL Emergence Gradient")
### the explanatory variables we can control (GPAS modelling factors):
###### where:
######        - each sampling methodology (i.e. across, within x individual, pool)
######          will be modelled separately because we have already established the proportion
######          of resources that needs to be allocated to each of these sampling methodology to
######          maximize GPAS performance during ABC Resouce Allocation Optimization steps.
######          And finally, the complexity makes the computation too computationally expensive!
MODFACTORS = c("MODEL", "MODEL_COVAR", "NPOP_TRAIN", "NPOP_TEST")
MODFACTOR_LABELS = c("Model", "Model Covariate", "Number of Training Populations", "Number of Validation Populations")
### numeric variables (will be used as both nester and nested) while categorical variables will be used as nesters alone!
NUMERICS = c(GENFACTORS[c(1,2,3)], MODFACTORS[c(3,4)]); NUMERICS_LABELS = c(GENFACTOR_LABELS[c(1,2,3)], MODFACTOR_LABELS[c(3,4)])
CATEGORICALS = c(GENFACTORS[c(4)], MODFACTORS[c(1,2)]); CATEGORICALS_LABELS = c(GENFACTOR_LABELS[c(4)], MODFACTOR_LABELS[c(1,2)])
### list the different combinations of subsetting parameters and random effect variable (nester) to use for mixed models
PARAMS_MIXED = list()
counter = 1
for (i in 1:nlevels(dat$GROUPING)){
  for (j in 1:nlevels(dat$ALGORITHM)){
    for (k in 1:length(METRICS)){
      for (l in 1:length(c(NUMERICS, CATEGORICALS))){
        PARAMS_MIXED[[counter]] = c(i,j,k,l)
        counter = counter + 1
      }
    }
  }
}
### list the different combinations of subsetting parameters to use for building the fixed models
PARAMS_FIXED = list()
counter = 1
for (i in 1:nlevels(dat$GROUPING)){
  for (j in 1:nlevels(dat$ALGORITHM)){
    for (k in 1:length(METRICS)){
      PARAMS_FIXED[[counter]] = c(i,j,k)
      counter = counter + 1
    }
  }
}
### define the function for data subsetting and mixed model building
MIXED_MODELING_FUNCTION = function(par){
  # par = c(2,1,1,1)
  i = par[1] # across vs. within populations sampling
  j = par[2] # inidividual vs. pool genotyping
  k = par[3] # GPAS metrics
  l = par[4] # nesting (random) variable
  across_or_within = as.character(levels(dat$GROUPING))[i]
  indi_vs_pool = as.character(levels(dat$ALGORITHM))[j]
  subdat = droplevels(subset(subset(dat, ALGORITHM==indi_vs_pool), GROUPING==across_or_within))
  metric = METRICS[k]
  DATA = eval(parse(text=paste0("data.frame(", paste(paste0("subdat$", c(metric, NUMERICS, CATEGORICALS)), collapse=","), ")")))
  colnames(DATA) = c(metric, NUMERICS, CATEGORICALS)
  DATA = droplevels(DATA[!is.na(DATA[,1]) & !is.infinite(DATA[,1]) & complete.cases(DATA), ])
  DATA[,2:(length(NUMERICS)+1)] = scale(DATA[,2:(length(NUMERICS)+1)], scale=TRUE, center=TRUE)
  random_var = c(NUMERICS, CATEGORICALS)[l]
  fixed_vars = c(NUMERICS, CATEGORICALS)[c(NUMERICS, CATEGORICALS) != random_var]
  random_interactions = paste0("(0+", NUMERICS[NUMERICS != random_var], "|", random_var, ")")
  # random_interactions = c(paste0("(0+", NUMERICS[NUMERICS != random_var], "|", random_var, ")"),
  #                         paste0("(1|", CATEGORICALS[CATEGORICALS != random_var], ":", random_var, ")"))
  random_interactions[1] = gsub("0\\+", "", random_interactions[1])
  if(across_or_within == "WITHIN"){
    random_var = random_var[(random_var != "NPOP_TRAIN") & (random_var != "NPOP_TEST")]
    fixed_vars = fixed_vars[(fixed_vars != "NPOP_TRAIN") & (fixed_vars != "NPOP_TEST")]
    random_interactions = random_interactions[!grepl("NPOP_TRAIN", random_interactions) & !grepl("NPOP_TEST", random_interactions)]
  }
  stopifnot(length(random_var) == 1) ### stops if we get the non-varying NPOP_TRAIN and NPOP_TEST variables
  MODEL_MIXED = tryCatch(lmer(as.formula(paste0(metric, "~", paste(c(fixed_vars, random_interactions), collapse="+"))), data=DATA),
                warning=function(e){
                  lmer(as.formula(paste0(metric, "~", paste(c(fixed_vars, random_interactions), collapse="+"))), data=DATA, control=lmerControl(optimizer="Nelder_Mead"))
                })
  if (isSingular(MODEL_MIXED)){
    print("Removing **TEST** population-related variables, and re-modeling!")
    MODEL_MIXED = tryCatch(lmer(as.formula(paste0(metric, "~", paste(c(fixed_vars, random_interactions[!grepl("TEST", random_interactions)]), collapse="+"))), data=DATA),
                  warning=function(e){
                    lmer(as.formula(paste0(metric, "~", paste(c(fixed_vars, random_interactions[!grepl("TEST", random_interactions)]), collapse="+"))), data=DATA, control=lmerControl(optimizer="Nelder_Mead"))
                  })
  }
  saveRDS(MODEL_MIXED, paste0("MODEL_MIXED-", across_or_within, "-", indi_vs_pool, "-", metric, "-", random_var, ".rds"))
  print(paste0("END: ", date()))
  print("######################################################################")
  return(0)
}
### define the function for data subsetting and fixed model building
FIXED_MODELING_FUNCTION = function(par){
  # par = c(2,1,1)
  i = par[1] # across vs. within populations sampling
  j = par[2] # inidividual vs. pool genotyping
  k = par[3] # GPAS metrics
  across_or_within = as.character(levels(dat$GROUPING))[i]
  indi_vs_pool = as.character(levels(dat$ALGORITHM))[j]
  subdat = droplevels(subset(subset(dat, ALGORITHM==indi_vs_pool), GROUPING==across_or_within))
  metric = METRICS[k]
  DATA = eval(parse(text=paste0("data.frame(", paste(paste0("subdat$", c(metric, NUMERICS, CATEGORICALS)), collapse=","), ")")))
  colnames(DATA) = c(metric, NUMERICS, CATEGORICALS)
  DATA = droplevels(DATA[!is.na(DATA[,1]) & !is.infinite(DATA[,1]) & complete.cases(DATA), ])
  # standardize numeric explanatory variables
  means_sds_list = list(MEAN=apply(DATA[,2:(length(NUMERICS)+1)], MARGIN=2, FUN=mean, na.rm=TRUE), SD=apply(DATA[,2:(length(NUMERICS)+1)], MARGIN=2, FUN=sd, na.rm=TRUE))
  saveRDS(means_sds_list, paste0("MEAN_SD-", across_or_within, "-", indi_vs_pool, "-", metric, ".rds"))
  DATA[,2:(length(NUMERICS)+1)] = scale(DATA[,2:(length(NUMERICS)+1)], scale=TRUE, center=TRUE)
  EXPLANATORIES = c(NUMERICS, CATEGORICALS)
  if (across_or_within == "WITHIN"){
    DATA = droplevels(DATA[, !grepl("TRAIN", colnames(DATA)) & !grepl("TEST", colnames(DATA))])
    EXPLANATORIES = EXPLANATORIES[!grepl("TRAIN", EXPLANATORIES) & !grepl("TEST", EXPLANATORIES)]
  }
  MODEL_FIXED_0_NULL = lm(as.formula(paste0(metric, "~ 1")), data=DATA)
  MODEL_FIXED_1_FULLADD = lm(as.formula(paste0(metric, "~ .")), data=DATA)
  MODEL_FIXED_2_FULLFAC = lm(as.formula(paste0(metric, "~", paste(EXPLANATORIES, collapse="*"))), data=DATA)
  saveRDS(MODEL_FIXED_0_NULL, paste0("MODEL_FIXED-", across_or_within, "-", indi_vs_pool, "-", metric, "_NULL.rds"))
  saveRDS(MODEL_FIXED_1_FULLADD, paste0("MODEL_FIXED-", across_or_within, "-", indi_vs_pool, "-", metric, "-FULLADD.rds"))
  saveRDS(MODEL_FIXED_2_FULLFAC, paste0("MODEL_FIXED-", across_or_within, "-", indi_vs_pool, "-", metric, "-FULLFAC.rds"))
  print(paste0("END: ", date()))
  print("######################################################################")
  return(0)
}
### perform mixed modelling
nCores = 9
registerDoParallel(nCores)
out_parallel = foreach(x=PARAMS_MIXED) %dopar% MIXED_MODELING_FUNCTION(x)
gc()
### perform fixed model building
nCores = 4
registerDoParallel(nCores)
foreach(x=PARAMS_FIXED[1:10]) %dopar% FIXED_MODELING_FUNCTION(x)
gc()
nCores = 1
registerDoParallel(nCores)
foreach(x=PARAMS_FIXED[11:20]) %dopar% FIXED_MODELING_FUNCTION(x)
gc()
### plotting
# define color schems for each explanatory variable
QTL_COLOURS = RColorBrewer::brewer.pal(length(unique(dat$QTL)), "GnBu")
MGR_COLOURS = RColorBrewer::brewer.pal(3, "PuRd")[2:3]
BGS_COLOURS = RColorBrewer::brewer.pal(length(unique(dat$BGS)), "YlGn")
NPOP_TRAIN_COLOURS = c(RColorBrewer::brewer.pal(9, "PuRd")[1:2], RColorBrewer::brewer.pal(9, "PuRd"))
NPOP_TEST_COLOURS = c(RColorBrewer::brewer.pal(9, "YlGnBu")[1:2], RColorBrewer::brewer.pal(9, "YlGnBu"))
GRAD_COLOURS = RColorBrewer::brewer.pal(length(unique(dat$GRAD)), "RdYlBu")
MODEL_COLOURS = RColorBrewer::brewer.pal(length(unique(dat$MODEL)), "Set3")
MODEL_COVAR_COLOURS = RColorBrewer::brewer.pal(length(unique(dat$MODEL_COVAR)), "Pastel1")
SOLID_CATEGORICAL_COLOURS = unique(c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Dark2")))
# clean-up
rm(dat)
gc()
# across the 4 sampling schemes
SAMPLING_SCHEMES = c("ACROSS-INDI", "ACROSS-POOL", "WITHIN-INDI", "WITHIN-POOL")
for (i in 1:length(SAMPLING_SCHEMES)){
  # i = 1
  sampling_scheme = SAMPLING_SCHEMES[i]
  # fixed null, full additive and full factorial across metrics
    for (j in 1:length(METRICS)){
    # j = 1
    metric = METRICS[j]
    MEAN_SD = readRDS(paste0("MEAN_SD-", sampling_scheme, "-", metric, ".rds"))
    mod0 = readRDS(paste0("MODEL_FIXED-", sampling_scheme, "-", metric, "-NULL.rds"))
    if ((sampling_scheme == "ACROSS-INDI") | (sampling_scheme == "ACROSS-POOL")){
      mod1 = readRDS(paste0("MODEL_FIXED-", sampling_scheme, "-", metric, "-FULLADD.rds"))
      mod2 = readRDS(paste0("MODEL_FIXED-", sampling_scheme, "-", metric, "-FULLFAC.rds"))
    }
    # mixed models per metric across nester
    for (k in 1:length(c(NUMERICS, CATEGORICALS))){
      # k = 7 ## k = 1
      nester = c(NUMERICS, CATEGORICALS)[k]
      mod3 = tryCatch(readRDS(paste0("MODEL_MIXED-", sampling_scheme, "-", metric, "-", nester, ".rds")),
                error=function(e){
                  print(paste0(sampling_scheme, "-", metric, "-", nester))
                  write.csv(c("ERROR: mixed model does not exist!"), file=paste0("ERROR-", sampling_scheme, "-", metric, "-", nester, ".csv"))
                })
      if(is.null(mod3)){
        next
      }
      # model assessment
      summary_mod = summary(mod3)
      V_fixed = vcov(mod3)
      V_random = VarCorr(mod3)
      if ((sampling_scheme == "ACROSS-INDI") | (sampling_scheme == "ACROSS-POOL")){
        AIC_TEST = AIC(mod0, mod1, mod2, mod3)
        BIC_TEST = BIC(mod0, mod1, mod2, mod3)
      } else {
        AIC_TEST = AIC(mod0, mod3)
        BIC_TEST = BIC(mod0, mod3)
      }
      # extracting data from the model
      mod_dat = mod3@frame
      col_mod_dat = colnames(mod_dat)
      FACTORS_LIST = col_mod_dat[!grepl(metric, col_mod_dat) & !grepl(nester, col_mod_dat)]
      FORMULA = as.formula(paste0(metric, "~", paste(FACTORS_LIST, collapse="+")))
      NESTER_LEVELS = levels(eval(parse(text=paste0("mod_dat$", nester))))
      if (length(NESTER_LEVELS)==0){
        # for numeric nesters
        if (length(eval(parse(text=paste0("unique(round(mod_dat$", nester, ", 4))")))) == length(eval(parse(text=paste0("unique(mod_dat$", nester, ")"))))) {
          # test if rounding the numerics to 4 decimal places sill gets us the same unique numbers as if they were not rounded-off
          eval(parse(text=paste0("mod_dat$", nester, " = as.factor(as.character(round(mod_dat$", nester, ", 4)))")))
        }
        NESTER_LEVELS = levels(eval(parse(text=paste0("mod_dat$", nester))))
      }
      # plotting mixed model results
      ## fixed variables:
      mod3_fixef = fixef(mod3)
      mod3_fixef_names = names(mod3_fixef)
      x = c(-1,0,1)
      y = lapply(mod3_fixef[2:length(mod3_fixef)], FUN=function(b){mod3_fixef[1] + (b*x)})
      mod3_fixef_pvals = 2 * (1 - pnorm(abs(data.frame(coef(summary_mod))$t.value)))
      mod3_fixef_sig = lapply(mod3_fixef_pvals, FUN=function(x){if(x<=0.05){"*"} else if(x<=0.01){"**"} else if(x<=0.001){"***"} else {"ns"}})
      svg(paste0("LINEPLOTS-", sampling_scheme, "-", metric, "-NESTER_", nester, "-FIXED_EFFECTS.svg"), width=10, height=7)
        layout(mat=matrix(c(1,1,2), nrow=1))
        plot(x=c(-1, 1), y=c(min(unlist(y)),max(unlist(y))), type="n", xlab="Standardized Numeric Factor", ylab=metric, main=paste0(sampling_scheme, "-", metric, "\nFixed Effects with nester: ", nester))
        for (xi in 2:length(mod3_fixef_names)){
          # xi = 2
          fixef_xi = mod3_fixef_names[xi]
          lines(x=x, y=eval(parse(text=paste0("y$", fixef_xi))), lty=1, lwd=2, col=SOLID_CATEGORICAL_COLOURS[xi])
        }
        mod3_fixef_order = order(matrix(unlist(y), ncol=3, byrow=TRUE)[,3], decreasing=TRUE)
        plot(x=1,y=1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        legend_test = paste0(mod3_fixef_names[2:length(mod3_fixef_names)][mod3_fixef_order], " (", mod3_fixef_sig[2:length(mod3_fixef_names)][mod3_fixef_order], ")")
        legend("center", legend=legend_test, col=SOLID_CATEGORICAL_COLOURS[2:length(mod3_fixef_names)][mod3_fixef_order], lty=1, lwd=2, bty="n")
        if (length(summary(mod3)$optinfo$conv$lme4$messages)>0){legend("bottomright", legend=summary(mod3)$optinfo$conv$lme4$messages)}
      dev.off()
      ## random variables:
      mod3_ranef = eval(parse(text=paste0("ranef(mod3)$", nester)))
      mod3_ranef_names = colnames(mod3_ranef)[2:ncol(mod3_ranef)]
      x = c(-1,0,1)
      Y1 = mod3_ranef[,1] + (mod3_ranef[,2:ncol(mod3_ranef)]*x[1])
      Y2 = mod3_ranef[,1] + (mod3_ranef[,2:ncol(mod3_ranef)]*x[2])
      Y3 = mod3_ranef[,1] + (mod3_ranef[,2:ncol(mod3_ranef)]*x[3])
      RANSD = data.frame(V_random)
      mod3_ranef_sig = RANSD$sdcor[(RANSD$var1!="(Intercept)") & !is.na(RANSD$var1)] / RANSD$sdcor[nrow(RANSD)]
      svg(paste0("LINEPLOTS-", sampling_scheme, "-", metric, "-NESTER_", nester, "-RANDOM_EFFECTS.svg"), width=10, height=7)
        len_nester = length(NESTER_LEVELS) + 1 ### added 1 for the legend
        plot_dims = c(round(sqrt(len_nester)), (floor(sqrt(len_nester)) + ceiling(sqrt(len_nester) %% floor(sqrt(len_nester)))))
        par(mfrow=plot_dims)
        for (i_nester_level in 1:length(NESTER_LEVELS)){
          # i_nester_level = 1
          nester_level = NESTER_LEVELS[i_nester_level]
          # convert the numeric nester level to it's non-standard normal form
          if (sum(grepl(nester, NUMERICS))>0){
            if ((nester == "MGR") | (nester == "BGS")){
              nester_level_name = round((as.numeric(nester_level) * MEAN_SD$SD[nester]) + MEAN_SD$MEAN[nester], 4)
            } else {
              nester_level_name = round((as.numeric(nester_level) * MEAN_SD$SD[nester]) + MEAN_SD$MEAN[nester])
            }
          } else {
            nester_level_name = nester_level
          }
          plot(x=c(-1,1), y=c(min(unlist(Y1,Y2,Y3)), max(unlist(Y1,Y2,Y3))), type="n", ylab=metric, xlab="Standardized Numeric Factor", main=paste0(sampling_scheme, "-", metric, "\nFixed Effects with nester: ", nester, ": ", nester_level_name))
          for (j_factor in 1:ncol(Y1)){
            # j_factor = 1
            # y = mod3_ranef[i_nester_level, 1] + (mod3_ranef[i_nester_level, j_factor] * x)
            y = c(Y1[i_nester_level, j_factor], Y2[i_nester_level, j_factor], Y3[i_nester_level, j_factor])
            lines(x=x, y=y, lty=1, lwd=2, col=SOLID_CATEGORICAL_COLOURS[j_factor])
          }
        }
        plot(x=1,y=1,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
        legend_text = paste0(mod3_ranef_names, " (sd/sderr = ", round(mod3_ranef_sig,2), ")")
        legend("center",legend=legend_text, col=SOLID_CATEGORICAL_COLOURS[1:length(mod3_ranef_names)], lwd=2, bty="n")
        if (length(summary(mod3)$optinfo$conv$lme4$messages)>0){legend("bottomright", legend=summary(mod3)$optinfo$conv$lme4$messages)}
      dev.off()
      # FIGURE BY NESTER LEVEL: iterate across nester leves
      for (l in 1:length(NESTER_LEVELS)){
        # l = 1
        nester_level = NESTER_LEVELS[l]
        sub_mod_dat = droplevels(eval(parse(text=paste0("subset(mod_dat,",  nester, "=='", nester_level, "')"))))
        XCATEGOR = FACTORS_LIST %in% CATEGORICALS
        REGRESSX = !XCATEGOR
        LOGX = FACTORS_LIST == "QTL"
        LOGX_BASE = rep(1, times=length(LOGX)); LOGX_BASE[LOGX] = 5
        VIOLIN_COLOURS = eval(parse(text=paste0("list(", paste(paste0(FACTORS_LIST, "_COLOURS"), collapse=","), ")")))
        # convert the numeric nester level to it's non-standard normal form
        if (sum(grepl(nester, NUMERICS))>0){
          if ((nester == "MGR") | (nester == "BGS")){
            nester_level_name = round((as.numeric(nester_level) * MEAN_SD$SD[nester]) + MEAN_SD$MEAN[nester], 4)
          } else {
            nester_level_name = round((as.numeric(nester_level) * MEAN_SD$SD[nester]) + MEAN_SD$MEAN[nester])
          }
        } else {
          nester_level_name = nester_level
        }
        svg(paste0("VIOLINPLOTS-", sampling_scheme, "-", metric, "-NESTER_", nester, "-LEVEL_", nester_level_name, "-ACROSS_FACTORS.svg"), width=10, height=15)
        tryCatch(
          violinplotter(formula=FORMULA, data=sub_mod_dat, LOGX=LOGX, LOGX_BASE=LOGX_BASE, XCATEGOR=XCATEGOR, REGRESSX=REGRESSX, VIOLIN_COLOURS=VIOLIN_COLOURS, TITLE=rep(paste0(nester, ":\n", nester_level_name), times=length(FACTORS_LIST))),
          error=function(e){
            new_fomula = as.formula(paste0(metric, "~", paste(FACTORS_LIST[FACTORS_LIST != "NPOP_TEST"], collapse="+")))
            tryCatch(violinplotter(formula=new_fomula, data=sub_mod_dat, LOGX=LOGX, LOGX_BASE=LOGX_BASE, XCATEGOR=XCATEGOR, REGRESSX=REGRESSX, VIOLIN_COLOURS=VIOLIN_COLOURS, TITLE=rep(paste0(nester, ":\n", nester_level_name), times=length(FACTORS_LIST))),
              error=function(e){
                new_fomula = as.formula(paste0(metric, "~", paste(FACTORS_LIST[(FACTORS_LIST != "NPOP_TEST") & (FACTORS_LIST != "NPOP_TRAIN")], collapse="+")))
                tryCatch(violinplotter(formula=new_fomula, data=sub_mod_dat, LOGX=LOGX, LOGX_BASE=LOGX_BASE, XCATEGOR=XCATEGOR, REGRESSX=REGRESSX, VIOLIN_COLOURS=VIOLIN_COLOURS, TITLE=rep(paste0(nester, ":\n", nester_level_name), times=length(FACTORS_LIST))),
                  error=function(e){
                    err_n = nrow(sub_mod_dat)
                    err_dat = cbind(sub_mod_dat, rep(c("A", "B"), each=ceiling(err_n/2))[1:err_n])
                    eval(parse(text=paste0("err_dat$", metric, " = rnorm(err_n)")))
                    colnames(err_dat)[ncol(err_dat)] = "ERROR"
                    err_formula = as.formula(paste0(metric, "~ ERROR"))
                    violinplotter(formula=err_formula, data=err_dat, LOGX=LOGX, LOGX_BASE=LOGX_BASE, XCATEGOR=XCATEGOR, REGRESSX=REGRESSX, VIOLIN_COLOURS=VIOLIN_COLOURS, TITLE=rep(paste0(nester, ":\n", nester_level_name), times=length(FACTORS_LIST)))
                    legend("center", legend="ERROR\n(dummy variable)")
                })
            })
        })
        dev.off()
      }
      # FIGURE BY FACTOR: iterate across factors and composite plot for each nester level
      for (m in 1:length(FACTORS_LIST)){
        # m = 1
        x_var = FACTORS_LIST[m]
        formula = as.formula(paste0(metric, "~", x_var))
        VIOLIN_COLOURS = eval(parse(text=paste0(x_var, "_COLOURS")))
        ### setting violinplotter parameters
        if (sum(grepl(x_var, CATEGORICALS))>0){
          # categoricals
          XCATEGOR=TRUE
          REGRESSX=FALSE
          LOGX=FALSE
        } else {
          # numerics
          XCATEGOR=FALSE
          REGRESSX=TRUE
          if(x_var=="QTL"){
            LOGX=TRUE
            LOGX_BASE=5
          } else {
            LOGX=FALSE
          }
        }
        len_nester = length(NESTER_LEVELS)
        plot_dims = c(round(sqrt(len_nester)), (floor(sqrt(len_nester)) + ceiling(sqrt(len_nester) %% floor(sqrt(len_nester)))))
        if(x_var=="MODEL"){
          svg(paste0("VIOLINPLOTS-", sampling_scheme, "-", metric, "-NESTER_", nester, "-FACTOR_", x_var, "-ACROSS_LEVELS.svg"), width=10, height=15)
        } else {
          svg(paste0("VIOLINPLOTS-", sampling_scheme, "-", metric, "-NESTER_", nester, "-FACTOR_", x_var, "-ACROSS_LEVELS.svg"), width=10, height=10)
        }
        par(mfrow=plot_dims)
        for (n in 1:length(NESTER_LEVELS)){
          # n = 1
          nester_level = NESTER_LEVELS[n]
          # convert the numeric nester level to it's non-standard normal form
          if (sum(grepl(nester, NUMERICS))>0){
            if ((nester == "MGR") | (nester == "BGS")){
              nester_level_name = round((as.numeric(nester_level) * MEAN_SD$SD[nester]) + MEAN_SD$MEAN[nester], 4)
            } else {
              nester_level_name = round((as.numeric(nester_level) * MEAN_SD$SD[nester]) + MEAN_SD$MEAN[nester])
            }
          } else {
            nester_level_name = nester_level
          }
          sub_mod_dat = droplevels(eval(parse(text=paste0("subset(mod_dat,",  nester, "=='", nester_level, "')"))))
          tryCatch(
            violinplotter(formula=formula, data=sub_mod_dat, TITLE=paste0(nester, ":\n", nester_level_name), XLAB=x_var, LOGX=LOGX, LOGX_BASE=LOGX_BASE, XCATEGOR=XCATEGOR, REGRESSX=REGRESSX, VIOLIN_COLOURS=VIOLIN_COLOURS),
            error=function(e){
              err_n = nrow(sub_mod_dat)
              err_dat = cbind(sub_mod_dat, rep(c("A", "B"), each=ceiling(err_n/2))[1:err_n])
              eval(parse(text=paste0("err_dat$", metric, " = rnorm(err_n)")))
              colnames(err_dat)[ncol(err_dat)] = "ERROR"
              err_formula = as.formula(paste0(metric, "~ ERROR"))
              violinplotter(formula=err_formula, data=err_dat, TITLE=paste0(nester, ":\n", nester_level_name), XLAB=x_var, LOGX=LOGX, LOGX_BASE=LOGX_BASE, XCATEGOR=XCATEGOR, REGRESSX=REGRESSX, VIOLIN_COLOURS=VIOLIN_COLOURS)
              legend("center", legend="ERROR\n(dummy variable)")
            })
        }
        dev.off()
      }
    }
    rm(mod_dat); gc()
  }
}
