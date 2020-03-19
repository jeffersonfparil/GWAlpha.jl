########################
### Linear modelling ###
########################

### (2) Sensible mixed modes ...
### ... where we are interested in the effects of:
###   - the number of training populations and
###   - the different GPAS models
###   - as well as the number of simulated QTL for within population data sets
### while controlling for the random (assumed random) effects of:
###   - the number of  simulated QTL (applicable only for across population datasets),
###   - migration rate,
###   - background selection intensity and
###   - resistance allele emergence gradient.
### We are trying to peer more closely with eyes for sensible models (create an Rscript ("sensible_mixed_modelling.R") with the following contents:)
args = commandArgs(trailingOnly=TRUE) ### csv file
# args = c("GPASim_data.rds")
### load libraries
library(doParallel)
library(lme4)
library(violinplotter)
### load data
dat = readRDS(args[1])
### lumping MODEL and MODEL_COVARIATE together into one categorical variable
dat$MODEL = as.factor(paste0(as.character(dat$MODEL), "_", as.character(dat$MODEL_COVARIATE)))
SAMPLING_SCHEMES = c("ACROSS-INDI", "ACROSS-POOL", "WITHIN-INDI", "WITHIN-POOL")
METRICS = c("CORRELATION", "LOG10_RMSD", "TRUE_POSITIVE_RATE", "FALSE_DISCOVERY_RATE", "TPR_OVER_TPR_PLUS_FDR")
# COLORS_MODELS = pals::tableau20()
colors_models = rep(c(RColorBrewer::brewer.pal(8,"Set2"), RColorBrewer::brewer.pal(8,"Set3")), times=2)[1:nlevels(dat$MODEL)]
COLORS_MODELS = data.frame(IDX=1:length(colors_models), MODEL=levels(dat$MODEL), COLORS=colors_models)
### modelling function
SENSIBLE_MIXED_MODELS = function(idx){
  # idx = c(3,1)
  i = idx[1]
  j = idx[2]
  sampling_scheme = SAMPLING_SCHEMES[i]
  metric = METRICS[j]
  idx =  (paste(dat$GROUPING, dat$ALGORITHM, sep="-") == sampling_scheme) & !is.na(eval(parse(text=paste0("dat$", metric)))) & !is.infinite(eval(parse(text=paste0("dat$", metric))))
  DATA = droplevels(dat[idx, ])
  if ((metric=="LOG10_RMSD") | (metric=="FALSE_DISCOVERY_RATE")){
    INTERCEPT="1"
  } else {
    INTERCEPT="0"
  }
  if (grepl("ACROSS", sampling_scheme)){
    ### We are interested on the effects of CONTROLLABLE FACTORS: (1) the number of training populations, and (2) GPAS model on GPAS performance
    ### and we want to control for the effects of UNCOTROLLABLE FACTORS: (1) the number of QTL, (2) migration rate, (3) QTL emergence gradient, and (4) background selection intensity
    ### and also we want to control for the effect of the number of validation populations which we assume to be a random effect
    MODEL_0_NULL = eval(parse(text=paste0("lm(", metric, " ~ 1, data=DATA)")))
    MODEL_1_FULL = eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + NPOP_TRAIN + MODEL + (1|QTL) + (1|MGR) + (1|GRAD) + (1|BGS) + (1|NPOP_TEST), data=DATA, control=lmerControl(optimizer='Nelder_Mead'))")))
    MODEL_2_QTL = eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + NPOP_TRAIN + MODEL + (1|QTL), data=DATA, control=lmerControl(optimizer='Nelder_Mead'))")))
    MODEL_3_MGR = eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + NPOP_TRAIN + MODEL + (1|MGR), data=DATA, control=lmerControl(optimizer='Nelder_Mead'))")))
    MODEL_4_GRAD = eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + NPOP_TRAIN + MODEL + (1|GRAD), data=DATA, control=lmerControl(optimizer='Nelder_Mead'))")))
    MODEL_5_BGS = eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + NPOP_TRAIN + MODEL + (1|BGS), data=DATA, control=lmerControl(optimizer='Nelder_Mead'))")))
    MODEL_6_NPOPTEST = eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + NPOP_TRAIN + MODEL + (1|NPOP_TEST), data=DATA, control=lmerControl(optimizer='Nelder_Mead'))")))
    ### Are the trends in GPAS performance as number of training populations increase differ across different GPAS models?
    ### To do this properly we need to scale the number of training populations
    scaled_DATA = DATA
    # scaled_DATA$NPOP_TRAIN = scale(DATA$NPOP_TRAIN, scale=TRUE, center=TRUE) + (mean(DATA$NPOP_TRAIN)/sd(DATA$NPOP_TRAIN)) ### add a non-centrality parameter that moves the minimum back to zero so we can fix the y-intercept at zero (equivalent to center=FALSE and scale=TRUE)
    scaled_DATA$NPOP_TRAIN = scale(DATA$NPOP_TRAIN, scale=TRUE, center=FALSE)
    MODEL = tryCatch(
            eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + (NPOP_TRAIN|MODEL) + (1|QTL) + (1|MGR) + (1|GRAD) + (1|BGS) + (1|NPOP_TEST), data=scaled_DATA, control=lmerControl(optimizer='Nelder_Mead'))"))),
            warning=function(e){
            eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + (NPOP_TRAIN|MODEL) + (1|QTL) + (1|MGR) + (1|GRAD) + (1|BGS) + (1|NPOP_TEST), data=scaled_DATA, REML=FALSE, control=lmerControl(optimizer='Nelder_Mead', check.conv.grad=.makeCC('warning', tol = 0.5, relTol = NULL)))")))
          })
    ### how good is this new model with the previous ones
    AIC_test = AIC(MODEL_0_NULL, MODEL_1_FULL, MODEL_2_QTL, MODEL_3_MGR, MODEL_4_GRAD, MODEL_5_BGS, MODEL_6_NPOPTEST, MODEL)
    BIC_test = BIC(MODEL_0_NULL, MODEL_1_FULL, MODEL_2_QTL, MODEL_3_MGR, MODEL_4_GRAD, MODEL_5_BGS, MODEL_6_NPOPTEST, MODEL)
  } else {
    ### reduce the size of WITHIN POPULATION DATASETS by random sampling 200,000 data rows
    set.seed(55)
    DATA = DATA[sample(1:nrow(DATA), size=200000), ]
    ### conert the number of QTL to log(QTL,5)
    DATA$QTL = log(DATA$QTL, 5)
    ### build the models
    MODEL_0_NULL = eval(parse(text=paste0("lm(", metric, " ~ 1, data=DATA)")))
    MODEL_1_FULL = eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + MODEL + QTL + (1|MGR) + (1|GRAD) + (1|BGS), data=DATA, control=lmerControl(optimizer='Nelder_Mead'))")))
    # MODEL_2_QTL = eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + MODEL + (1|QTL), data=DATA, control=lmerControl(optimizer='Nelder_Mead'))")))
    MODEL_3_MGR = eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + MODEL  + QTL + (1|MGR), data=DATA, control=lmerControl(optimizer='Nelder_Mead'))")))
    MODEL_4_GRAD = eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + MODEL  + QTL + (1|GRAD), data=DATA, control=lmerControl(optimizer='Nelder_Mead'))")))
    MODEL_5_BGS = eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + MODEL  + QTL + (1|BGS), data=DATA, control=lmerControl(optimizer='Nelder_Mead'))")))
    # MODEL_6_NPOPTEST = eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + MODEL + (1|NPOP_TEST), data=DATA, control=lmerControl(optimizer='Nelder_Mead'))")))
    ### Are the trends in GPAS performance as number of QTL increase differ across different GPAS models? ### Just wanted the models for across and within populations to be as symmetric as possible ;-P
    ### To do this properly we need to scale the number of QTL
    scaled_DATA = DATA
    scaled_DATA$QTL = scale(DATA$QTL, scale=TRUE, center=FALSE)
    MODEL = tryCatch(
            eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + (QTL|MODEL) + (1|MGR) + (1|GRAD) + (1|BGS), data=scaled_DATA, control=lmerControl(optimizer='Nelder_Mead'))"))),
            warning=function(e){
            eval(parse(text=paste0("lmer(", metric, " ~ ", INTERCEPT, " + (QTL|MODEL) + (1|MGR) + (1|GRAD) + (1|BGS), data=scaled_DATA, REML=FALSE, control=lmerControl(optimizer='Nelder_Mead', check.conv.grad=.makeCC('warning', tol = 0.5, relTol = NULL)))")))
          })
    ### how good is this new model with the previous ones
    AIC_test = AIC(MODEL_0_NULL, MODEL_1_FULL, MODEL_3_MGR, MODEL_4_GRAD, MODEL_5_BGS, MODEL)
    BIC_test = BIC(MODEL_0_NULL, MODEL_1_FULL, MODEL_3_MGR, MODEL_4_GRAD, MODEL_5_BGS, MODEL)
  }
  print("###########################################")
  print(paste0("SENSIBLE_MODEL-", sampling_scheme, "-", metric))
  print(AIC_test)
  print(BIC_test)
  print(paste0(rownames(BIC_test)[BIC_test$BIC == min(BIC_test$BIC)], " is the best model!"))
  saveRDS(BIC_test, file=paste0("BIC_TEST-", sampling_scheme, "-", metric, ".rds"))
  saveRDS(MODEL_1_FULL, file=paste0("SENSIBLE_MODEL-", sampling_scheme, "-", metric, ".rds"))
  saveRDS(MODEL, file=paste0("SENSIBLE_MODEL_NPOP_TRAIN_SLOPE-", sampling_scheme, "-", metric, ".rds"))
  print("###########################################")
  ### plot the effects of the number of training populations (or QTL) and GPAS model (fixed effects from SENSIBLE_MODEL-*)
  if (grepl("ACROSS", sampling_scheme)){
    nested = "NPOP_TRAIN"
  } else {
    nested = "QTL"
  }
  svg(paste0("SENSIBLE_MODEL-SCATTERBAR-", sampling_scheme, "-", metric, ".svg"), width=20, height=10)
    old_par = par(no.readonly=TRUE)
    layout(matrix(c(1,1,2,2,2), byrow=TRUE, nrow=1))
    par(mar=c(8.5, 4, 8.5, 3), cex=1)
    # plot mean points with sd bars
    agg_means = merge(eval(parse(text=paste0("aggregate(", metric, " ~ MODEL + ", nested, ", data=DATA, FUN=mean)"))), COLORS_MODELS, by="MODEL")
    agg_stdev = merge(eval(parse(text=paste0("aggregate(", metric, " ~ MODEL + ", nested, ", data=DATA, FUN=sd)"))), COLORS_MODELS, by="MODEL")
    Y0 = eval(parse(text=paste0("agg_means$", metric))) - eval(parse(text=paste0("agg_stdev$", metric)))
    Y1 = eval(parse(text=paste0("agg_means$", metric))) + eval(parse(text=paste0("agg_stdev$", metric)))
    if (INTERCEPT=="1"){
      y_intercept = fixef(MODEL_1_FULL)[names(fixef(MODEL_1_FULL))=="(Intercept)"]
    } else {
      y_intercept = 0.00
    }
    newx = sort(unique(eval(parse(text=paste0("DATA$", nested)))))                        ### for regression line plotting and setting ylim
    newy = y_intercept + (newx * fixef(MODEL_1_FULL)[names(fixef(MODEL_1_FULL))==nested]) ### for regression line plotting and setting ylim
    if (nested == "NPOP_TRAIN"){
      plot(x=eval(parse(text=paste0("agg_means$", nested))), y=eval(parse(text=paste0("agg_means$", metric))), ylim=c(min(c(Y0,Y1, newy)), max(c(Y0,Y1, newy))), col=paste0(agg_means$COLORS, "80"), pch=19, cex=2,
          xlab="Number of Training Populations", ylab=metric)
    } else {
      plot(x=eval(parse(text=paste0("agg_means$", nested))), y=eval(parse(text=paste0("agg_means$", metric))), ylim=c(min(c(Y0,Y1)), max(c(Y0,Y1))), col=paste0(agg_means$COLORS, "80"), pch=19, cex=2,
          xlab="log5(Number of QTL)", ylab=metric)
    }
    mtext(side=3, line=3, adj=0, cex=1, paste0(nested, " (Fixed Effect)"))
    mtext(side=3, line=2, adj=0, cex=0.75, as.character(MODEL_1_FULL@call)[2])
    arrows(x0=eval(parse(text=paste0("agg_means$", nested))), y0=Y0, y1=Y1, angle=90, code=3, lwd=1, length=0.1, col=paste0(agg_means$COLORS, "80"))
    grid()
    # plot regression line
    lines(x=newx, y=newy, lty=2, lwd=2, col="grey55")
    coef = data.frame(summary(MODEL_1_FULL)$coefficients)
    rownames(coef) = gsub("^MODEL", "", rownames(coef))
    tval = coef$t.val[rownames(coef)==nested]
    if(tval < 0){
      pval = 2*(pt(tval, df=1, lower.tail=TRUE))
    } else {
      pval = 2*(pt(tval, df=1, lower.tail=FALSE))
    }
    legend("bottomright", legend=c(paste0("y-intercept = ", round(y_intercept,4)),
                                  paste0("slope = ", round(fixef(MODEL_1_FULL)[names(fixef(MODEL_1_FULL))==nested], 4)),
                                  paste0("p-value = ", round(pval, 4))
                                  ), bty="n")
    # metric ~ MODEL: barplot with sd bars
    fixed_effects = coef[(rownames(coef) != nested) & (rownames(coef) != "(Intercept)"), ]
    fixed_effects$p.value = unlist(lapply(fixed_effects$t.value, FUN=function(x){if (x < 0){out = 2*pt(x, df=1, lower.tail=TRUE)} else {out = 2*pt(x, df=1, lower.tail=FALSE)};return(out)}))
    fixed_effects$MODEL = rownames(fixed_effects)
    fixed_effects = merge(fixed_effects, COLORS_MODELS, by="MODEL")
    fixed_effects = fixed_effects[order(fixed_effects$Estimate, decreasing=TRUE), ]
    old_par = par(no.readonly=TRUE)
    par(mar=c(20, 5, 5, 2), cex=1)
    bp = barplot(height=fixed_effects$Estimate,
                 border=FALSE,
                 horiz=FALSE,
                 names.arg=fixed_effects$MODEL,
                #  ylim=c(min(c(0, fixed_effects$Estimate - (2*fixed_effects$Std..Error))), max(fixed_effects$Estimate + (2*fixed_effects$Std..Error))),
                 ylim=c(min(c(0, fixed_effects$Estimate-fixed_effects$Std..Error, fixed_effects$Estimate+fixed_effects$Std..Error)), max(c(0, fixed_effects$Estimate-fixed_effects$Std..Error, fixed_effects$Estimate+fixed_effects$Std..Error))),
                 col=as.character(fixed_effects$COLORS),
                 las=2,
                 ylab=metric)
    mtext(side=3, line=3, adj=0, cex=1, "MODEL (Fixed Effect)")
    mtext(side=3, line=2, adj=0, cex=0.75, as.character(MODEL_1_FULL@call)[2])
    if (INTERCEPT == "1"){
      ### specify the contrasting categorical variable level which lost it's identity when we have a non-zero intercept
      mtext(side=3, line=1, adj=0, cex=0.75, paste0("y-intercept = ", levels(DATA$MODEL)[1]))
    }
    arrows(x0=bp, y0=fixed_effects$Estimate-fixed_effects$Std..Error, y1=fixed_effects$Estimate+fixed_effects$Std..Error, angle=90, code=3, lwd=2, length=0.1, col="grey55")
    par(old_par)
  dev.off()
  ### plot the effects of the number of training populations (or QTL) nested within GPAS models (random slopes from SENSIBLE_MODEL_NPOP_TRAIN_SLOPE-*)
  svg(paste0("SENSIBLE_MODEL-SLOPES_PER_MODEL-", sampling_scheme, "-", metric, ".svg"), width=20, height=10)
    U = as.data.frame(ranef(MODEL)$MODEL)
    # U$`(Intercept)` = eval(parse(text=paste0("U$", nested))) * sd(eval(parse(text=paste0("DATA$", nested))))
    # eval(parse(text=paste0("U$", nested, " = U$", nested, " * sd(DATA$", nested, ")")))
    U$`(Intercept)` = eval(parse(text=paste0("U$", nested)))

    U$MODEL = as.factor(rownames(U))
    U = merge(U, COLORS_MODELS, by="MODEL")
    # newx = sort(unique(c(0, eval(parse(text=paste0("scaled_DATA$", nested))))))
    newx = sort(unique(eval(parse(text=paste0("scaled_DATA$", nested)))))
    U$MIN_POINT = U$`(Intercept)` + (eval(parse(text=paste0("U$", nested))) * min(newx))
    U$MAX_POINT = U$`(Intercept)` + (eval(parse(text=paste0("U$", nested))) * max(newx))
    U = U[order(U$MAX_POINT, decreasing=TRUE), ]
    xlim = c(min(eval(parse(text=paste0("DATA$", nested)))), max(eval(parse(text=paste0("DATA$", nested)))))
    ylim = c(min(c(U$MIN_POINT, U$MAX_POINT)), max(c(U$MIN_POINT, U$MAX_POINT)))
    old_par = par(no.readonly=TRUE)
    layout(matrix(c(1,1,1,2), byrow=TRUE, nrow=1))
    par(cex=1.5, mar=c(5,5,5,1))
    if (nested == "NPOP_TRAIN"){
      plot(x=xlim, y=ylim, type="n", xlab="Number of Training Populations", ylab=metric)
    } else {
      plot(x=xlim, y=ylim, type="n", xlab="log5(Number of QTL)", ylab=metric)
    }
    mtext(side=3, line=3, adj=0, cex=1, paste0(nested, " slopes per MODEL (Random Slopes)"))
    mtext(side=3, line=1, adj=0, cex=0.75, as.character(MODEL@call)[2])
    for (k in 1:nrow(U)){
      # k = 1
      lines(x=xlim, y=c(U$MIN_POINT[k], U$MAX_POINT[k]), lty=1, lwd=3, col=as.character(U$COLORS)[k])
    }
    grid()
    par(mar=c(1,1,1,1))
    plot(x=0,y=0,type="n", xaxt="n", yaxt="n",bty="n",xlab="",ylab="")
    legend("center", legend=U$MODEL, col=as.character(U$COLORS), lty=1, lwd=3, bty="n", cex=0.5)
    par(old_par)
  dev.off()
  ### violin plots of the trend in GPAS performance as number of training populations (or QTL) increases nested within GPAS models
  svg(paste0("SENSIBLE_MODEL-VIOLINPLOTS_PER_MODEL-", sampling_scheme, "-", metric, ".svg"), width=20, height=20)
    MODEL_NAMES_ORDERED = data.frame(MODEL_NAMES=fixed_effects$MODEL, ORDER=c(1:nrow(fixed_effects)))
    MODEL_NAMES_UNORDERED = data.frame(MODEL_NAMES=as.character(unique(DATA$MODEL)), FACTOR_MODEL_NAMES=unique(DATA$MODEL))
    MODEL_NAMES_MERGED = merge(MODEL_NAMES_ORDERED, MODEL_NAMES_UNORDERED, by="MODEL_NAMES")
    MODEL_NAMES_MERGED = MODEL_NAMES_MERGED[order(MODEL_NAMES_MERGED$ORDER, decreasing=FALSE),]
    model_names = MODEL_NAMES_MERGED$FACTOR_MODEL_NAMES
    len_nester = length(model_names)
    plot_dims = c(round(sqrt(len_nester)), (floor(sqrt(len_nester)) + ceiling(sqrt(len_nester) %% floor(sqrt(len_nester)))))
    par(mfrow=plot_dims)
    for (k in 1:length(model_names)){
      # k = 1
      model = model_names[k]
      sub_DATA = droplevels(DATA[DATA$MODEL == model, ])
      sub_DATA = merge(sub_DATA, COLORS_MODELS, by="MODEL")
      f = as.formula(paste0(metric, " ~ ", nested))
      tryCatch(violinplotter(formula=f, data=sub_DATA, XCATEGOR=FALSE, HSDX=TRUE, LOGX=FALSE, REGRESSX=TRUE, TITLE=paste0("MODEL:\n", model_names[k]), VIOLIN_COLOURS=as.character(sub_DATA$COLORS)),
        error=function(e){
               eval(parse(text=paste0("sub_DATA$", metric, " = rnorm(nrow(sub_DATA))")))
               eval(parse(text=paste0("sub_DATA$", nested, " = rep(c('A', 'B'), each=ceiling(nrow(sub_DATA)/2))[1:nrow(sub_DATA)]")))
               violinplotter(formula=f, data=sub_DATA, XCATEGOR=FALSE, HSDX=TRUE, LOGX=FALSE, REGRESSX=TRUE, TITLE=paste0("MODEL:\n", model_names[k]), VIOLIN_COLOURS=sub_DATA$COLORS)
               legend("center", legend=c("ERROR!", "Dummy variables!"))
        })
    }
  dev.off()
  rm(DATA)
  gc()
  return(0)
}
### parallel execution
params_idx = list(); counter=1
for (i in 1:length(SAMPLING_SCHEMES)){
  for (j in 1:length(METRICS)){
    params_idx[[counter]] = c(i,j)
    counter = counter + 1
  }
}
nCores = 32
registerDoParallel(nCores)
out_parallel = foreach(x=params_idx) %dopar% SENSIBLE_MIXED_MODELS(x)
gc()
### execution:
# time nohup Rscript sensible_mixed_modelling.R &
