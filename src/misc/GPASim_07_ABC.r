#################################
### ABC Optimization Analysis ###
#################################

#####################################################################################################
### Load data
args = commandArgs(trailingOnly=TRUE)
# args = c("ABC_MERGED_RESOURCE_OPTIM.csv") ### merged results of each ABC optimization for GPAS performance per simulated landscape
dat = read.csv(args[1])
str(dat)
# dat = dat[dat$TPR_OVER_TPR_PLUS_FDR <= 1.0, ] ### for some reason I'm gettingTPR greater than 1.00 even after removing ERRORED datasets!
dat = dat[dat$TPR_OVER_TPR_PLUS_FDR < 1.0, ] ### for some reason I'm gettingTPR greater than 1.00 even after removing ERRORED datasets!

#####################################################################################################
### Load libraries
# library(vioplot)
# tryCatch({
#   library(agricolae); ### hard to install - will use git repo and load specific functions required for HSD.test()
#   "R-installed library loaded!"
#   },
#   error=function(e){LIB_PATH =  .libPaths()[1];
#                   source(paste0(LIB_PATH, "/agricolae/R/order.group.R"));
#                   source(paste0(LIB_PATH, "/agricolae/R/orderPvalue.R"));
#                   source(paste0(LIB_PATH, "/agricolae/R/lastC.R"));
#                   source(paste0(LIB_PATH, "/agricolae/R/tapply.stat.R"));
#                   source(paste0(LIB_PATH, "/agricolae/R/HSD.test.R"));
#                   "Git repo agricolae scripts loaded!"
#                 }
# )
library("lme4")
tryCatch(
  {
    remotes::install_github("jeffersonfparil/violinplotter")
    library(violinplotter)
  }, error=function(e){
    install.packages("remotes")
    remotes::install_github("jeffersonfparil/violinplotter")
    library(violinplotter)
  }
)

#####################################################################################################
### Labels and COLOURS5
METRIC = c("LOG10_RMSD", "TPR_OVER_TPR_PLUS_FDR")
# METRIC_LABEL = c(expression(log[10] ~ (RMSD+1)), "True Positive Rate")
METRIC_LABEL = c("log10(RMSD+1)", "True Positive Rate / \n(True Positive Rate + False Discovery Rate)")
ABC_CLASS = c("NULL", "GPAS", "GP", "GWAS")
ALLOCATION_CLASS = c("ACROSS_POOL", "ACROSS_INDI", "WITHIN_POOL", "WITHIN_INDI")
ALLOCATION_LABEL = c("Across Populations &\nPool Genotyping", "Across Populations &\nIndividual Genotyping", "Within Populations &\nPool Genotyping", "Within Populations &\nIndividual Genotyping")
COLOURS1 = RColorBrewer::brewer.pal(4, "Pastel1")
COLOURS2 = RColorBrewer::brewer.pal(4, "Set1")
COLOURS3 = rep(COLOURS1, times=length(unique(dat$QTL)))
COLOURS4 = rep(COLOURS2, times=length(unique(dat$QTL)))
COLOURS5 = c("#68AED4", "#C5D9EE", "#ABDC8C", "#F5FAB8")
COLOURS6 = c("#B3E2CD", "#FBB4AE", "#FDCDAC", "#CCCCCC")

# # ###################
# ### Violin plots
# violinplot_func = function(FOR_VIOPLOT, response, explanatory, colors, label, xlab='', ylab='', las=2, legend=TRUE){
#   explanatory_vars = unlist(strsplit(explanatory, "\\+"))
#   nclasses = 1
#   for (i in explanatory_vars){
#     nclasses = nclasses * eval(parse(text=paste0("nlevels(as.factor(FOR_VIOPLOT$", i, "))")))
#   }
#   means = eval(parse(text=paste0("aggregate(", response ,"~", explanatory, ", FUN=mean, data=FOR_VIOPLOT)")))
#   stds = eval(parse(text=paste0("aggregate(", response ,"~", explanatory, ", FUN=sd, data=FOR_VIOPLOT)")))
#   plot_violin = eval(parse(text=paste0("vioplot(jitter(", response, ")~", explanatory, ", data=FOR_VIOPLOT, col=c(", paste(rep(colors, times=nclasses), collapse=","), "), border=NA, drawRect=FALSE, xlab='", xlab, "', ylab='", ylab, "', main='", label, "', cex.main=1.0, cex.axis=1.0, las=", las, ")")))
#   CI = qnorm(((0.95)/2)+0.50) * (stds[,ncol(stds)] / sqrt(nrow(FOR_VIOPLOT)-1)) ### 95% confidence interval
#   points(1:nclasses, means[,ncol(means)], pch=20)
#   arrows(x0=1:nclasses, y0=means[,ncol(means)]+CI, y1=means[,ncol(means)]-CI, angle=90, code=3, length=0.1)
#   grid()
#   if (legend==TRUE){
#     legend("bottomright", legend="95% CI", bty="o")
#   }
# }

#####################################################################################################
### Prelimiary plots
### Partition
# svg("ABC_optimization_allocation_vioplots.svg", width=10, height=7)
# par(mfrow=c(2,2), mar=c(7,5,4,1))
for (i in 1:length(ALLOCATION_CLASS)){
  # i = 1
  response = ALLOCATION_CLASS[i]
  explanatory = "ABC_CLASS*QTL"
  label = ALLOCATION_LABEL[i]
  svg(paste0("ABC_optimization_allocation_vioplots_", response,".svg"), width=16, height=4)
  # jpeg(paste0("ABC_optimization_allocation_vioplots_", response,".jpg"), width=2000, height=700, quality=100)
  # violinplot_func(FOR_VIOPLOT=dat, response=response, explanatory=explanatory, colors="COLOURS1", ylab="Resources Allocated", label=paste0(label, " X ", explanatory))
  violinplotter(formula=as.formula(paste0(response, " ~ ", explanatory)), dat, LOGX=c(FALSE, TRUE), LOGX_BASE=c(1, 5), XCATEGOR=c(TRUE, FALSE), REGRESSX=c(FALSE,TRUE),
                VIOLIN_COLOURS=COLOURS1, TITLE=paste0(label, "\nX\n", c("ABC_CLASS", "QTL", "ABC_CLASS:QTL Interactions")))
  dev.off()
}
svg("ABC_optimization_allocation_histograms.svg", width=10, height=7)
par(cex=2, mfrow=c(2,2), mar=c(5,5,4,1))
for (i in 1:length(ALLOCATION_CLASS)){
  response = ALLOCATION_CLASS[i]
  label = ALLOCATION_LABEL[i]
  hist(eval(parse(text=paste0("dat$", response))), col=COLOURS1[i], border=COLOURS2[i], xlim=c(0, 1), xlab="Fraction of Resources Allocated", ylab="Frequency", main=label)
}
dev.off()
### Metrics
# svg("ABC_optimization_metrics_vioplots.svg", width=7, height=9)
# par(mfrow=c(2,1), mar=c(7,5,4,1))
for (i in 1:length(METRIC)){
  # i = 1
  response = METRIC[i]
  explanatory = "ABC_CLASS*QTL"
  label = METRIC_LABEL[i]
  svg(paste0("ABC_optimization_metrics_vioplots_", response,".svg"), width=16, height=4)
  # violinplot_func(FOR_VIOPLOT=dat, response=response, explanatory=explanatory, colors="COLOURS3", ylab="Resources Allocated", label=paste0(label, " X ", explanatory))
  violinplotter(formula=as.formula(paste0(response, " ~ ", explanatory)), dat, LOGX=c(FALSE, TRUE), LOGX_BASE=c(1, 5), XCATEGOR=c(TRUE, FALSE), REGRESSX=c(FALSE,TRUE),
                VIOLIN_COLOURS=COLOURS3, TITLE=paste0(label, "\nX\n", c("ABC_CLASS", "QTL", "ABC_CLASS:QTL Interactions")))
  dev.off()
}

ABC_modelling_and_plotting = function(dat, id_name=""){
  #####################################################################################################
  ### Linear modelling for mean comparisons between optimization methods: NULL vs GPAS vs GWAS vs GP
  ### converting QTL, MGR, BGS and GRAD into categorical factors for simplicity of mean comparisons
  dat$QTL = as.factor(dat$QTL)
  dat$MGR = as.factor(dat$MGR)
  dat$BGS = as.factor(dat$BGS)
  dat$GRAD = as.factor(dat$GRAD)
  ### iterate across allocation fraction and expected metric values
  rm(HSD_OUT)
  ALLOCATION_AND_METRICS = colnames(dat)[8:ncol(dat)]
  for (i in ALLOCATION_AND_METRICS) {
    # i = "TPR_OVER_TPR_PLUS_FDR"
    mod_null = tryCatch(eval(parse(text=paste0("lm(", i, "~ 1, data=dat)"))), error=function(e){NULL})
    mod_saturated = tryCatch(eval(parse(text=paste0("lm(", i, "~ ABC_CLASS*QTL*MGR*BGS*GRAD, data=dat)"))), error=function(e){NULL})
    mod_saturated_randABC = tryCatch(eval(parse(text=paste0("lmer(", i, "~ (1|ABC_CLASS)*QTL*MGR*BGS*GRAD, data=dat)"))), error=function(e){NULL})
    mod_saturated_randQTL = tryCatch(eval(parse(text=paste0("lmer(", i, "~ ABC_CLASS*(1|QTL)*MGR*BGS*GRAD, data=dat)"))), error=function(e){NULL})
    mod_saturated_randMGR = tryCatch(eval(parse(text=paste0("lmer(", i, "~ ABC_CLASS*QTL*(1|MGR)*BGS*GRAD, data=dat)"))), error=function(e){NULL})
    mod_saturated_randBGS = tryCatch(eval(parse(text=paste0("lmer(", i, "~ ABC_CLASS*QTL*MGR*(1|BGS)*GRAD, data=dat)"))), error=function(e){NULL})
    mod_saturated_randGRAD = tryCatch(eval(parse(text=paste0("lmer(", i, "~ ABC_CLASS*QTL*MGR*BGS*(1|GRAD), data=dat)"))), error=function(e){NULL})
    mod_fullAdditive = tryCatch(eval(parse(text=paste0("lm(", i, "~ ABC_CLASS+QTL+MGR+BGS+GRAD, data=dat)"))), error=function(e){NULL})
    mod_fullAdditive_randABC = tryCatch(eval(parse(text=paste0("lmer(", i, "~ (1|ABC_CLASS)+QTL+MGR+BGS+GRAD, data=dat)"))), error=function(e){NULL})
    mod_fullAdditive_randQTL = tryCatch(eval(parse(text=paste0("lmer(", i, "~ ABC_CLASS+(1|QTL)+MGR+BGS+GRAD, data=dat)"))), error=function(e){NULL})
    mod_fullAdditive_randMGR = tryCatch(eval(parse(text=paste0("lmer(", i, "~ ABC_CLASS+QTL+(1|MGR)+BGS+GRAD, data=dat)"))), error=function(e){NULL})
    mod_fullAdditive_randBGS = tryCatch(eval(parse(text=paste0("lmer(", i, "~ ABC_CLASS+QTL+MGR+(1|BGS)+GRAD, data=dat)"))), error=function(e){NULL})
    mod_fullAdditive_randGRAD = tryCatch(eval(parse(text=paste0("lmer(", i, "~ ABC_CLASS+QTL+MGR+BGS+(1|GRAD), data=dat)"))), error=function(e){NULL})
    mod_abcQTL = tryCatch(eval(parse(text=paste0("lm(", i, "~ ABC_CLASS*QTL, data=dat)"))), error=function(e){NULL})
    mod_abcQTL_randABC = tryCatch(eval(parse(text=paste0("lmer(", i, "~ (1|ABC_CLASS)*QTL, data=dat)"))), error=function(e){NULL})
    mod_abcQTL_randQTL = tryCatch(eval(parse(text=paste0("lmer(", i, "~ ABC_CLASS*(1|QTL), data=dat)"))), error=function(e){NULL})
    mod_abc = tryCatch(eval(parse(text=paste0("lm(", i, "~ ABC_CLASS, data=dat)"))), error=function(e){NULL})
    mod_abc_randABC = tryCatch(eval(parse(text=paste0("lmer(", i, "~ (1| ABC_CLASS), data=dat)"))), error=function(e){NULL})
    ### model comparisons (using BIC only)
    rm(models_list)
    models_list = c()
    models_list_initial = ls()[grep("mod_", ls())]
    for (mod in models_list_initial){
      if (!is.null(eval(parse(text=mod)))){
        models_list = c(models_list, mod)
      }
    }
    CHISQ_TEST = tryCatch({out=eval(parse(text=paste0("anova(", paste(models_list, collapse=","), ", test='Chisq')"))); out}, error=function(e){"ERROR!"})
    AIC_TEST = eval(parse(text=paste0("AIC(", paste0(models_list, collapse=","), ")")))
    BIC_TEST = eval(parse(text=paste0("BIC(", paste0(models_list, collapse=","), ")")))
    BIC_TEST = BIC_TEST[(rownames(BIC_TEST) != "mod_null"), ]
    mod_selected = rownames(BIC_TEST[(BIC_TEST$BIC == min(BIC_TEST$BIC)), ])[1]
    ### just peeking into the best among the models tested
    print(mod_selected)
    # print(summary(eval(parse(text=mod_selected))) )
    print(anova(eval(parse(text=mod_selected))) )
    # hsd = eval(parse(text=paste0("HSD.test(", mod_selected, ", trt='ABC_CLASS')$groups")))
    formula = as.formula(eval(parse(text=mod_selected))$call[[2]])
    hsd = mean_comparison_HSD(formula=formula, data=dat, explanatory_variable_name="ABC_CLASS")
    # hsd = hsd[order(rownames(hsd)), ]
    hsd = hsd[c(4,3,1,2), ]
    # HSD = data.frame(ALLOCATION_AND_METRICS=rep(i, times=nrow(hsd)), ABC_CLASS=rownames(hsd), VALUE=hsd[,1], GROUP=hsd[,2])
    HSD = data.frame(ALLOCATION_AND_METRICS=rep(i, times=nrow(hsd)), ABC_CLASS=hsd$LEVELS, VALUE=hsd$MEANS, GROUP=hsd$GROUPING)
    if (exists("HSD_OUT")) {
      HSD_OUT = rbind(HSD_OUT, HSD)
    } else {
      HSD_OUT = HSD
    }
    rm(mod_selected)
  }
  write.table(HSD_OUT, paste0("ABC_optimization_HSD", id_name, ".csv"), row.names=FALSE, sep=",")
  ### plot bar and violin plots stacked for allocation and side-by-side for the metrics
  svg(paste0("ABC_optimization_comparisons_barplots", id_name, ".svg"), width=11, height=9)
  rm(DATA_COMPARISON)
  for (abc_class in ABC_CLASS) {
    if (exists("DATA_COMPARISON")) {
      DATA_COMPARISON = droplevels(rbind(DATA_COMPARISON, subset(HSD_OUT, ABC_CLASS==abc_class)))
    } else {
      DATA_COMPARISON = droplevels(subset(HSD_OUT, ABC_CLASS==abc_class))
    }
  }
  par(mfrow=c(1, 2), mar=c(5, 5, 1, 1))
  ### stacked resource allocation barplot
  ### violin plots
  plot_layout_matrix = matrix(c(1,2,1,3), ncol=2, byrow=TRUE)
  layout(plot_layout_matrix)
  ### stacked resource allocation barplot
  FOR_STACKED_BARPLOT = matrix(DATA_COMPARISON$VALUE, byrow=FALSE, ncol=length(ABC_CLASS))[1:(length(ALLOCATION_AND_METRICS)-2), ]
  COLNAMES = unique(DATA_COMPARISON$ABC_CLASS)
  ROWNAMES = unique(DATA_COMPARISON$ALLOCATION_AND_METRICS)[1:4]
  colnames(FOR_STACKED_BARPLOT) = COLNAMES
  rownames(FOR_STACKED_BARPLOT) = ROWNAMES
  FOR_STACKED_BARPLOT = FOR_STACKED_BARPLOT[, order(colnames(FOR_STACKED_BARPLOT))]
  ylim=c(0, 1.5)
  par(mar=c(5, 5, 1, 1))
  barplot(FOR_STACKED_BARPLOT[order(rownames(FOR_STACKED_BARPLOT), decreasing=TRUE), ], col=COLOURS5, bord=COLOURS5, xlab="Optimization Method", ylab="Fraction Allocated",
          ylim=ylim, legend.text=TRUE, args.legend=list(y=ylim[2], bty="n", cex=1.0), cex.axis=1.0, cex.names=1.0, cex.lab=1.0)
  ### cloud plots
  FOR_VIOPLOT = data.frame(ABC_CLASS=dat$ABC_CLASS, LOG10_RMSD=dat$LOG10_RMSD, TPR_OVER_TPR_PLUS_FDR=dat$TPR_OVER_TPR_PLUS_FDR)
  FOR_VIOPLOT$ABC_CLASS = factor(FOR_VIOPLOT$ABC_CLASS, levels=ABC_CLASS) ### arrange by ABC_CLASS. i..e. null, gpas, gp and finally gwas
  MEANS = list(aggregate(LOG10_RMSD ~ ABC_CLASS, FUN=mean, data=FOR_VIOPLOT), aggregate(TPR_OVER_TPR_PLUS_FDR ~ ABC_CLASS, FUN=mean, data=FOR_VIOPLOT))
  STDS = list(aggregate(LOG10_RMSD ~ ABC_CLASS, FUN=sd, data=FOR_VIOPLOT), aggregate(TPR_OVER_TPR_PLUS_FDR ~ ABC_CLASS, FUN=sd, data=FOR_VIOPLOT))
  par(mar=c(3, 3, 3, 1))
  for (i in 1:length(METRIC)){
    # i = 1
    metric = METRIC[i]
    metric_label = METRIC_LABEL[i]
    explanatory_variable_name="ABC_CLASS"
    # violinplot_func(FOR_VIOPLOT=FOR_VIOPLOT, response=metric, explanatory="ABC_CLASS", colors="COLOURS6", label=metric_label, las=1, HSD=TRUE)
    plot_violin_1x(dat=FOR_VIOPLOT, response_variable_name=metric, explanatory_variable_name=explanatory_variable_name, COLOURS=COLOURS6, title=metric_label)
    mean_comparison_HSD(formula=as.formula(paste0(metric, " ~ ", explanatory_variable_name)), data=FOR_VIOPLOT, explanatory_variable_name=explanatory_variable_name, PLOT=TRUE)
  }
  dev.off()
  return(HSD)
}

### all data
ABC_modelling_and_plotting(dat=dat, id_name="")
### extremes: QTL=5, QTL=100, MGR=0.0001, GRAD=0
ABC_modelling_and_plotting(dat=droplevels(dat[dat$QTL==5, ]), id_name="-EXTREMES-QTL5")
ABC_modelling_and_plotting(dat=droplevels(dat[dat$QTL==100, ]), id_name="-EXTREMES-QTL100")
ABC_modelling_and_plotting(dat=droplevels(dat[dat$MGR==0.0001, ]), id_name="-EXTREMES-MGR0.0001")
ABC_modelling_and_plotting(dat=droplevels(dat[dat$GRAD==0, ]), id_name="-EXTREMES-GRAD0")
