module GWAlpha

### NOTE:
### Also known as the PoolGPAS_module.jl in https://gitlab.com/jeffersonfparil/genomic_prediction

#####################
### load packages ###
#####################
using DelimitedFiles
using Distributions
using LinearAlgebra
using Optim
using GLMNet
using RCall
using ProgressMeter
using DataFrames
using CSV
using ColorBrewer
# JULIA_SCRIPT_HOME = @__DIR__
# # # test:
# # JULIA_SCRIPT_HOME = "/data/Lolium/Softwares/genomic_prediction/src"
# push!(LOAD_PATH, JULIA_SCRIPT_HOME)
# using sync_parsing_module
# using filter_sync_module
# using LMM_module
# using GP_module
# using pval_heuristic_module
include("sync_parsing_module.jl")
include("filter_sync_module.jl")
include("poolFST_module.jl")
include("LMM_module.jl")
include("GP_module.jl")
include("pval_heuristic_module.jl")
### for parallel computation of GWAlpha_ML_parallel_module.jl
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
include("GWAlpha_ML_parallel_module.jl")
@everywhere using GWAlpha.GWAlpha_ML_parallel_module

############################
### function definitions ###
############################
function neg_log_likelihood_cdfbeta(beta::Array{Float64,1}, data_A::Array{Float64,1}, data_B::Array{Float64,1})
	-sum(
		 log.(10,
		 	 (Distributions.cdf.(Distributions.Beta(beta[1], beta[2]), data_A)) .-
		 	 (Distributions.cdf.(Distributions.Beta(beta[1], beta[2]), append!(zeros(1), data_A[1:(length(data_A)-1)])))
		     )
		 )     -sum(
		 log.(10,
		 	 (Distributions.cdf.(Distributions.Beta(beta[3], beta[4]), data_B)) .-
		 	 (Distributions.cdf.(Distributions.Beta(beta[3], beta[4]), append!(zeros(1), data_B[1:(length(data_B)-1)])))
		     )
		 )
end

function cost_function_RR(lambda::Array{Float64,1}, y::Array{Float64,1}, X::Array{Float64,2})
	# beta_hat = LMM_module.inverse((X' * X) + (lamda * I)) * X' * y ### approximates very closely with the less costly:
	beta_hat = X' * LMM_module.inverse((X * X') + (lambda[1] * I)) * y
	y_dev = y - (X * beta_hat)
	COST = y_dev' * y_dev
	return(COST)
end

function GWAlpha_ML(filename_sync::String, filename_phen_py::String, MAF::Float64)
	### load the sync and phenotype files
	sync = DelimitedFiles.readdlm(filename_sync, '\t')
	phen = DelimitedFiles.readdlm(filename_phen_py)

	### gather phenotype specifications
	NPOOLS = length(split(phen[5], ['=', ',', '[', ']', ';'])) - 3 #less the first leading and trailing elemenets
	if length(split(phen[1], ['=', '\"'])) < 3
		global NAME = split(phen[1], ['=', '\''])[3]
	else
		global NAME = split(phen[1], ['=', '\"'])[3]
	end
	SD = parse.(Float64, split(phen[2], ['=',';'])[2])
	MIN = parse.(Float64, split(phen[3], ['=',';'])[2])
	MAX = parse.(Float64, split(phen[4], ['=',';'])[2])
	PERC = parse.(Float64, split(phen[5], ['=', ',', '[', ']', ';'])[3:(NPOOLS+1)])
	QUAN = parse.(Float64, split(phen[6], ['=', ',', '[', ']', ';'])[3:(NPOOLS+1)])
	BINS = append!([x for x in PERC], 1) - append!(zeros(1), PERC)

	### gather genotype (allele frequency) specificications
	NSNP = size(sync)[1]
	n_pools_sync = size(sync)[2] - 3
	if NPOOLS != n_pools_sync
		println("The number of pools with phenotype data does not match the number of pools with allele frequency data!")
		println("Please check you input files :-)")
		println("Remove leading and intervening whitespaces in the phenotype file.")
		exit()
	else
		n_pools_sync = nothing #clear out contents of this redundant n_pools variable
	end

	### iterate across SNPs
	ALPHA_OUT = []
	ALLELE_FREQ = []
	ALLELE_ID = []
	LOCUS_ID = []
	LOCUS_W_EFF = [] # (1 for with effects, -999 for no effects of filtered out alleles) multiplier of ALPHA_OUT and LOD to exclude alleles with zero effects --> prevents the distortion of the distribution of alpha as affected by the zero effect alleles
	COUNTS = zeros(Int64, NPOOLS, 6) #nrow=n_pools and ncol=A,T,C,G,N,DEL
	progress_bar = ProgressMeter.Progress(NSNP, dt=1, desc="GWAlpha_ML Progress: ",  barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow) #progress bar
	for snp in 1:NSNP
		# println(snp)
		#parse allele counts from the sync file
		for i in 1:NPOOLS
			COUNTS[i,:] = [parse(Int64, x) for x in split.(sync[snp, 4:(NPOOLS+3)], [':'])[i]]
		end
		#convert to frequencies per pool
		FREQS = COUNTS ./ ( sum(COUNTS, dims=2) .+ 1e-10 ) #added 1e-10 to the denominator to avoid NAs in pools with no allele counts (zero depth; which should actually have been filtered out after mpileup using awk)
		allele_freqs = sum(FREQS .* BINS, dims=1)
		#iterate across alleles while filtering by MAF
		if (sum(COUNTS) != 0.0)
			if (minimum(allele_freqs[allele_freqs .!= 0.0]) >= MAF) & (maximum(allele_freqs) < (1.0 - MAF)) #locus filtering by mean MAF
				for allele in 1:6
					if (allele_freqs[allele] > 0.0) & (maximum(FREQS[:,allele]) < 0.999999)  #allele filtering remove alleles with no counts and that the number of pools with allele frequency close to one should not occur even once!
					# if (sum(FREQS[:,allele] .== 0.0) < NPOOLS) & (sum(FREQS[:,allele] .> 0.999999) < 1) #filter-out alleles with at least 1 pool fixed for that allele because it causes a failure in the optimization
						freqA = FREQS[:, allele]
						pA = sum(freqA .* BINS)
						pB = 1 - pA
						BINA = (freqA .* BINS) ./ pA
						BINB = ( (1 .- freqA) .* BINS ) ./ (1-pA)
						percA = cumsum(BINA)
						percB = cumsum(BINB)

						### optimize (minimize) -log-likelihood of these major allele frequencies modelled as a beta distribution
						# using Nelder-Mead optimization or Box minimisation (try-catch if one or the other fails with preference to Nelder-Mead)
						lower_limits = [1e-20, 1e-20, 1e-20, 1e-20]
						upper_limits = [1.0, 1.0, 1.0, 1.0]
						initial_values = [0.1, 0.1, 0.1, 0.1]
						BETA = try
							Optim.optimize(beta->neg_log_likelihood_cdfbeta(beta, percA, percB), initial_values, NelderMead())
						catch
							try
								Optim.optimize(beta->neg_log_likelihood_cdfbeta(beta, percA, percB), lower_limits, upper_limits, initial_values)
							catch ### lower limits of 1e-20 to 1e-6 causes beta dist parameter values to shrink to zero somehow - so we'r'e setting lower limits to 1e-5 instead
								lower_limits = [1e-5, 1e-5, 1e-5, 1e-5]
								Optim.optimize(beta->neg_log_likelihood_cdfbeta(beta, percA, percB), lower_limits, upper_limits, initial_values)
							end
						end
						MU_A = MIN + ((MAX-MIN)*BETA.minimizer[1]/(BETA.minimizer[1]+BETA.minimizer[2]))
						MU_B = MIN + ((MAX-MIN)*BETA.minimizer[3]/(BETA.minimizer[3]+BETA.minimizer[4]))

						### compute alpha
						W_PENAL = 2*sqrt(pA*pB)
						ALPHA = W_PENAL*(MU_A - MU_B) / SD
						append!(ALPHA_OUT, ALPHA)
						append!(ALLELE_ID, allele)
						append!(LOCUS_ID, snp)
						append!(LOCUS_W_EFF, 1)
						append!(ALLELE_FREQ, pA)
					else
						### for explicit zero effects for null (low to none) frequency alleles
						# append!(ALPHA_OUT, 0.0)
						# append!(ALLELE_ID, allele)
						# append!(LOCUS_ID, snp)
						# append!(LOCUS_W_EFF, 0)
						# append!(ALLELE_FREQ, allele_freqs[allele])
					end
				end
			end
		end
		# println("$snp")
		ProgressMeter.update!(progress_bar, snp)
	end

	### estimate heuristic p-values
	# alpha_mean = mean(ALPHA_OUT[LOCUS_W_EFF .== 1])
	# alpha_sd = std(ALPHA_OUT[LOCUS_W_EFF .== 1])
	# P_VALUES = ones(length(LOCUS_ID))
	# P_VALUES[LOCUS_W_EFF .== 1] .= [pval_Normal(x, alpha_mean, alpha_sd) for x in ALPHA_OUT[LOCUS_W_EFF .== 1]]
	# LOD = -Distributions.log.(10, P_VALUES)
	P_VALUES, LOD = pval_heuristic_module.estimate_PVAL_and_LOD(convert(Array{Float64,1}, ALPHA_OUT))

	### output
	for i in 1:length(ALLELE_ID) #convert int allele ID into corresponting A, T, C, G, N, DEL
		if ALLELE_ID[i] == 1; ALLELE_ID[i] = "A"
		elseif ALLELE_ID[i] == 2; ALLELE_ID[i] = "T"
		elseif ALLELE_ID[i] == 3; ALLELE_ID[i] = "C"
		elseif ALLELE_ID[i] == 4; ALLELE_ID[i] = "G"
		elseif ALLELE_ID[i] == 5; ALLELE_ID[i] = "N"
		elseif ALLELE_ID[i] == 6; ALLELE_ID[i] = "DEL"
		end
	end
	OUT = DataFrames.DataFrame(LOCUS_ID=LOCUS_ID, CHROM=sync[LOCUS_ID,1], POS=sync[LOCUS_ID,2], ALLELE=ALLELE_ID, FREQ=ALLELE_FREQ, ALPHA=ALPHA_OUT, PVALUES=P_VALUES, LOD=LOD)
	return(OUT)
end

function GWAlpha_GP(filename_sync::String, filename_phen_csv::String, MAF::Float64, DEPTH::Int64; MODEL="FIXED_RR", COVARIATE=nothing)
	###########################################################################
	### parse the sync file into allele frequencies, filter by MAF, and load ###
	############################################################################
	geno = try
		DelimitedFiles.readdlm(string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv"), ',')
	catch
		try
			rm(string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv"))
		catch
			println(string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv: does not exist!"))
		end
		sync_parsing_module.sync_parse(filename_sync); #output will be string(split(filename_sync, ".")[1], "_ALLELEFREQ.csv")
		DelimitedFiles.readdlm(string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv"), ',')
		#FORMAT: CHROM,POS,REF,ALLELE,ALLELE_FREQ_POOL1,ALLELE_FREQ_POOL2,...ALLELE_FREQ_POOLN
	end
	X_raw = convert(Array{Float64}, geno[:, 4:end]') #npools x nSNP_alleles
	n = size(X_raw)[1]
	l = size(X_raw)[2]
	LOCI = collect(1:l)[ filter_sync_module.filter_sync(filename_sync=filename_sync, MAF=MAF, DEPTH=DEPTH) ]
	X_filtered = X_raw[:, LOCI]
	### build the intercept and covariate matrix to be concatenated (hcat) to the allele count vector
	if COVARIATE == nothing
		INTERCEPT_AND_COVARIATE = ones(n, 1)
	else
		INTERCEPT_AND_COVARIATE = hcat(ones(n), COVARIATE)
	end
	### concatenate the intercept or/and  covariates with the genotype datas
	X = hcat(INTERCEPT_AND_COVARIATE, X_filtered)
	#######################################################
	### load phenotype file (GWAlpha.py-cogenompatible) ###
	#######################################################
	phen = DelimitedFiles.readdlm(filename_phen_csv, ',') #we need a new y file format with explicit pool phenotype mean without normality assumption of the distribution
	y = phen[:,end]
	#####################
	### MODEL FITTING ###
	#####################
	if MODEL == "FIXED_LS"
		#####################
		### Least Squares ###
		#####################
		b = convert(Array{Float64,1}, X' * LMM_module.inverse(X * X')  * y)
	elseif MODEL == "FIXED_RR"
		########################
		### Ridge regression ###
		########################
		# lower_limit=[0.00]
		# upper_limit=[1.0e+10]
		# initial_value=[1.00]
		# lambda = Optim.optimize(par->cost_function_RR(par, y, X), lower_limit, upper_limit, initial_value).minimizer[1] # optimse the lambda for the L2 regularisation
		# b = X' * LMM_module.inverse((X * X') + (lambda * I)) * y
		ALPHA=0.00
		GLMNET_cv = try
				GLMNet.glmnetcv(X[:,2:end], y, alpha=ALPHA)
			catch
				GLMNet.glmnetcv(X[:,2:end], y, alpha=ALPHA, nfolds=(size(X, 1)-1))
			end
		idx_lambda = argmin(GLMNET_cv.meanloss)
		lambda = GLMNET_cv.lambda[idx_lambda]
		m = GLMNet.glmnet(X[:,2:end], y, standardize=false, lambda=[lambda], alpha=ALPHA)
		b = vcat(m.a0, m.betas[:,1])
	elseif MODEL == "FIXED_GLMNET"
		##############
		### GLMNET ###
		##############
		ALPHA=0.50
		GLMNET_cv = try
				GLMNet.glmnetcv(X[:,2:end], y, alpha=ALPHA)
			catch
				GLMNet.glmnetcv(X[:,2:end], y, alpha=ALPHA, nfolds=(size(X, 1)-1))
			end
		# beta = GLMNET_cv.path.betas[:, argmin(GLMNET_cv.meanloss)]
		# b = vcat(GLMNET_cv.a0[idx_optim], beta)
		idx_lambda = argmin(GLMNET_cv.meanloss)
		lambda = GLMNET_cv.lambda[idx_lambda]
		m = GLMNet.glmnet(X[:,2:end], y, standardize=false, lambda=[lambda], alpha=ALPHA)
		b = vcat(m.a0, m.betas[:,1])
	elseif MODEL == "FIXED_LASSO"
		#############
		### LASSO ###
		#############
		ALPHA=1.00
		GLMNET_cv = try
				GLMNet.glmnetcv(X[:,2:end], y, alpha=ALPHA)
			catch
				GLMNet.glmnetcv(X[:,2:end], y, alpha=ALPHA, nfolds=(size(X, 1)-1))
			end
		# beta = GLMNET_cv.path.betas[:, argmin(GLMNET_cv.meanloss)]
		# b = vcat(GLMNET_cv.a0[idx_optim], beta)
		idx_lambda = argmin(GLMNET_cv.meanloss)
		lambda = GLMNET_cv.lambda[idx_lambda]
		m = GLMNet.glmnet(X[:,2:end], y, standardize=false, lambda=[lambda], alpha=ALPHA)
		b = vcat(m.a0, m.betas[:,1])
	elseif MODEL == "MIXED_RR"
		INTCOVAR_EFF, EFF = try
				GP_module.GGMIX(X=X_filtered, y=y, Z=COVARIATE, alfa=0.0)
			catch
				(repeat([0.0], inner=(1+size(COVARIATE, 2))), repeat([0.0], inner=size(X_filtered, 2)))
			end
		b = vcat(INTCOVAR_EFF, EFF)
	elseif MODEL == "MIXED_GLMNET"
		INTCOVAR_EFF, EFF = try
				GP_module.GGMIX(X=X_filtered, y=y, Z=COVARIATE, alfa=0.5)
			catch
				(repeat([0.0], inner=(1+size(COVARIATE, 2))), repeat([0.0], inner=size(X_filtered, 2)))
			end
		b = vcat(INTCOVAR_EFF, EFF)
	elseif MODEL == "MIXED_LASSO"
		INTCOVAR_EFF, EFF = try
				GP_module.GGMIX(X=X_filtered, y=y, Z=COVARIATE, alfa=1.0)
			catch
				(repeat([0.0], inner=(1+size(COVARIATE, 2))), repeat([0.0], inner=size(X_filtered, 2)))
			end
		b = vcat(INTCOVAR_EFF, EFF)
	else
		println(string("AwWwwW! SowWwWyyyy! We have not implemented the model: ", MODEL, " yet."))
		println("¯\\_(๑❛ᴗ❛๑)_/¯ ʚ(´◡`)ɞ")
	end

	#############
	## OUTPUT ###
	#############
	CHROM = vcat(["Intercept"], geno[:,1][LOCI])
	POS = vcat(NaN, geno[:,2][LOCI])
	ALLELE_ID = vcat(["N"], geno[:,3][LOCI])
	ALLELE_FREQ = vcat(NaN, maximum(X_raw[:, LOCI], dims=1)[1,:])
	if COVARIATE == nothing
		COVAR_EFF = nothing
		INTERCEPT = b[1]
		EFF = b[2:end]
	else
		ncovar = size(COVARIATE, 2)
		COVAR_EFF = b[2:(ncovar+1)]
		INTERCEPT = b[1]
		EFF = b[(ncovar+2):end]
	end
	P_VALUES, LOD = pval_heuristic_module.estimate_PVAL_and_LOD(EFF)
	OUT = DataFrames.DataFrame(LOCUS_ID=vcat([0], LOCI), CHROM=CHROM, POS=POS, ALLELE=ALLELE_ID, FREQ=ALLELE_FREQ, BETA=vcat([INTERCEPT], EFF), PVALUES=vcat([NaN], P_VALUES), LOD=vcat([NaN], LOD))
	# return(OUT)
	return(OUT, COVAR_EFF) ### random effects (BLUPs)
end

function write_output(OUT::DataFrames.DataFrame, filename_phe::String, MODEL::String)
	filename = basename(filename_phe)
	dirname(filename_phe) == "" ? dir = "" : dir = string(dirname(filename_phe), "/")
	if endswith(filename, ".py")
		alphas_fname = string(dir, replace(filename, ".py" => string("-", MODEL, "_Alphas.csv")))
	else
		alphas_fname = string(dir, replace(filename, ".csv" => string("-", MODEL, "_Alphas.csv")))
	end
	CSV.write(alphas_fname, OUT)
	return(alphas_fname)
end

function plot_manhattan(OUT::DataFrames.DataFrame, filename_phe::String, MODEL::String, FPR::Float64=0.01)
	### set missing values to zero
	OUT.LOD[isnan.(OUT.LOD)] .= 0.0
	### plot in R
	@rput OUT;
	@rput filename_phe;
	@rput MODEL;
	@rput FPR;
	R"OUT = droplevels(OUT[OUT$CHROM != 'Intercept', ])";
	R"OUT$LOCUS_ID = as.numeric(unlist(OUT$LOCUS_ID))";
	R"OUT$CHROM = as.factor(unlist(OUT$CHROM))";
	R"OUT$POS = as.numeric(unlist(OUT$POS))";
	R"OUT$ALLELE = as.character(unlist(OUT$ALLELE))";
	R"OUT$FREQ = as.numeric(unlist(OUT$FREQ))";
	R"OUT$PVALUES = as.numeric(unlist(OUT$PVALUES))";
	R"OUT$LOD = as.numeric(unlist(OUT$LOD))";
	R"fname_split = unlist(strsplit(basename(filename_phe), '[.]'))";
	R"title = paste(fname_split[1:(length(fname_split)-1)], collapse='[.]')";
	R"fname_png = paste0(dirname(filename_phe), '/', title, '-', MODEL, '_Manhattan.png')";
	R"contigs = levels(OUT$CHROM)";
	R"colours = RColorBrewer::brewer.pal(6, 'Paired')";
	R"colours_points = rep(colours[c(1,3)], times=ceiling(length(contigs)/2))";
	R"colours_labels = rep(colours[c(2,4)], times=ceiling(length(contigs)/2))";
	R"bonferroni_threshold = -log10(FPR/nrow(OUT))";
	### plot -log10(pvalue) points
	R"png(fname_png, width=2000, height=800)";
	R"layout(matrix(c(1,1,1,2), nrow=1, byrow=TRUE))";
	R"par(cex=1.5)";
	R"plot(x=OUT$LOCUS_ID, y=OUT$LOD, col=colours_points[as.numeric(OUT$CHROM)], ylim=c(0, max(c(OUT$LOD, bonferroni_threshold))), pch=20, main=paste0(title, '\nBonferroni Threshold at ', round(bonferroni_threshold, 2)), xaxt='n', xlab='Chromosome or Scaffold', ylab=expression(log[10]~(pvalue)), las=1)";
	R"xaxis_pos = aggregate(OUT$LOCUS_ID ~ OUT$CHROM, FUN=median)";
	R"mtext(side=1, at=xaxis_pos[,2], text=xaxis_pos[,1], cex=1.5)";
	R"abline(h=bonferroni_threshold, lty=2, lwd=2, col=colours[5])";
	R"grid()";
	# QQ plot (p-values are assumed to be uniformly distributed)
	# observed_pval = runif(100); observed_pval = observed_pval[order(observed_pval, decreasing=FALSE)] ### TEST
	# sort the observed p-values
	R"observed_pval = OUT$PVALUES[order(OUT$PVALUES, decreasing=FALSE)]";
	# calculate the pdf of the observed p-values, i.e. the probability density for a p-value interval which corresponds to each observed p-value
	R"observed_pval_density = density(observed_pval, n=length(observed_pval), from=0, to=1) ### calculate the density of the observed p-values between 0 and 1 because we are dealing with p-values which are probabilities!";
	# calculate the cummulative probabilities of the observed p-values based on the pdf: where Prob(pvalue) = Density(pvalue | pvalue_interval) * pvalue_interval
	R"observed_pval_cumprob = cumsum(observed_pval_density$y * (observed_pval_density$x[2]-observed_pval_density$x[1]))";
	# calculate the expected quantiles based on the cummulative probabilities of the observed p-values
	R"expected_pval = qunif(p=observed_pval_cumprob, min=0, max=1) ### calculate the expected quantiles based on the observed cummulative probabilities";
	# transform into -log10 scale
	R"observed_lod = -log10(observed_pval + 1e-200)";
	R"expected_lod = -log10(expected_pval + 1e-200)";
	# plot
	R"plot(x=c(min(observed_lod), max(observed_lod)), y=c(min(observed_lod), max(observed_lod)), type='n', , main='QQ Plot', xlab=expression(Expected~~log[10]~(pvalue)), ylab=expression(Observed~~log[10]~(pvalue)))";
	R"points(x=expected_lod, y=observed_lod, type='p', pch=20, col=colours[5])";
	R"lines(x=c(0, max(observed_lod)), y=c(0, max(observed_lod)), lty=2, lwd=2, col='gray')";
	R"grid()";
	R"dev.off()";
	@rget fname_png;
	return(fname_png)
end

##########################
###					   ###
### EXECUTIVE FUNCTION ###
###					   ###
##########################
"""
# __________________________________________________________________
# Genomic prediction and genome-wide association using Pool-seq data

`PoolGPAS(filename_sync::String, filename_phen::String, MAF::Float64, DEPTH::Int64; MODEL="FIXED_GWAlpha", COVARIATE=nothing)`

Build genomic prediction models and perform genome-wide association (GPAS) on quantitative traits by inferring additive genetic effects
using Pool sequencing (Pool-seq) data.

# Input
1. [synchronized pileup filename](https://sourceforge.net/p/popoolation2/wiki/Manual/)
2. phenotype data filename
- **.py** extension for iterative maximum likelihood estimation i.e. `MODEL="FIXED_GWAlpha"`, e.g.:
```
	Pheno_name='Phenotype Name';
	sig=0.06724693662723039;		# standard deviation
	MIN=0.0;						# minimum phenotype value
	MAX=0.424591738712776;			# maximum phenotype value
	perc=[0.2,0.4,0.6,0.8];			# cummulative pool sizes percentiles excluding the last pool
	q=[0.16,0.20,0.23,0.27,0.42];	# phenotype values corresponding to each percentile
```
- **.csv** extension for comma-separated headerless poolsizes and corresponding mean phenotype values, e.g.:
```
	200.0,0.11988952929875112
	200.0,0.18030259365994225
	200.0,0.21548030739673382
	200.0,0.24966378482228616
	200.0,0.31328530259365983
```
3. *MAF*: minimum allele frequency threshold (default=0.01)
4. *DEPTH*: minimum sequencing depth threshold (default=10)
5. *MODEL*: GPAS model to use (default="FIXED_GWAlpha")
	- FIXED_GWAlpha
	- FIXED_LS
	- FIXED_RR (alpha=0.0)
	- FIXED_GLMNET (alpha=0.5)
	- FIXED_LASSO (alpha=1.0)
	- MIXED_RR (alpha=0.0)
	- MIXED_GLMNET (alpha=0.5)
	- MIXED_LASSO (alpha=1.0)
6. *COVARIATE*: array of covariate/s to use (default=nothing; currently not applicable for FIXED_GWAlpha model)
7. *FPR*: False positive rate or the significance level to use to define the Bonferroni threshold (default=0.01)

# Output
1. DataFrames.DataFrame of additive allele effects with the corresponding identification (CHROM, POS, ALLELE, FREQ)
2. Array of covariate effects
3. Additive allele effects file (csv format): `string(dir, replace(filename, ".py" => string("-", MODEL, "_Alphas.csv")))` or `string(dir, replace(filename, ".csv" => string("-", MODEL, "_Alphas.csv")))`
4. Manhattan plot (png format): `string(dir, replace(filename, ".py" => string("-", MODEL, "_Manhattan.png")))` or `string(dir, replace(filename, ".csv" => string("-", MODEL, "_Manhattan.png")))`

# Examples
```
### Genome-wide association study
filename_sync = "test/test.sync"
filename_phen_py = "test/test.py"
filename_phen_csv = "test/test.csv"
@time OUT_GWAS = GWAlpha.PoolGPAS(filename_sync, filename_phen_py, MAF=0.01, DEPTH=10, MODEL="FIXED_GWAlpha", COVARIATE=nothing, FPR=0.01, PARALLEL=true)
### Genomic prediction model building
GWAlpha.poolFST_module.Fst_pairwise(sync_fname=filename_sync, window_size=100000, pool_sizes=[20,20,20,20,20], METHOD="WeirCock")
filename_covar_csv= string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_COVARIATE_FST.csv")
using DelimitedFiles
COVARIATE = DelimitedFiles.readdlm(filename_covar_csv, ',')
@time OUT_GP   = GWAlpha.PoolGPAS(filename_sync, filename_phen_csv, MAF=0.01, DEPTH=10, MODEL="FIXED_RR", COVARIATE=COVARIATE)
```
"""
function PoolGPAS(filename_sync::String, filename_phen::String; MAF::Float64=0.01, DEPTH::Int64=10, MODEL::String="FIXED_GWAlpha", COVARIATE::Any=nothing, FPR::Float64=0.01, PARALLEL::Bool=false)
	if MODEL == "FIXED_GWAlpha"
		if (PARALLEL == false)
			OUT = GWAlpha_ML(filename_sync, filename_phen, MAF)
		else
			OUT = GWAlpha_ML_parallel_module.GWAlpha_ML_parallel(filename_sync, filename_phen, MAF)
		end
		COVAR_EFF = nothing
	else
		### Load covariate if presented with a filename
		if (typeof(COVARIATE)==String)
			println(string("Please manually load your covariate file: \"", COVARIATE, "\" and use the proper delimiter."))
			println(string("Example:"))
			println(string("	using DelimitedFiles"))
			println(string("	COVARIATE = DelimitedFiles.readdlm(\"", COVARIATE, "\", ',') ### for comma-separated files or..."))
			println(string("	COVARIATE = DelimitedFiles.readdlm(\"", COVARIATE, "\", '\\t') ### for tab-delimited files"))
			return("Error!")
		else
			OUT, COVAR_EFF = GWAlpha_GP(filename_sync, filename_phen, MAF, DEPTH, MODEL=MODEL, COVARIATE=COVARIATE)
		end
	end
	### output and plotting
	alphas_fname = write_output(OUT, filename_phen, MODEL)
	plot_fname = plot_manhattan(OUT, filename_phen, MODEL, FPR)
	println("===============================================================")
	println("Everything went well. Please check the output files:")
	println(alphas_fname)
	println(plot_fname)
	println("===============================================================")
	return(OUT, COVAR_EFF)
end

end
