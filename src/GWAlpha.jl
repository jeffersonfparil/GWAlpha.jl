module GWAlpha

#####################
### load packages ###
#####################
using DelimitedFiles
using Distributions
using LinearAlgebra
using Optim
using GLMNet
using Plots; Plots.pyplot()
using ProgressMeter
using DataFrames
using CSV
using ColorBrewer

include("sync_parsing_module.jl")
include("filter_sync_module.jl")
include("LMM_module.jl")
include("GP_module.jl")
include("pval_heuristic_module.jl")

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
		n_pools_sync = nothing #clear out contents of this redundun n_pools variable
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
			rm(string(split(filename_sync, ".")[end-1], "_ALLELEFREQ.csv"))
		catch
			sync_parsing_module.sync_parse(filename_sync); #output will be string(split(filename_sync, ".")[1], "_ALLELEFREQ.csv")
			DelimitedFiles.readdlm(string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv"), ',')
			#FORMAT: CHROM,POS,REF,ALLELE,ALLELE_FREQ_POOL1,ALLELE_FREQ_POOL2,...ALLELE_FREQ_POOLN
		end
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
		b = X' * LMM_module.inverse(X * X')  * y
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

function plot_manhattan(OUT::DataFrames.DataFrame, filename_phe::String, MODEL::String)
	### set missing values to zero
	OUT.LOD[isnan.(OUT.LOD)] .= 0.0
	### plot
	LOD_plot = Plots.plot([0, length(OUT.LOD)], [0, maximum(OUT.LOD[.!isnan.(OUT.LOD) .& .!isinf.(OUT.LOD)])],
		seriestype=:scatter,
		marker=(:circle, 5, 0.5, :White, Plots.stroke(0, :white)),
		xlabel="SNP ID",
		ylabel="LOD",
		legend=false,
		size=(1500, 400),
		xrotation=00,
		yrotation=90);
	contigs = unique(OUT.CHROM)
	nColours = length(contigs)
	# idx_remove_yellow = collect(1:9)[collect(1:9) .!= 6]
	# colours = repeat(ColorBrewer.palette("Pastel1",  nColours < 9 ? nColours : 9)[idx_remove_yellow], outer=convert(Int, ceil(nColours/8)))
	# colours_lab = repeat(ColorBrewer.palette("Set1", nColours < 9 ? nColours : 9)[idx_remove_yellow], outer=convert(Int, ceil(nColours/8)))
	colours = repeat(ColorBrewer.palette("Pastel1", 9 > nColours ? nColours : 9), outer=convert(Int, ceil(nColours/9)))
	colours_lab = repeat(ColorBrewer.palette("Set1", 9 > nColours ? nColours : 9), outer=convert(Int, ceil(nColours/9)))
	i_loc = [0, 0] # counter and locus position
	for contig in contigs
		subset_LOD = OUT[OUT.CHROM .== contig, :LOD]
		x0 = (i_loc[2]+1)
		x1 = (i_loc[2]+length(subset_LOD))
		Plots.plot!(LOD_plot, x0:x1, subset_LOD,
					seriestype=:scatter,
					marker=(:circle, 5, 0.5, colours[i_loc[1]+1], Plots.stroke(0, :white)));
		length(contigs) <= 30 ? Plots.annotate!([(x0+((x1-x0)/2), 0, text(contig, 10, colours_lab[i_loc[1]+1], :center))]) : nothing
		i_loc[1] = i_loc[1]+1
		i_loc[2] = i_loc[2]+length(subset_LOD)
	end
	bonferroni_threshold = -log10(0.05 / length(OUT.LOD))
	Plots.plot!(LOD_plot, [1,length(OUT.LOD)], [bonferroni_threshold, bonferroni_threshold], line=(:line, :dash, 0.5, 2, :red), label="Bonferroni threshold");
	filename = basename(filename_phe)
	dirname(filename_phe) == "" ? dir = "" : dir = string(dirname(filename_phe), "/")
	if endswith(filename, ".py")
		plot_fname = string(dir, replace(filename, ".py" => string("-", MODEL, "_Manhattan.png")))
	else
		plot_fname = string(dir, replace(filename, ".csv" => string("-", MODEL, "_Manhattan.png")))
	end
	Plots.savefig(plot_fname)
	return(plot_fname)
end

##########################
###					   ###
### EXECUTIVE FUNCTION ###
###					   ###
##########################
function PoolGPAS(filename_sync::String, filename_phen::String, MAF::Float64, DEPTH::Int64; MODEL="FIXED_GWAlpha", COVARIATE=nothing)
	# ################################
	# ### TESTS:
	# filename_sync = "test_1kloci_g1000_p01_POOLS_GENO.sync"
	# filename_phen = "test_1kloci_g1000_p01_POOLS_PHENO.py"
	# MAF = 0.01
	# MODEL = "GWAlpha"
	# COVARIATE = nothing
	# ################################

	if MODEL == "FIXED_GWAlpha"
		OUT = GWAlpha_ML(filename_sync, filename_phen, MAF)
		COVAR_EFF = nothing
	else
		OUT, COVAR_EFF = GWAlpha_GP(filename_sync, filename_phen, MAF, DEPTH, MODEL=MODEL, COVARIATE=COVARIATE)
	end

	## output and plotting
	alphas_fname = write_output(OUT, filename_phen, MODEL)
	plot_fname = plot_manhattan(OUT, filename_phen, MODEL)
	println("===============================================================")
	println("Everything went well. Please check the output files:")
	println(alphas_fname)
	println(plot_fname)
	println("===============================================================")
	return(OUT, COVAR_EFF)
end

end
