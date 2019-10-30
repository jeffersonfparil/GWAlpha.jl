###############################################
###											###
### Simple iterative genomewide association ###
###											###
############################################### same module name and filename

module GWAS_module

#####################
###				  ###
### load packages ###
###				  ###
#####################
using Distributions
using Optim
using LinearAlgebra
using ProgressMeter
JULIA_SCRIPT_HOME = @__DIR__
# JULIA_SCRIPT_HOME = "/data/Lolium/Softwares/genomic_prediction/src"
push!(LOAD_PATH, JULIA_SCRIPT_HOME)
using LMM_module

# ############################################################################################
# ### TESTS
# using JLD2
# using DataFrames
# # JLD2.@save "ZZZ_TEST_INDI.jld2" id_merge X_train y_train PC_train K_train QTL_SPEC
# JLD2.@load "ZZZ_TEST_INDI.jld2"
# X_raw = X_train
# y = y_train
# MAF = 1.0 / (2.0*length(y))
# COVARIATE = K_train
# ############################################################################################

#####################
###				  ###
### main function ###
###				  ###
#####################
function GWAS(X_raw::Array{Float64,2}, y::Array{Float64,1}, MAF::Float64; COVARIATE=nothing, MODEL_TYPE="FIXED_LS", TEST="LRT")
	### MODEL_TYPE = ["FIXED_LS", "MIXED_ML", "MIXED_FAST", "MIXED_EMMAX", "MIXED_SVD"]
	### TEST = ["LRT", "WALD"]
	n = size(X_raw)[1]
	l = size(X_raw)[2]
	### filter loci by MAF
	LOC = collect(1:l)[ ((mean(X_raw, dims=1) ./ 2) .> MAF)[:] .& ((mean(X_raw, dims=1) ./ 2) .< (1.0-MAF))[:] ]
	X_filtered = X_raw[:, LOC]
	l = size(X_filtered)[2]
	### build the intercept and covariate matrix to be concatenated (hcat) to the allele count vector
	if COVARIATE == nothing
		INTERCEPT_AND_COVARIATE = ones(n, 1)
	else
		INTERCEPT_AND_COVARIATE = hcat(ones(n), COVARIATE)
	end
	### Variance components optimization method (per locus ("MIXED") vs once for all loci ("FAST_MIXED")); and
	### Phenotype tranformation and intercepts via EMMAx
	FAST_OPTIM = []
	if MODEL_TYPE == "MIXED_FAST"
		println("Calculating a common set of variance components for all loci. This will take a few minutes.")
		X = hcat(INTERCEPT_AND_COVARIATE[:,1], X_filtered) 	#fixed effects variables: intercept and genotype data
		Z = INTERCEPT_AND_COVARIATE[:,2:end]				#random effects variables: PC or K covariates
		lower_limits = [1.0e-10, 1.0e-10]
		upper_limits = [1.0e+10, 1.0e+10]
		initial_values = [1.0, 1.0] ### NOTE: setting initial values for the variances at 1.00 prevents zero slope errors during minimization 20190921
		# @time od = OnceDifferentiable(pars->LMM_module.LMM(pars, X, Z, y; OPTIMIZE=true, METHOD="ML"), initial_values)
		# @time test = optimize(od, initial_values, lower_limits, upper_limits, Fminbox{GradientDescent}())
		# @time test1 = Optim.optimize(pars->LMM_module.LMM(pars, X, Z, y; OPTIMIZE=true, METHOD="ML"), lower_limits, upper_limits, initial_values)
		# append!(FAST_OPTIM, Optim.optimize(pars->LMM_module.LMM(pars, X, Z, y; OPTIMIZE=true, METHOD="ML"), lower_limits, upper_limits, initial_values).minimizer)
		append!(FAST_OPTIM, Optim.optimize(pars->LMM_module.LMM_optim(pars, X, Z, y; METHOD="ML"), lower_limits, upper_limits, initial_values).minimizer)
	elseif MODEL_TYPE =="MIXED_EMMAX"
		println("Calculating EMMAx-derived (FOR SQUARE KINSHIP MATRIX COVARIATE) common set of variance components for all loci. This will take a few minutes.")
		X = hcat(INTERCEPT_AND_COVARIATE[:,1], X_filtered) 	#fixed effects variables: intercept and genotype data
		Z = INTERCEPT_AND_COVARIATE[:,2:end]				#random effects variables: PC or K covariates #must be square!
		lower_limits = [1.0e-10, 1.0e-10]
		upper_limits = [1.0e+10, 1.0e+10]
		initial_values = [1.0, 1.0]
		# append!(FAST_OPTIM, Optim.optimize(pars->LMM_module.LMM(pars, X, Z, y; OPTIMIZE=true, METHOD="REML_EMMAX"), lower_limits, upper_limits, initial_values).minimizer)
		append!(FAST_OPTIM, Optim.optimize(pars->LMM_module.LMM_optim(pars, X, Z, y; METHOD="REML_EMMAX"), lower_limits, upper_limits, initial_values).minimizer)
	elseif MODEL_TYPE == "MIXED_SVD"
		### NOTE: NOT WORKING! 2019-09-11
		println("Calculating a common set of variance components for all loci. This will take a few minutes.")
		X = hcat(INTERCEPT_AND_COVARIATE[:,1], X_filtered) 	#fixed effects variables: intercept ONLY
		Z = INTERCEPT_AND_COVARIATE[:,2:end]				#random effects variables: PC or K covariates
		lower_limits = [1.0e-10, 1.0e-10]
		upper_limits = [1.0e+10, 1.0e+10]
		initial_values = [1.0, 1.0]
		# σ2e, σ2u = Optim.optimize(pars->LMM_module.LMM(pars, ones(length(y)), Z, y; OPTIMIZE=true, METHOD="REML_SVD"), lower_limits, upper_limits, initial_values).minimizer
		σ2e, σ2u = Optim.optimize(pars->LMM_module.LMM_optim(pars, ones(length(y)), Z, y; METHOD="REML_SVD"), lower_limits, upper_limits, initial_values).minimizer
		U, S, Vt = LinearAlgebra.svd(Z)
		# Vs = sqrt.(1 ./ (S .+ (σ2e/(σ2e+σ2u))))
		# X_rotated = (U' * X) .* reshape(repeat(Vs, outer=size(X, 2)), size(X))
		# y_rotated = (U' * y) .* Vs
		nPC = 10
		U_r = U[:,1:nPC]
		S_r = S[1:nPC]
		Vs = sqrt.(1 ./ (S_r .+ (σ2e/(σ2e+σ2u))))
		X_rotated = (U_r' * X) .* reshape(repeat(Vs, outer=size(X,2)), (nPC, size(X,2)))
		y_rotated = (U_r' * y) .* Vs
		append!(FAST_OPTIM, [X_rotated, y_rotated, σ2e, σ2u])
	end
	### initialize output arrays
	EFF = [] ### effects of each allele
	INT_COVAR_EFF = Array{Float64}(undef, size(INTERCEPT_AND_COVARIATE,2), l) ### effects of intercept and covariates across loci x alleles
	PVAL = [] ### p-value of the effects of each allele
	println(string("Iterating across ",  l, " loci."))
	progress_bar = ProgressMeter.Progress(l, dt=1, desc="GWAS Progress: ",  barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow) #progress bar
	for i in 1:l
		x = hcat(INTERCEPT_AND_COVARIATE, X_filtered[:,i]) ### the allele x's are at the last column
		if MODEL_TYPE == "FIXED_LS"
			######################
			### LS/FIXED MODEL ###
			######################
			n = size(x, 1)
			b_hat = LMM_module.inverse(x' * x) * (x' * y)
			y_hat = x * b_hat
			y_err = y - y_hat
			Vy_err = 1/n * (y_err' * y_err)
			b_ster = abs( sqrt(Vy_err * LMM_module.inverse(x' * x))[end,end] ) #absolute value of the square root to prevent compplex numbers
			t_val = abs(b_hat[end]) / b_ster
			t_distn = Distributions.TDist(n-1)
			p_val = 2*(1.00 - Distributions.cdf(t_distn, t_val))
			append!(EFF, b_hat[end])
			append!(PVAL, p_val)
			INT_COVAR_EFF[:, i] .= b_hat[1:(end-1)]
		elseif MODEL_TYPE == "MIXED_ML"
			###################
			### MIXED MODEL ### via MAXIMUM LIKELIHOOD ESTIMATION
			###################
			X = x[:, [1, end]]	#fixed
			Z = x[:, 2:(end-1)]	#random
			n = size(X, 1)
			### variance components Vu and Ve computed for each locus
			lower_limits = [1.0e-10, 1.0e-10]
			upper_limits = [1.0e+10, 1.0e+10]
			initial_values = [0.1, 0.1]
			# var_e, var_u = Optim.optimize(parameters->LMM_module.LMM(parameters, X, Z, y; OPTIMIZE=true, METHOD="ML"), lower_limits, upper_limits, initial_values).minimizer
			var_e, var_u = Optim.optimize(parameters->LMM_module.LMM_optim(parameters, X, Z, y; METHOD="ML"), lower_limits, upper_limits, initial_values).minimizer
			# fixef, ranef, var_fixef = LMM_module.LMM([var_e, var_u], X, Z, y; OPTIMIZE=false)
			fixef, ranef, var_fixef = LMM_module.LMM_estimate([var_e, var_u], X, Z, y; VARFIXEF=true)
			if TEST == "LRT"
				### likelihood ratio test
				# L_null = exp(-LMM_module.LMM([var_e, var_u], X[:,1], Z, y; OPTIMIZE=true, METHOD="ML"))
				L_null = exp(-LMM_module.LMM_optim([var_e, var_u], X[:,1], Z, y; METHOD="ML"))
				# L_SNP = exp(-LMM_module.LMM([var_e, var_u], X, Z, y; OPTIMIZE=true, METHOD="ML"))
				L_SNP = exp(-LMM_module.LMM_optim([var_e, var_u], X, Z, y; METHOD="ML"))
				likelihood_ratio_test_statistic = 2*log(L_SNP / L_null)
				p_val = 1.0 - Distributions.cdf(Distributions.Chisq(1), likelihood_ratio_test_statistic)
			elseif TEST =="WALD"
				### Wald test
				wald_test_statistic = (fixef[end]^2) / var_fixef[end]
				p_val = 1.0 - Distributions.cdf(Distributions.Chisq(1), wald_test_statistic)
			end
			append!(EFF, fixef[end]) ### allele fixed effect/s: exclude the intercept
			append!(PVAL, p_val)
			INT_COVAR_EFF[:, i] .= vcat(fixef[1], ranef) #intercept and covariance random effects
		elseif MODEL_TYPE == "MIXED_FAST" || MODEL_TYPE == "MIXED_EMMAX"
			########################
			### FAST MIXED MODEL ### via MAXIMUM LIKELIHOOD ESTIMATION
			########################
			### AND
			#########################
			### EMMAx MIXED MODEL ### via EMMAX ALGORITHM AND RESTRICTED MAXIMUM LIKELIHOOD ESTIMATION
			#########################
			### variance components computed once for all loci (ML - MIXED_FAST) or
			### variance components computed once using Kinship matrix-transformed phenotype and intercepts (REML - MIXED_EMMAX), or
			var_e, var_u = FAST_OPTIM ### using the pre-computed variance components via ML or REML (EMMAX)
			X = x[:, [1, end]]	#fixed: intercept and ith-SNP
			Z = x[:, 2:(end-1)]	#random
			n = size(X, 1)
			# fixef, ranef, var_fixef = LMM_module.LMM([var_e, var_u], X, Z, y; OPTIMIZE=false)
			fixef, ranef, var_fixef = LMM_module.LMM_estimate([var_e, var_u], X, Z, y; VARFIXEF=true)
			if TEST == "LRT"
				### likelihood ratio test
				# L_null = exp(-LMM_module.LMM([var_e, var_u], X[:,1], Z, y; OPTIMIZE=true, METHOD="ML"))
				L_null = exp(-LMM_module.LMM_optim([var_e, var_u], X[:,1], Z, y; METHOD="ML"))
				# L_SNP = exp(-LMM_module.LMM([var_e, var_u], X, Z, y; OPTIMIZE=true, METHOD="ML"))
				L_SNP = exp(-LMM_module.LMM_optim([var_e, var_u], X, Z, y; METHOD="ML"))
				likelihood_ratio_test_statistic = 2*log(L_SNP / L_null)
				p_val = 1.0 - Distributions.cdf(Distributions.Chisq(1), likelihood_ratio_test_statistic)
			elseif TEST =="WALD"
				### Wald test
				wald_test_statistic = (fixef[end]^2) / var_fixef[end]
				p_val = 1.0 - Distributions.cdf(Distributions.Chisq(1), wald_test_statistic)
			end
			append!(EFF, fixef[end]) ### allele fixed effect/s: exclude the intercept
			append!(PVAL, p_val)
			INT_COVAR_EFF[:, i] .= vcat(fixef[1], ranef) #intercept and covariance random effects
		elseif MODEL_TYPE == "MIXED_SVD"
			#######################
			### SVD MIXED MODEL ### via LIMIX's SVD DECOMPOSITION AND RESTRICTED MAXIMUM LIKELIHOOD ESTIMATION
			#######################
			### NOTE: NOT WORKING! 2019-09-11
			### variance components computed once using Kinship matrix-transformed phenotype and intercepts (REML - MIXED_SVD)
			### the BLUPs and BLUEs solved on the rotated X and y: X_rot and y_rot
			X_rotated, y_rotated, var_e, var_u = FAST_OPTIM ### using the SVD-rotated LIMIX X and y
			n = size(X_rotated, 1)
			x_rotated = hcat(ones(n), X_rotated[:,i])
			b_hat = LMM_module.inverse(x_rotated' * x_rotated) * (x_rotated' * y_rotated)
			# y_hat = x[:, [1, end]] * b_hat
			y_hat = x_rotated * b_hat
			y_err = y_rotated - y_hat
			Vy_err = 1/n * (y_err' * y_err)
			b_ster = abs( sqrt(Vy_err * LMM_module.inverse(x_rotated' * x_rotated))[end,end] ) #absolute value of the square root to prevent compplex numbers
			t_val = abs(b_hat[end]) / b_ster
			t_distn = Distributions.TDist(n-1)
			p_val = 2*Distributions.ccdf(t_distn, t_val)
			append!(EFF, b_hat[end])
			append!(PVAL, p_val)
			INT_COVAR_EFF[:, i] .= b_hat[1:(end-1)]
		else
			println("Improper MODEL_TYPE value. Set MODEL_TYPE = {FIXED_LS, MIXED_ML, MIXED_FAST, MIXED_EMMAX}.")
			exit()
		end
		ProgressMeter.update!(progress_bar, i)
	end
	INTCOVAR_EFF = mean(INT_COVAR_EFF, dims=2) #mean of each covariate across loci since the intercept and covariance effects are ~normally distributed
	LOD = -log.(10, PVAL .+ 1.0e-20)
	return(LOC, INTCOVAR_EFF, EFF, PVAL, LOD)
end

# ############################################################################################
# ### SAMPLE EXECUTION
# @time LOCI_FIXED, INTCOVAR_EFF_FIXED, EFF_FIXED, PVAL_FIXED, LOD_FIXED = GWAS(X_train, y_train, 0.001; COVARIATE=PC_train[:, 1:10], MODEL_TYPE="FIXED_LS") ### iterative fixed model
# @time LOCI_FASTM, INTCOVAR_EFF_FASTM, EFF_FASTM, PVAL_FASTM, LOD_FASTM = GWAS(X_train, y_train, 0.001; COVARIATE=PC_train[:, 1:10], MODEL_TYPE="MIXED_FAST") ### one-time variance component estimation followed by iterative effect estimations
# @time LOCI_RLMIX, INTCOVAR_EFF_RLMIX, EFF_RLMIX, PVAL_RLMIX, LOD_RLMIX = GWAS(X_train, y_train, 0.001; COVARIATE=PC_train,          MODEL_TYPE="MIXED_EMMAX") ### one-time variance component estimation with y-transformation  | Z; NOTE: EMMAX requires a square covariate (Kinship matrix)
# @time LOCI_SVMIX, INTCOVAR_EFF_SVMIX, EFF_SVMIX, PVAL_SVMIX, LOD_SVMIX = GWAS(X_train, y_train, 0.001; COVARIATE=PC_train,          MODEL_TYPE="MIXED_SVD") ### one-time variance component estimation with y-transformation  | Z; NOTE: SVD requires a square covariate (Kinship matrix) but the Kinship-derived PC matrix do not need to be square
# # @time LOCI_MLMIX, INTCOVAR_EFF_MLMIX, EFF_MLMIX, PVAL_MLMIX, LOD_MLMIX = GWAS(X_train, y_train, 0.001; COVARIATE=PC_train[:, 1:10], MODEL_TYPE="MIXED_ML") ### iterative mixed model maximum likelihood estimation
# ### PLOT TESTS
# using Plots; Plots.pyplot()
# using ColorBrewer
# idx = []
# for i in 1:size(QTL_SPEC, 1)
# 	for j in 1:nrow(id_merge)
#     	if (QTL_SPEC.CHROM[i] == id_merge.CHROM[j]) & (QTL_SPEC.POS[i] == id_merge.POS[j])
#       		push!(idx, j)
# 	  		continue
#     	end
#  	end
# end
# num_rows = convert(Int, round(length(idx) / nrow(QTL_SPEC)))
# num_cols = nrow(QTL_SPEC)
# IDX = reshape(idx, num_rows, num_cols)
# idx_QTL = IDX[1,sortperm(QTL_SPEC.EFF_FRAC, rev=true)]
# COLORS = ColorBrewer.palette("RdYlGn", num_cols)
# FASTLM_PLOT =	Plots.plot(LOCI_FASTM, LOD_FASTM, seriestype=:scatter); #manhattan plot
# 				Plots.vline!(FASTLM_PLOT, unique(idx_QTL), line=(4, :solid, 0.6, COLORS)); #actual QTL
# 				Plots.title!(FASTLM_PLOT, "FAST-LMM")
# REMLMX_PLOT	=	Plots.plot(LOCI_RLMIX, LOD_RLMIX, seriestype=:scatter);
# 				Plots.vline!(REMLMX_PLOT, unique(idx_QTL), line=(4, :solid, 0.6, COLORS));
# 				Plots.title!(REMLMX_PLOT, "REML-EMMAX")
# SVDLMM_PLOT = 	Plots.plot(LOCI_SVMIX, LOD_SVMIX, seriestype=:scatter);
# 				Plots.vline!(SVDLMM_PLOT, unique(idx_QTL), line=(4, :solid, 0.6, COLORS));
# 				Plots.title!(SVDLMM_PLOT, "SVD-LIMIX")
# FIXEDM_PLOT = 	Plots.plot(LOCI_FIXED, LOD_FIXED, seriestype=:scatter);
# 				Plots.vline!(FIXEDM_PLOT, unique(idx_QTL), line=(4, :solid, 0.6, COLORS));
# 				Plots.title!(FIXEDM_PLOT, "FIXED MODEL")
# plot(FASTLM_PLOT, REMLMX_PLOT, SVDLMM_PLOT, FIXEDM_PLOT, layout=(4,1))
# ############################################################################################

end #end of GWAS module
