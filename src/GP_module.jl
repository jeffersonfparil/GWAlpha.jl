#####################################
###								  ###
###   Genomic prediction models	  ###
### 			 and			  ### NOTE: MIXED MODEL VIA SVD X AND Y ROTATION HAVE TERRIBLE PREDICTIVE POWERS!!!! 20190915
### Non-iterative GWAS approaches ###
###								  ###
##################################### same module name and filename

module GP_module

#####################
###				  ###
### load packages ###
###				  ###
#####################
using Statistics
using LinearAlgebra
using Distributions
using Optim
using GLMNet
using RCall
# JULIA_SCRIPT_HOME = @__DIR__
# # JULIA_SCRIPT_HOME = "/data/Lolium/Softwares/genomic_prediction/src"
# push!(LOAD_PATH, JULIA_SCRIPT_HOME)
# using LMM_module
# using pval_heuristic_module
include("LMM_module.jl")
include("pval_heuristic_module.jl")

# ############################################################################################
# ### TESTS
# using JLD2
# using DataFrames
# # JLD2.@save "ZZZ_TESTS_VARIABLES.jld2" id_merge X_train y_train PC_train K_train QTL_SPEC
# JLD2.@load "ZZZ_TESTS_VARIABLES.jld2"
# X_raw = X_train
# y = y_train
# MAF = 1.0 / (2.0*length(y))
# COVARIATE = K_train
# ############################################################################################

#####################
###				  ###
### sub functions ###
###				  ###
#####################

# ### log-likelihood function of mean and sd of normally distributed effects
# function loglik_Normal(mu_sd::Array{Float64,1}, EFF::Array{Float64,1})
# 	-sum(Distributions.logpdf.(Distributions.Normal(mu_sd[1], mu_sd[2]), EFF))
# end

### using ggmix
function GGMIX(;X::Array{Float64,2}, y::Array{Float64,1}, Z::Array{Float64,2}, alfa::Float64=1.0)
	@rput X;
	@rput y;
	@rput Z;
	@rput alfa;
	R"library(ggmix)";
	R"fit <- ggmix(x=X, y=y, kinship=Z, standardize=FALSE, alpha=alfa)";
	R"hdbic <- gic(fit)";
	# ### forcing at least 1 fixed effect other than the intercept!
	# R"df_gic = data.frame(LAMBDA=hdbic$lambda, DF=hdbic$ggmix_fit$df, GIC=hdbic$gic)"
	# R"df_gic = df_gic[df_gic$DF != 1, ]"
	# R"lambda = df_gic$LAMBDA[df_gic$GIC == min(df_gic$GIC)][1]"
	# R"beta_hat = as.vector(coef(hdbic, s=lambda))";
	# R"u_hat = as.vector(ranef(hdbic, s=lambda))";
	R"beta_hat = as.vector(coef(hdbic))";
	R"u_hat = as.vector(ranef(hdbic))";
	@rget beta_hat;
	@rget u_hat;
	b0_hat = beta_hat[1];
	b_hat = beta_hat[2:end];
	INTCOVAR_EFF = vcat(b0_hat, u_hat)
	EFF = b_hat
	return(INTCOVAR_EFF, EFF)
end

#####################
###				  ###
### main function ###
###				  ###
#####################
function GP(X_raw::Array{Float64,2}, y::Array{Float64,1}, MAF::Float64; COVARIATE=nothing, MODEL="FIXED_LS")
	# ### test
	# using JLD2
	# JLD2.@load("ZZZ_TEST_ACROSS_IND_DATA.jld2")
	# ###########
	n = size(X_raw)[1]
	l = size(X_raw)[2]
	### filter loci by MAF
	LOC = collect(1:l)[ ((Statistics.mean(X_raw, dims=1) ./ 2) .> MAF)[:] .& ((Statistics.mean(X_raw, dims=1) ./ 2) .< (1.0-MAF))[:] ]
	X_filtered = X_raw[:, LOC]
	n, l = size(X_filtered)
	### build the intercept and covariate matrix to be concatenated (hcat) to the allele count vector
	if COVARIATE == nothing
		INTERCEPT_AND_COVARIATE = ones(n, 1)
	else
		INTERCEPT_AND_COVARIATE = hcat(ones(n), COVARIATE)
	end
	### concatenate the intercept or/and  covariates with the genotype datas
	X = hcat(INTERCEPT_AND_COVARIATE, X_filtered)
	if MODEL == "FIXED_LS"
		### Simple least squares estimation
		beta = X' * LMM_module.inverse(X * X') * y
		INTCOVAR_EFF = beta[1:size(INTERCEPT_AND_COVARIATE, 2)]
		EFF = beta[(size(INTERCEPT_AND_COVARIATE, 2)+1):end]
	elseif MODEL == "FIXED_RR"
		ALPHA=0.00
		GLMNET_cv = try
				GLMNet.glmnetcv(X[:,2:end], y, standardize=false, alpha=ALPHA)
			catch
				GLMNet.glmnetcv(X[:,2:end], y, standardize=false, alpha=ALPHA, nfolds=(size(X, 1)-1))
			end
		idx_lambda = argmin(GLMNET_cv.meanloss)
		lambda = GLMNET_cv.lambda[idx_lambda]
		m = GLMNet.glmnet(X[:,2:end], y, standardize=false, lambda=[lambda], alpha=ALPHA)
		beta = m.betas[:,1]
		INTCOVAR_EFF = vcat(m.a0, beta[1:(size(INTERCEPT_AND_COVARIATE, 2)-1)])
		EFF = beta[size(INTERCEPT_AND_COVARIATE, 2):end]
	elseif MODEL == "FIXED_GLMNET"
		ALPHA=0.50
		GLMNET_cv = try
				GLMNet.glmnetcv(X[:,2:end], y, standardize=false, alpha=ALPHA)
			catch
				GLMNet.glmnetcv(X[:,2:end], y, standardize=false, alpha=ALPHA, nfolds=(size(X, 1)-1))
			end
		idx_lambda = argmin(GLMNET_cv.meanloss)
		lambda = GLMNET_cv.lambda[idx_lambda]
		m = GLMNet.glmnet(X[:,2:end], y, standardize=false, lambda=[lambda], alpha=ALPHA)
		beta = m.betas[:,1]
		INTCOVAR_EFF = vcat(m.a0, beta[1:(size(INTERCEPT_AND_COVARIATE, 2)-1)])
		EFF = beta[size(INTERCEPT_AND_COVARIATE, 2):end]
	elseif MODEL == "FIXED_LASSO"
		ALPHA=1.00
		GLMNET_cv = try
				GLMNet.glmnetcv(X[:,2:end], y, standardize=false, alpha=ALPHA)
			catch
				GLMNet.glmnetcv(X[:,2:end], y, standardize=false, alpha=ALPHA, nfolds=(size(X, 1)-1))
			end
		idx_lambda = argmin(GLMNET_cv.meanloss)
		lambda = GLMNET_cv.lambda[idx_lambda]
		m = GLMNet.glmnet(X[:,2:end], y, standardize=false, lambda=[lambda], alpha=ALPHA)
		beta = m.betas[:,1]
		INTCOVAR_EFF = vcat(m.a0, beta[1:(size(INTERCEPT_AND_COVARIATE, 2)-1)])
		EFF = beta[size(INTERCEPT_AND_COVARIATE, 2):end]
	elseif MODEL == "MIXED_RR"
		INTCOVAR_EFF, EFF = try
								GGMIX(X=X_filtered, y=y, Z=COVARIATE, alfa=0.00)
							catch
								[repeat([0.0], inner=size(INTERCEPT_AND_COVARIATE, 2)), repeat([0.0], inner=size(X_filtered, 2))]
							end
	elseif MODEL == "MIXED_GLMNET"
		INTCOVAR_EFF, EFF = try
								GGMIX(X=X_filtered, y=y, Z=COVARIATE, alfa=0.50)
							catch
								[repeat([0.0], inner=size(INTERCEPT_AND_COVARIATE, 2)), repeat([0.0], inner=size(X_filtered, 2))]
							end
	elseif MODEL == "MIXED_LASSO"
		INTCOVAR_EFF, EFF = try
								GGMIX(X=X_filtered, y=y, Z=COVARIATE, alfa=1.00)
							catch
								[repeat([0.0], inner=size(INTERCEPT_AND_COVARIATE, 2)), repeat([0.0], inner=size(X_filtered, 2))]
							end
	elseif MODEL == "MIXED_LS" #mixed model least squares solved via maximum likelihood (ML)
		X_fixed = hcat(INTERCEPT_AND_COVARIATE[:,1], X_filtered)
		Z_random = INTERCEPT_AND_COVARIATE[:,2:end]
		lower_limits = [1.0e-10, 1.0e-10]
		upper_limits = [1.0e+10, 1.0e+10]
		initial_values = [0.1, 0.1]
		var_e, var_u = Optim.optimize(params->LMM_module.LMM(params, X_fixed, Z_random, y; OPTIMIZE=true, METHOD="ML"), lower_limits, upper_limits, initial_values).minimizer
		fixef, ranef, var_fixef = LMM_module.LMM([var_e, var_u], X_fixed, Z_random, y; OPTIMIZE=false, METHOD="ML", VARFIXEF=false)
		INTCOVAR_EFF = vcat(fixef[1], ranef)
		EFF = fixef[2:end]
	elseif MODEL == "MIXED_EMMAX" #mixed model EMMAX solved via restricted maximum likelihood (REML)
		X_fixed = hcat(INTERCEPT_AND_COVARIATE[:,1], X_filtered)
		Z_random = INTERCEPT_AND_COVARIATE[:,2:end]
		lower_limits = [1.0e-10, 1.0e-10]
		upper_limits = [1.0e+10, 1.0e+10]
		initial_values = [0.1, 0.1]
		var_e, var_u = Optim.optimize(params->LMM_module.LMM(params, X_fixed, Z_random, y; OPTIMIZE=true, METHOD="REML_EMMAX"), lower_limits, upper_limits, initial_values).minimizer
		fixef, ranef, var_fixef = LMM_module.LMM([var_e, var_u], X_fixed, Z_random, y; OPTIMIZE=false, METHOD="REML_EMMAX", VARFIXEF=false)
		INTCOVAR_EFF = vcat(fixef[1], ranef)
		EFF = fixef[2:end]
	# elseif MODEL == "MIXED_RR" || MODEL == "MIXED_GLMNET" || MODEL == "MIXED_LASSO"
	# 	### RR via y and X rotation
	# 	X_fixed = X_filtered 							#fixed effects variables: NO intercept and genotype data
	# 	Z_random = INTERCEPT_AND_COVARIATE[:,2:end]		#random effects variables: PC or K covariates #must be square!
	# 	lower_limits = [1.0e-10, 1.0e-10]
	# 	upper_limits = [1.0e+10, 1.0e+10]
	# 	initial_values = [0.1, 0.1]
	# 	println("Variance components estimation")
	# 	# ### SVD rotation (LMM LASSO, Rakitsch, 2013) ################################################################## start
	# 	σ2e, σ2u = Optim.optimize(pars->LMM_module.LMM(pars, ones(n), Z_random, y; OPTIMIZE=true, METHOD="ML"), lower_limits, upper_limits, initial_values).minimizer
	# 	su = σ2u / (σ2u + σ2e)
	# 	U, s, Vt = LinearAlgebra.svd(Z_random)
	# 	s_su = sqrt.( 1/(su .* s) )'
	# 	X_rotated = (U' * X_fixed ) .* reshape(repeat(s_su, outer=size(X_fixed, 2)), size(X_fixed))
	# 	y_rotated = ((U' * y) .* s_su)[:,1]
	# 	### GLMNet fitting for RR, LASSO and alpha=0.5
	# 	println("Model fitting rotated data")
	# 	if MODEL == "MIXED_RR"
	# 		GLMNET_cv = GLMNet.glmnetcv(X_rotated, y_rotated, alpha=0.00)
	# 	elseif MODEL == "MIXED_GLMNET"
	# 		GLMNET_cv = GLMNet.glmnetcv(X_rotated, y_rotated, alpha=0.50)
	# 	elseif MODEL == "MIXED_LASSO"
	# 		GLMNET_cv = GLMNet.glmnetcv(X_rotated, y_rotated, alpha=1.00)
	# 	end
	# 	beta_rotated = GLMNET_cv.path.betas[:, argmin(GLMNET_cv.meanloss)]
	# 	# println("Unrotating")
	# 	# Xhi = X_fixed .* reshape(repeat(s_su, outer=size(X_fixed, 2)), size(X_fixed))
	# 	# b0_hat, u_hat, σ2_b0_hat = LMM_module.LMM([σ2e, σ2u], ones(n, 1), Z_random, y; OPTIMIZE=false, VARFIXEF=false)
	# 	# y_pred = ( (Xhi*beta_rotated) ./ s_su) + (Z_random * u_hat)
	# 	# beta = Xhi' * LMM_module.inverse(Xhi * Xhi') * (y_pred .* s_su)
	# 	prinln("Predicting")
	# 	wu = LMM_module.inverse(Z_random) * (y - (X_fixed * beta_rotated)) ### this is greedy! Takes all effects from X and just give it to Z
	# 	y_pred = (X_fixed * beta_rotated) + (Z_random * wu)
	# 	plot(y, y_pred, seriestype=:scatter)
	# 	INTCOVAR_EFF = vcat(0.0, u_hat)
	# 	EFF = beta[:,1]
	# 	# y_pred = (X_fixed * beta) + (Z_random * u_hat)
	# 	# plot(y, y_pred, seriestype=:scatter)
	# 	### SVD rotation (LMM LASSO, Rakitsch, 2013) ################################################################## end
	else
		println("Improper MODEL value. Set MODEL = {FIXED_LS, FIXED_RR, FIXED_GLMNET, FIXED_LASSO, MIXED_RR, MIXED_GLMNET, MIXED_LASSO, MIXED_LS, MIXED_EMMAX}.")
	end
	### Estimate p-values
	# # println("Estimating p-values")
	# μ_eff, σ_eff = try
	# 				Optim.optimize(mu_sd->loglik_Normal(mu_sd, EFF), [0.0, 1.0e-10], [1.0e+10, 1.0e+10], [0.1, 0.1]).minimizer
	# 			catch
	# 				try
	# 					Optim.optimize(mu_sd->loglik_Normal(mu_sd, EFF), [0.0, 1.0e-5], [1.0e+5, 1.0e+5], [0.1, 0.1]).minimizer
	# 				catch
	# 					[0.0, 0.0] ### for no SNP effects!
	# 				end
	# 			end
	# PVAL = ones(length(LOC))
	# 	idx_less_than_mean = EFF .<= μ_eff
	# 	idx_more_than_mean = EFF .> μ_eff
	# 	PVAL[idx_less_than_mean] = Distributions.cdf.(Distributions.Normal(μ_eff, σ_eff), EFF[idx_less_than_mean])
	# 	PVAL[idx_more_than_mean] = Distributions.cdf.(Distributions.Normal(μ_eff, σ_eff), EFF[idx_more_than_mean])
	# LOD = -log10.(PVAL .+ 1.0e-20) ### avoiding infinities
	P_VALUES, LOD = pval_heuristic_module.estimate_PVAL_and_LOD(EFF)
	return(LOC, INTCOVAR_EFF, EFF, P_VALUES, LOD)
end


# ############################################################################################
# ### SAMPLE EXECUTION
# @time LOCI_LS, INTCOVAR_EFF_LS, EFF_LS, PVAL_LS, LOD_LS = GP(X_raw, y, MAF; COVARIATE=K_train, MODEL="FIXED_LS");
# @time LOCI_RR, INTCOVAR_EFF_RR, EFF_RR, PVAL_RR, LOD_RR = GP(X_raw, y, MAF; COVARIATE=K_train, MODEL="FIXED_RR");
# @time LOCI_GL, INTCOVAR_EFF_GL, EFF_GL, PVAL_GL, LOD_GL = GP(X_raw, y, MAF; COVARIATE=K_train, MODEL="FIXED_GLMNET");
# @time LOCI_LA, INTCOVAR_EFF_LA, EFF_LA, PVAL_LA, LOD_LA = GP(X_raw, y, MAF; COVARIATE=K_train, MODEL="FIXED_LASSO");
# @time LOCI_RM, INTCOVAR_EFF_RM, EFF_RM, PVAL_RM, LOD_RM = GP(X_raw, y, MAF; COVARIATE=K_train, MODEL="MIXED_RR");
# @time LOCI_GM, INTCOVAR_EFF_GM, EFF_GM, PVAL_GM, LOD_GM = GP(X_raw, y, MAF; COVARIATE=K_train, MODEL="MIXED_GLMNET");
# @time LOCI_LM, INTCOVAR_EFF_LM, EFF_LM, PVAL_LM, LOD_LM = GP(X_raw, y, MAF; COVARIATE=K_train, MODEL="MIXED_LASSO");
# # @time LOCI_ML, INTCOVAR_EFF_ML, EFF_ML, PVAL_ML, LOD_ML = GP(X_raw, y, MAF; COVARIATE=nothing, MODEL="MIXED_LS");
# # @time LOCI_ME, INTCOVAR_EFF_ME, EFF_ME, PVAL_ME, LOD_ME = GP(X_raw, y, MAF; COVARIATE=K_train, MODEL="MIXED_EMMAX");
# ### PLOT TESTS
# using Plots; Plots.pyplot()
# # using ColorBrewer
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
# # # COLORS = ColorBrewer.palette("RdYlGn", num_cols)
# PLOT_LS1 =			Plots.plot(LOCI_LS, LOD_LS, seriestype=:scatter); #manhattan plot
# 					Plots.vline!(PLOT_LS1, unique(idx_QTL), line=(4, :solid, 0.6, :green)); #actual QTL
# 					Plots.title!(PLOT_LS1, "LS (no penalisation)");
# PLOT_RR =			Plots.plot(LOCI_RR, LOD_RR, seriestype=:scatter);
# 					Plots.vline!(PLOT_RR, unique(idx_QTL), line=(4, :solid, 0.6, :green));
# 					Plots.title!(PLOT_RR, "RR (alpha=0)");
# PLOT_GL =			Plots.plot(LOCI_GL, LOD_GL, seriestype=:scatter);
# 					Plots.vline!(PLOT_GL, unique(idx_QTL), line=(4, :solid, 0.6, :green));
# 					Plots.title!(PLOT_GL, "GLMNET (alpha=0.5)");
# PLOT_LA =			Plots.plot(LOCI_LA, LOD_LA, seriestype=:scatter);
# 					Plots.vline!(PLOT_LA, unique(idx_QTL), line=(4, :solid, 0.6, :green));
# 					Plots.title!(PLOT_LA, "LASSO (alpha=1)");
#
# PLOT_LS2 =			Plots.plot(LOCI_LS, LOD_LS, seriestype=:scatter); #manhattan plot
# 					Plots.vline!(PLOT_LS2, unique(idx_QTL), line=(4, :solid, 0.6, :green)); #actual QTL
# 					Plots.title!(PLOT_LS2, "LS (no penalisation)");
# PLOT_RM =			Plots.plot(LOCI_RM, LOD_RM, seriestype=:scatter);
# 					Plots.vline!(PLOT_RM, unique(idx_QTL), line=(4, :solid, 0.6, :green));
# 					Plots.title!(PLOT_RM, "MIXED RR");
# PLOT_GM =			Plots.plot(LOCI_GM, LOD_GM, seriestype=:scatter);
# 					Plots.vline!(PLOT_GM, unique(idx_QTL), line=(4, :solid, 0.6, :green));
# 					Plots.title!(PLOT_GM, "MIXED GLMNET");
# PLOT_LM =			Plots.plot(LOCI_LM, LOD_LM, seriestype=:scatter);
# 					Plots.vline!(PLOT_LM, unique(idx_QTL), line=(4, :solid, 0.6, :green));
# 					Plots.title!(PLOT_LM, "MIXED LASSO");
# # PLOT_RM =			Plots.plot(LOCI_RM, LOD_RM, seriestype=:scatter);
# # 					Plots.vline!(PLOT_RM, unique(idx_QTL), line=(4, :solid, 0.6, :green));
# # 					Plots.title!(PLOT_RM, "RR Mixed");
# # PLOT_GM =			Plots.plot(LOCI_GM, LOD_GM, seriestype=:scatter);
# # 					Plots.vline!(PLOT_GM, unique(idx_QTL), line=(4, :solid, 0.6, :green));
# # 					Plots.title!(PLOT_GM, "GLMNET Mixed");
# # PLOT_LM =			Plots.plot(LOCI_LM, LOD_LM, seriestype=:scatter);
# # 					Plots.vline!(PLOT_LM, unique(idx_QTL), line=(4, :solid, 0.6, :green));
# # 					Plots.title!(PLOT_LM, "LASSO Mixed");
# # plot(PLOT_LS, PLOT_RR, PLOT_GL, PLOT_LA, PLOT_LS, PLOT_RM, PLOT_GM, PLOT_LM, layout=(4,2))
# plot(PLOT_LS1, PLOT_LS2, PLOT_RR, PLOT_RM, PLOT_GL, PLOT_GM, PLOT_LA, PLOT_LM, layout=(4,2))
# ############################################################################################

end #end of GWAS module


#### MISC TESTS
# function loglik_linear_mixed_model(parameters, X, Z, y; OPTIMIZE=true)
# 	σ2e = parameters[1] # variance of the error effects (assuming homoscedasticity)
# 	σ2u = parameters[2] # variance of the other random effects (assuming homoscedasticity)
# 	n = length(y)		# number of individual samples
# 	l = size(X, 2)		# number of fixed effects
# 	q = size(Z, 2)		# number of random effects excluding error effect
# 	### Variance components
# 	R = σ2e * I					# error variance-covariance matrix
# 	D = σ2u * I					# random effects variance-covariance matrix
# 	V = (Z * D * Z') + R		# y variance-covariance matrix
# 	inv_V = LMM_module.inverse(V)
# 	M = try
# 			LMM_module.inverse(LinearAlgebra.cholesky(σ2u*Z + R)) ### Cholesky decomposition
# 		catch
# 			LMM_module.inverse(LinearAlgebra.lu(σ2u*Z + R).L) ### LU decomposition
# 		end
# 	y_new = M' * y
# 	intercept_new = sum(M', dims=2) #M' * ones(n)
# 	V_new = M' * V * M
#  	loglik = 	0.5 * ( log(abs(det(V_new))) .+ ((y_new - intercept_new)' * LMM_module.inverse(V_new) * (y_new - intercept_new)) .+ (n*log(2*pi)) )[1,1] ### The negative log-likelihood function y given σ2e and σ2u
# 	if OPTIMIZE == true
# 		OUT = loglik
# 	else
# 		OUT = (y_rotated=y_new, intercept_rotated=intercept_new)
# 	end
# 	return(OUT)
# end

# ### Lambda optimisation function (minimise)
# function lambda_optim_func(lambda_vec, y, X)
# 	# b = X' * inverse((X * X') + (lambda_vec[1]*I)) * y
# 	b = LMM_module.inverse((X' * X) + (lambda_vec[1]*I)) * X' * y
# 	dy = sum((y - (X * b)).^2)/length(y) ### SS error
# 	return(dy)
# end

# ##############################################
# ### drafting a julia discourse question
# using Distributions
# using LinearAlgebra
# using Optim
# using ProgressMeter
#
# n = 100
# # p = 100*n
# p = 1000
# X = Distributions.rand([0.0, 1.0, 2.0], n, p)
# nQTL = 10
# idxQTL = sample(collect(1:p), nQTL)
# β = zeros(p)
# β[idxQTL] = rand(Distributions.Normal(2.0, 0.5), nQTL)
# # β = rand(Distributions.Chisq(1), p)
# y = (X * β) .+ rand(Distributions.Normal(0, 0.01))
#
# function glmnet_cost_func(β_vec, λ, α, y, X_withIntercept)
# 	β0 = β_vec[1]
# 	β = β_vec[2:end]
# 	X = X_withIntercept[:, 2:end]
# 	n, l = size(X)
# 	y_dev = y - (β0 .* ones(n) + (X * β))
# 	cost = 	( (1/n) * (y_dev' * y_dev) ) +
# 			( λ * ( (((1-α)/2)*(β' * β)) + (α * norm(β, 1)) ) )
# 	return(cost)
# end
# # test the glmnet_cost_func()
# b_0 = zeros(p+1);
# b_1 = ones(p+1);
# b_true = vcat([0.0], β);
# @time glmnet_cost_func(b_0, 0.5, 0.5, y, hcat(ones(n), X))
# @time glmnet_cost_func(b_1, 0.5, 0.5, y, hcat(ones(n), X))
# @time glmnet_cost_func(b_true, 0.5, 0.5, y, hcat(ones(n), X))
#
# # Nelder-Mead optimisation
# # X = copy(X_train)
# # y = copy(y_train)
# LAMBDA = 10.00;
# ALPHA = 1.00;
# function optim_layered_NelderMead(X, y, LAMBDA, ALPHA, nIter=1000)
# 	X_with_int = hcat(ones(size(X,1)), X);
# 	n, p = size(X_with_int)
# 	LS_solution = X_with_int \ y;
# 	# OPTIM = Optim.optimize(pars->glmnet_cost_func(pars, LAMBDA, ALPHA, y, X_with_int), LS_solution, NelderMead(), Optim.Options(iterations = nIter))
# 	# beta_temp = OPTIM.minimizer;
# 	# soft thresholding thing-o
# 	beta_temp = LS_solution
# 	SCALE_SD = 3.00 ### NOTE: Find a way to translate lambda ito this SD scaling factor for the thresholding!
# 	MU_beta = mean(beta_temp);
# 	SD_beta = std(beta_temp);
# 	idx = (beta_temp .> (MU_beta + (SCALE_SD*SD_beta))) .| (beta_temp .< (MU_beta - (SCALE_SD*SD_beta)));
# 	X_new = hcat(ones(n), X_with_int[:, idx]);
# 	LS_new = X_new \ y;
# 	# OPTIM = Optim.optimize(pars->glmnet_cost_func(pars, LAMBDA, ALPHA, y, X_new), LS_new, NelderMead(), Optim.Options(iterations = nIter*5))
# 	# od = Optim.OnceDifferentiable(pars->glmnet_cost_func(pars, LAMBDA, ALPHA, y, X_new), LS_new; autodiff = :forward);
# 	# OPTIM = Optim.optimize(od, LS_new, BFGS());
# 	td = Optim.TwiceDifferentiable(pars->glmnet_cost_func(pars, LAMBDA, ALPHA, y, X_new), LS_new; autodiff = :forward);
# 	OPTIM = Optim.optimize(td, LS_new, Newton());
# 	beta_hat = zeros(p)
# 	beta_hat[1] = OPTIM.minimizer[1]
# 	beta_hat[idx] = OPTIM.minimizer[2:end]
# 	INTERCEPT = beta_hat[1]
# 	SNP_EFFECTS = beta_hat[2:end]
# 	# snp_effects_fraction = SNP_EFFECTS ./ sum(SNP_EFFECTS)
# 	# OUT = SNP_EFFECTS .+ (INTERCEPT .* snp_effects_fraction)
# 	# return(OUT)
# 	return(INTERCEPT, SNP_EFFECTS)
# end
# @time INTERCEPT, SNP_EFFECTS = optim_layered_NelderMead(X, y, LAMBDA, ALPHA)
#
# #test looking for best lambda
# list_lambda = collect(1:10:100)
# nlambda = length(list_lambda)
# PREDICTORS = zeros(p+1, nlambda)
# RMSE = zeros(nlambda)
# for i in 1:nlambda
# 	println(i)
# 	INTERCEPT, SNP_EFFECTS = optim_layered_NelderMead(X, y, list_lambda[i], 1.00)
# 	PREDICTORS[:, i] .= vcat(INTERCEPT, SNP_EFFECTS)
# 	n, p = size(X)
# 	dev = y - (hcat(ones(n), X) * PREDICTORS[:,i])
# 	RMSE[i] = sqrt((dev' * dev)/n)
# end
# BETAS = PREDICTORS[:, argmin(RMSE)]
# INTERCEPT = BETAS[1]
# SNP_EFFECTS = BETAS[2:end]
# 	#stats
# 	# CORR_beta = cor(β, OUT)
# 	y_pred = (X * SNP_EFFECTS) .+ INTERCEPT;
# 	CORR_y = cor(y, y_pred)
# 	RMSE_y = sqrt((y .- y_pred)' * (y .- y_pred))/n
# 	PVAL = 2 .* Distributions.ccdf(Distributions.Normal(mean(SNP_EFFECTS), std(SNP_EFFECTS)), abs.(SNP_EFFECTS));
# 	LOD = -log.(10, PVAL);
# 	using Plots; Plots.pyplot();
# 	#### GLMNET
# 	using GLMNet
# 	@time GLMNET_cv = GLMNet.glmnetcv(X, y, alpha=1.00);
# 	beta_glmnet = GLMNET_cv.path.betas[:, argmin(GLMNET_cv.meanloss)];
# 	PVAL_glmnet = 2 .* Distributions.ccdf(Distributions.Normal(mean(beta_glmnet), std(beta_glmnet)), abs.(beta_glmnet));
# 	LOD_glmnet = -log.(10, PVAL_glmnet);
# 	y_pred_glmnet = (X * beta_glmnet) .+  GLMNET_cv.path.a0[argmin(GLMNET_cv.meanloss)];
# 	#### LS
# 	@time beta_LS = X \ y;
# 	PVAL_LS = 2 .* Distributions.ccdf(Distributions.Normal(mean(beta_LS), std(beta_LS)), abs.(beta_LS));
# 	LOD_LS = -log.(10, PVAL_LS);
# 	y_pred_LS = X * beta_LS;
# 	####
# 	p1 = plot(LOD_glmnet, seriestype=:scatter);
# 		vline!(p1, idxQTL, line=(4, :solid, 0.6, :red));
# 	p2 = plot(y, y_pred_glmnet, seriestype=:scatter);
# 	p3 = plot(LOD, seriestype=:scatter);
# 		vline!(p3, idxQTL, line=(4, :solid, 0.6, :red));
# 	p4 = plot(y, y_pred, seriestype=:scatter);
# 	p5 = plot(LOD_LS, seriestype=:scatter);
# 		vline!(p5, idxQTL, line=(4, :solid, 0.6, :red));
# 	p6 = plot(y, y_pred_LS, seriestype=:scatter);
# 	plot(p1, p2, p3, p4, p5, p6, layout=(3,2))
#
#
# # ### TESTING GRADIENT DESCENT:
# # function soft_thresholding_operator(z, gamma)
# # 	if (z > 1.0e-10) & (gamma < abs(z))
# # 		out = z - gamma
# # 	elseif (z < 1.0e-10) & (gamma < abs(z))
# # 		out = z + gamma
# # 	# elseif gamma >= abs(z)
# # 	else
# # 		out = 0
# # 	# else
# # 	# 	out = nothing
# # 	end
# # 	return(out)
# # end
# # function coordinate_descent(X_std, y, b, l, a, niter=100)
# # 	n, p = size(X_std);
# # 	b[b .<= 1.0e-10] .= 0.0;
# # 	B = Array{Float64}(undef, p, niter);
# # 	progress_bar = ProgressMeter.Progress(niter, dt=1, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
# # 	for iter in 1:niter
# # 		beta_update = zeros(p);
# # 		#iterate across loci
# # 		for j in 1:p
# # 			### COVARIANCE UPDATES: WRONGITY-WRONG!!!
# # 			idx = b .>= 1.0e-10
# # 			z = ( (1/n) * ((y' * X_std[:,j]) - ((X_std[:,j]' * X_std[:,idx]) * b[idx])) ) + b[j]
# # 			gamma = l * a
# # 			beta_update[j] = soft_thresholding_operator(z, gamma) / (1 + (l*(1-a)))
# # 			# ### RAW GRADIENT FUNCTION
# # 			# idx = collect(1:p)
# # 			# idx_comp = idx[idx .!= j]
# # 			# y_sj = X_std[:, idx_comp] * b[idx_comp]
# # 			# dCOST/dbj = -( ((y - y_sj)' * X_std[:,j])/n ) + ( (l*(1-a)*b[j])+(l*a) )
# # 		end
# # 		b = copy(beta_update)
# # 		B[:,iter] .= b
# # 		ProgressMeter.update!(progress_bar, iter)
# # 	end
# # 	return(B)
# # end
# # b = X \ y;
# # l = 100.0;
# # a = 0.5;
# # X_std = (X .- mean(X)) ./ std(X);
# # beta = coordinate_descent(X_std, y, b, l, a, 1000);
# # myPlot = plot(beta, seriestype=:scatter);
# # idx_QTL = collect(1:p)[abs.(β) .>= 10.0];
# # vline!(myPlot, idx_QTL, line=(4, :solid, 0.6, :red));
# # plot(myPlot)
# ##############################################
