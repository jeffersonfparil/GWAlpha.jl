#####################################
###								  ###
###   Genomic prediction models	  ###
### 			 and			  ### NOTE: MIXED MODEL VIA SVD X AND Y ROTATION HAVE TERRIBLE PREDICTIVE POWERS!!!! 20190915
### Non-iterative GWAS approaches ###
###								  ###
##################################### same module name and filename

module LMM_module

#####################
###				  ###
### load packages ###
###				  ###
#####################
using Statistics
using LinearAlgebra
using Distributions
using Optim
using RCall
R"library(glmnet)"

#####################
###				  ###
### sub functions ###
###				  ###
#####################
### polyalgorithmic matrix inverse function
function inverse(A::Array{Float64,2})
	out = try
			LinearAlgebra.inv(A)
		catch
			LinearAlgebra.pinv(A)
		end
	return(out)
end
### least squares fixed effects estimation
function BLUE_least_squares(;X::Array{Float64,2}, y::Array{Float64,1}, inv_V::Array{Float64,2})
	if isnothing(inv_V)
		b_hat = (X') * inverse(X * X') * y
	else
		b_hat = (X' * inv_V) * inverse(X * X' * inv_V) * y
	end
	return(b_hat)
end
### glmnet fixed effects estimation
function BLUE_glmnet(;X::Array{Float64,2}, y::Array{Float64,1}, inv_V::Array{Float64,2}, alfa::Float64=0.0, lambda::Float64=1.00)
	### Jerome Friedman, Trevor Hastie, Rob Tibshirani that fits entire Lasso or ElasticNet regularization paths for linear, logistic, multinomial, and Cox models using cyclic coordinate descent
	### alpha = 0.00 is L2-norm penalization (ridge) [DEFAULT]
	### alpha = 1.00 is L1-norm penalization (LASSO)
	### Rotating X to account for the random effects variances
	X = copy(X)
	if isnothing(inv_V) == false
		X[:,:] = inv_V * X
	end
	### using the lambda that generates the least error and fitting without an intercept
	@rput X;
	@rput y;
	@rput alfa;
	@rput lambda;
	R"glmnet_out = glmnet(x=X, y=y, alpha=alfa, intercept=FALSE, lambda=lambda)";
	R"b_hat = as.numeric(glmnet_out$beta)";
	# R"y_pred = X%*%b_hat"
	# R"cbind(y, y_pred)"
	# R"sqrt(mean(y-y_pred)^2)"
	# R"plot(b_hat)"
	@rget b_hat;
	return(b_hat)
end
### mixed linear model parameter optimization and effects estimation
function LMM_optim(parameters::Array{Float64,1}; X::Array{Float64,2}, y::Array{Float64,1}, Z::Array{Float64,2}, METHOD_VAR_EST::String=["ML", "REML"][1], METHOD_FIX_EST::String=["LS", "GLMNET"][1], alfa::Float64=0.0, lambda::Float64=1.00, OPTIMIZE::Bool=true, VARFIXEF::Bool=false)
	### y = Xb + Zu + e; where Xb=fixed and Za=random
	σ2e = parameters[1] # variance of the error effects (assuming homoscedasticity)
	σ2u = parameters[2] # variance of the other random effects (assuming homoscedasticity)
	n = length(y)		# number of individual samples
	l = size(X, 2)		# number of fixed effects
	q = size(Z, 2)		# number of random effects excluding error effect
	### Variance components
	R = σ2e * I					# error variance-covariance matrix
	D = σ2u * I					# random effects variance-covariance matrix
	V = (Z * D * Z') + R		# y variance-covariance matrix
	inv_V = inverse(V)
	### Mixed model equations
	###### Fixed effects:
	if METHOD_FIX_EST == "LS"
		b_hat = BLUE_least_squares(X=X, y=y, inv_V=inv_V)
	elseif METHOD_FIX_EST == "GLMNET"
		b_hat = BLUE_glmnet(X=X, y=y, inv_V=inv_V, alfa=alfa, lambda=lambda)
	else
		println(string("Sorry. ", METHOD_FIX_EST, " is not a valid method of estimating the fixed effects. Please pick LS or GLMNET."))
		exit()
	end
	###### Random effects:
	u_hat = D * Z' * inv_V * (y - (X*b_hat))
	### Calculate the log-likelihood of the parameters (variances of the error and random effects: σ2e and σ2u) given the data
	### Or estimate the fixed and random effects as well as the fixed effects variance given the optimum parameters
	if OPTIMIZE
		################
		### OPTIMIZE ###
		################
		### Calculate the -log-likelihood given the error variance (σ2e) and random effects variance (σ2u) for maximum likelihood optimization to find the best fit σ2e and σ2u
		OUT = []
		if METHOD_VAR_EST == "ML"
			mu_y = (X * b_hat)
			n = length(y)
			push!(OUT, 0.5 * ( log(abs(det(V))) + ((y - mu_y)' * inv_V * (y - mu_y)) + (n*log(2*pi)) )) ### The negative log-likelihood function y given σ2e and σ2u
		elseif METHOD_VAR_EST == "REML"
			### NOTE: Z MUST BE SQUARE!
			M = try
					inverse(LinearAlgebra.cholesky(σ2u*Z + R)) ### Cholesky decomposition
				catch
					inverse(LinearAlgebra.lu(σ2u*Z + R).L) ### LU decomposition
				end
			y_new = M' * y
			intercept_new = sum(M', dims=2) #M' * ones(n)
			V_new = M' * V * M
			n = length(y_new)
			push!(OUT, 0.5 * ( log(abs(det(V_new))) .+ ((y_new - intercept_new)' * inverse(V_new) * (y_new - intercept_new)) .+ (n*log(2*pi)) )[1,1]) ### The negative log-likelihood function y given σ2e and σ2u
		else
			println(string("Sorry. ", METHOD_VAR_EST, " is not a valid method of estimating the variances of the random effects effects. Please pick ML or REML."))
			exit()
		end
	else
		################
		### ESTIMATE ###
		################
		### Return the fixed effects (b_hat), random effects (u_hat), and optional fixed effects variance (σ2_b_hat) estimates
		OUT = []
		if VARFIXEF
			### Computationally intensive!!!
			V_b_hat = inverse(X' * inverse(V) * X)
			σ2_b_hat = diag(V_b_hat)
		else
			σ2_b_hat = nothing
		end
		# push!(OUT, (fixef=b_hat, ranef=u_hat, var_fixef=σ2_b_hat))
		push!(OUT, [b_hat, u_hat, σ2_b_hat])
	end
	return(OUT[1])
end

#####################
###				  ###
### main function ###
###				  ###
#####################
"""
# __________________
# Genomic prediction

`GP()`

Genomic prediction modelling by estimating additive allelic effects

# Input
1.
# Output
1.

# Examples
```
using Statistics
n=10; m=100
X = convert(Array{Float64, 2}, reshape(rand([0,1], n*m), n, m))
Z = convert(Array{Float64, 2}, reshape(rand(collect(1:n), n^2), n, n))
y = rand(n)
```
"""
function LMM(;X::Array{Float64,2}, y::Array{Float64,1}, Z::Array{Float64,2}, METHOD_VAR_EST::String=["ML", "REML"][1], METHOD_FIX_EST::String=["LS", "GLMNET"][1], alfa::Float64=0.0)
	### MODEL: y = Xb + Zu + e, where: u~N(0,Vu) and e~(0,Ve)
	### Generate copies of the input vector and matrices so that we don't modify the original ones
	X = copy(X); y = copy(y); Z = copy(Z);
	### Center y vector, X matrix and Z matrix (This way we don't need to estimate the intercept! It will be at zero!)
	uy = mean(y); uX = mean(X); uZ = mean(Z);
	sy = std(y); sX = std(X); sZ = std(Z);
	y[:] = (y .- uy)
	X[:,:] = (X .- uX)
	Z[:,:] = (Z .- uZ)
	### Maximum likelihood estimation of the error variance (Ve) and random covariate variance (Vu)
	# lower_limit = [0.0, 0.0]
	# upper_limit = [Inf, Inf]
	lower_limit = [1.0e-50, 1.0e-50]
	upper_limit = [1.0e+50, 1.0e+50]
	initial_VeVu = [1.0, 1.0]
	if METHOD_FIX_EST == "GLMNET"
		@rput X;
		@rput y;
		@rput alfa;
		R"lambda = cv.glmnet(x=X, y=y, alpha=alfa, intercept=FALSE)$lambda.min";
		@rget lambda; ### expediting computations with a fixed lambda for GLMNET
	else
		lambda = 1.00
	end
	# optim_out = Optim.optimize(VeVu -> LMM_optim(VeVu, X=X, y=y, Z=Z, METHOD_VAR_EST=METHOD_VAR_EST, METHOD_FIX_EST=METHOD_FIX_EST, alfa=alfa, lambda=lambda, OPTIMIZE=true), lower_limit, upper_limit, initial_VeVu)
	@rput inverse;
	@rput BLUE_least_squares;
	@rput BLUE_glmnet;
	@rput LMM_optim;
	@rput X;
	@rput y;
	@rput Z;
	@rput METHOD_VAR_EST;
	@rput METHOD_FIX_EST;
	@rput alfa;
	@rput lambda;
	@rput lower_limit;
	@rput upper_limit;
	@rput initial_VeVu;
	# R"LMM_optim(c(0.0,0.0), X=X, y=y, Z=Z, METHOD_VAR_EST=METHOD_VAR_EST, METHOD_FIX_EST=METHOD_FIX_EST, alfa=alfa, lambda=lambda, OPTIMIZE=TRUE)"
	R"optim_out = optim(par=initial_VeVu, fn=LMM_optim, X=X, y=y, Z=Z, METHOD_VAR_EST=METHOD_VAR_EST, METHOD_FIX_EST=METHOD_FIX_EST, alfa=alfa, lambda=lambda, OPTIMIZE=TRUE, method='L-BFGS-B', lower=lower_limit, upper=upper_limit)"
	@rget optim_out;
	Ve, Vu = optim_out[:par]
	### Fixed effects (b), random effects (u), and fixed effects varaince (Vb) estimation (NOTE: No need to back-transform these estimated effects because we simply centered the explanatory variables meaning they're still on the same space)
	b, u, Vb = LMM_optim([Ve, Vu], X=X, y=y, Z=Z, METHOD_VAR_EST=METHOD_VAR_EST, METHOD_FIX_EST=METHOD_FIX_EST, alfa=alfa, lambda=lambda, OPTIMIZE=false)
	### Output
	return (b0=uy, b=b, u=u)
end
end #end of LMM module

# ############################################################################################
# ### Test using genomic_prediction/src/GPASim_*
# load modules and libraries
GWAlpha_HOME = "/home/student.unimelb.edu.au/jparil/Documents/QUANTITATIVE_GENETICS-combine_SNP_and_transcript_data/SRC/GWAlpha.jl/src"
# GWAlpha_HOME = "/data/Lolium/Softwares/GWAlpha.jl/src"
cd(GWAlpha_HOME)
include("sync_parsing_module.jl")
include("filter_sync_module.jl")
using DelimitedFiles
using Statistics
using LinearAlgebra
# setwd
DIR = "/home/student.unimelb.edu.au/jparil/Downloads/test_GWAlpha.jl"
# DIR = "/data/Lolium/Quantitative_Genetics/LOLSIM2020/LOLSIM_2rep_10qtl_0.00mr_0.25fgs_0.00bgs_1grad"
cd(DIR)
# define input files
prefix = "POP_03"
fname_sync = string("INPUT/", prefix, ".sync")
fname_pheno_py = string("INPUT/", prefix, ".py")
fname_pheno_csv = string("INPUT/", prefix, ".csv")
fname_fst_hivert = string("INPUT/", prefix, "_HIVERT.fst")
fname_fst_wericock = string("INPUT/", prefix, "_WEIRCOCK.fst")
### GENOTYPE
	# filter by MAF and depth
	MAF = 0.01; DEPTH=10
	idx = filter_sync_module.filter_sync(filename_sync=fname_sync, MAF=MAF, DEPTH=DEPTH)
	fname_sync_filtered = string(join(split(fname_sync, ".")[1:(end-1)], '.'), "_MAF", MAF, "_DEPTH", DEPTH, ".sync")
	# # parse unfiltered sync
	# sync_parsing_module.sync_parse(fname_sync)
	# GENO_unfiltered = DelimitedFiles.readdlm(string(join(split(fname_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv"), ',')
	# GENO_manual_filter = GENO_unfiltered[idx, :]
	# parse filtered sync
	sync_parsing_module.sync_parse(fname_sync_filtered)
	fname_geno_csv = string(join(split(fname_sync_filtered, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv")
	# load allele frequency data
	GENO = DelimitedFiles.readdlm(fname_geno_csv, ',')
	LOCI_INFO = GENO[Statistics.mean(GENO[:, 4:end], dims=2)[:,1] .!= 0, 1:3]
	X = convert(Array{Float64,2}, GENO[Statistics.mean(GENO[:, 4:end], dims=2)[:,1] .!= 0, 4:end]')
### PHENOTYPE
	PHENO = DelimitedFiles.readdlm(fname_pheno_csv, ',')
	pool_sizes = convert(Array{Int64}, PHENO[:, 1])
	y = convert(Array{Float64}, PHENO[:, 2])
### COVARIATE
	Z = DelimitedFiles.readdlm(fname_fst_wericock, ',')
### MODEL BUILDING
	@time b0, b, u = LMM_module.LMM(X=X, y=y, Z=Z, METHOD_VAR_EST="ML", METHOD_FIX_EST="LS")
	@time b0, b, u = LMM_module.LMM(X=X, y=y, Z=Z, METHOD_VAR_EST="REML", METHOD_FIX_EST="LS")
	@time b0, b, u = LMM_module.LMM(X=X, y=y, Z=Z, METHOD_VAR_EST="ML", METHOD_FIX_EST="GLMNET", alfa=0.0)
	@time b0, b, u = LMM_module.LMM(X=X, y=y, Z=Z, METHOD_VAR_EST="REML", METHOD_FIX_EST="GLMNET", alfa=0.0)
	@time b0, b, u = LMM_module.LMM(X=X, y=y, Z=Z, METHOD_VAR_EST="ML", METHOD_FIX_EST="GLMNET", alfa=0.5)
	@time b0, b, u = LMM_module.LMM(X=X, y=y, Z=Z, METHOD_VAR_EST="REML", METHOD_FIX_EST="GLMNET", alfa=0.5)
	@time b0, b, u = LMM_module.LMM(X=X, y=y, Z=Z, METHOD_VAR_EST="ML", METHOD_FIX_EST="GLMNET", alfa=1.0)
	@time b0, b, u = LMM_module.LMM(X=X, y=y, Z=Z, METHOD_VAR_EST="REML", METHOD_FIX_EST="GLMNET", alfa=1.0)
### PREDICTION
	y_pred = b0 .+ (X*b) + (Z*u)
	RMSE = sqrt(mean((y.-y_pred).^2))
	CORR = cor(y, y_pred)
	using UnicodePlots
	UnicodePlots.scatterplot(b)
	UnicodePlots.scatterplot(u)
	UnicodePlots.scatterplot(y, y_pred)
