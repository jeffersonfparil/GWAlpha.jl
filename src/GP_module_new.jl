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
# using GLMNet
# using RCall
# JULIA_SCRIPT_HOME = "/home/student.unimelb.edu.au/jparil/Documents/QUANTITATIVE_GENETICS-combine_SNP_and_transcript_data/SRC/GWAlpha.jl/src"
# push!(LOAD_PATH, JULIA_SCRIPT_HOME)
# using LMM_module
# using pval_heuristic_module
include("LMM_module.jl")
include("pval_heuristic_module.jl")

# ############################################################################################
# ### Test using genomic_prediction/src/GPASim_*
# load modules and libraries
JULIA_SCRIPT_HOME = "/home/student.unimelb.edu.au/jparil/Documents/QUANTITATIVE_GENETICS-combine_SNP_and_transcript_data/SRC/GWAlpha.jl/src"
push!(LOAD_PATH, JULIA_SCRIPT_HOME)
using sync_parsing_module
using filter_sync_module
using DelimitedFiles
# setwd
DIR = "/home/student.unimelb.edu.au/jparil/Downloads/test_GWAlpha.jl"
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
	X = convert(Array{Float64,2}, GENO[Statistics.mean(GENO[:, 4:end], dims=2)[:,1] .!= 0, 4:end])'
### PHENOTYPE
	PHENO = DelimitedFiles.readdlm(fname_pheno_csv, ',')
	pool_sizes = convert(Array{Int64}, PHENO[:, 1])
	y = convert(Array{Float64}, PHENO[:, 2])
### COVARIATE
	COVARIATE = DelimitedFiles.readdlm(fname_fst_wericock, ',')
############################################################################################

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

### METHOD1: FOR OPTIMIZATION
"""
# ________________________________
# Linear mixed model cost function

`LMM_optim(parameters::Array{Float64,1}, X::Array{Float64,2}, Z::Array{Float64,2}, y::Array{Float64,1}; METHOD="ML")`

Calculates the -log-likelihood of the observed data **y** given the parameters, error variance **σ2e** and  random effects variance **σ2u**.
Used for solving the mixed model `y = Xβ + Zu + ε` where β are the fixed effects, u are the random effects, and ε are the homoscedatic residual effects.

# METHOD
- **ML**: maximum likelihood estimation
- **REML_EMMAX**: restricted maximum likelihood estimation using the EMMAX approach which requires Z to be square
- **REML_SVD**: restricted maximum likelihood estimation using sigular value decomposition

# Examples
```
using Statistics
n=10; m=100; r=3
X = convert(Array{Float64, 2}, reshape(rand([0,1], n*m), n, m))
Z1 = convert(Array{Float64, 2}, reshape(rand(collect(1:r), n*r), n, r))
Z2 = convert(Array{Float64, 2}, reshape(rand(collect(1:n), n*n), n, n))
y = rand(n)
param_test = [1.0, 1.0]
LMM_module.LMM_optim(param_test, X, Z1, y, METHOD="ML")
LMM_module.LMM_optim(param_test, X, Z2, y, METHOD="REML_EMMAX")
LMM_module.LMM_optim(param_test, X, Z1, y, METHOD="REML_SVD")
```
"""
function LMM_optim(parameters::Array{Float64,1}, X::Array{Float64,2}, Z::Array{Float64,2}, y::Array{Float64,1}; METHOD="ML")
	### y = Xb + Zu + e; where Xb=fixed and Za=random
	### METHOD = ["ML", "REML_EMMAX", "REML_SVD"]; where "ML" == "REMML_SVD"
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
	b_hat = (X' * inv_V) * inverse(X * X' * inv_V) * y
	u_hat = D * Z' * inv_V * (y - (X*b_hat))
	### Optimize the -log-likelihood function or estimate fixed effects, random effects and fixed effect estimate standard error
	OUT = []
	if METHOD == "ML" || METHOD == "REML_SVD"
	# if METHOD == "ML"
		##########################
		### maximum likelihood ###
		##########################
		mu_y = (X * b_hat)
		n = length(y)
		push!(OUT, 0.5 * ( log(abs(det(V))) + ((y - mu_y)' * inv_V * (y - mu_y)) + (n*log(2*pi)) )) ### The negative log-likelihood function y given σ2e and σ2u
	elseif METHOD == "REML_EMMAX"
		#####################################
		### restricted maximum likelihood ### via EMMAX
		##################################### #NOTE: Z MUST BE SQUARE!
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
	end
		return(OUT[1])
end

### METHOD2: FOR FIXEF AND RANEF ESIMATION
"""
# _________________________________
# Estimate fixed and random effects

`LMM_estimate(parameters::Array{Float64,1}, X::Array{Float64,2}, Z::Array{Float64,2}, y::Array{Float64,1}; VARFIXEF=true)`

Estimate the fixed and random effects given the ML- or REML-estimated parameters (error variance **σ2e** and  random effects variance **σ2u**).
With an optional fixed effects variance estimates (default `VARFIX=true`) for Wald's test.

# EXAMPLES
```
using Statistics
n=10; m=100; r=3
X = convert(Array{Float64, 2}, reshape(rand([0,1], n*m), n, m))
Z = convert(Array{Float64}, reshape(rand(collect(1:r), n*r), n, r))
y = rand(n)
params = [1.0, 1.0]
LMM_module.LMM_estimate(params, X, Z, y, VARFIXEF=true)
LMM_module.LMM_estimate(params, X, Z, y, VARFIXEF=false)
```
"""
function LMM_estimate(parameters::Array{Float64,1}, X::Array{Float64,2}, Z::Array{Float64,2}, y::Array{Float64,1}; VARFIXEF=true)
	### y = Xb + Zu + e; where Xb=fixed and Za=random
	### VARFIXEF = [true, false]; calculate variance of fixed effects (for wald's test or not because it is too computationally expensive to calculate the inverse of them huge SNP matrix!!!)
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
	b_hat = (X' * inv_V) * inverse(X * X' * inv_V) * y
	u_hat = D * Z' * inv_V * (y - (X*b_hat))
	### Optimize the -log-likelihood function or estimate fixed effects, random effects and fixed effect estimate standard error
	OUT = []
	if VARFIXEF == true
		V_b_hat = inverse(X' * inverse(V) * X)
		σ2_b_hat = diag(V_b_hat)
	else
		σ2_b_hat = nothing
	end
	push!(OUT, (fixef=b_hat, ranef=u_hat, var_fixef=σ2_b_hat))
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
function GP(;X::Array{Float64,2}, y::Array{Float64,1}, FST::Array{Float64,2}, MODEL="FIXED_LS")
### model: y = Ab + g + e; where A is atandardized X; and  g ~ N(0, FST*Vg)

### Standardize X matrix and remove correlated columns
X_corr = sum((1.00 .- abs.(cor(X))) .< 1e-5, dims=2)[:,1] .== 2
A = (X .- mean(X)) ./ std(X)

lower_limit = [1e-100, 1e-100]
upper_limit = [1000.00, 1000.00]
se, su = Optim.optimize(sesu -> LMM_optim(sesu, A, COVARIATE, y, METHOD="ML"), lower_limit, upper_limit, [1.0,1.0]).minimizer
fixef, ranef, varfixef = LMM_estimate([se, su], A, COVARIATE, y, VARFIXEF=true)
using UnicodePlots
UnicodePlots.scatterplot(fixef)
UnicodePlots.scatterplot(ranef)
# A--> random and COVARIATE --> fixed
se, su = Optim.optimize(sesu -> LMM_optim(sesu, COVARIATE, A, y, METHOD="ML"), lower_limit, upper_limit, [1.0,1.0]).minimizer
fixef, ranef, varfixef = LMM_estimate([se, su], COVARIATE, A, y, VARFIXEF=true)
UnicodePlots.scatterplot(fixef)
UnicodePlots.scatterplot(ranef)

end
end #end of GP module
