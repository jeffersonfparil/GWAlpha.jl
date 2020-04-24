module LMM_module

using Statistics
using LinearAlgebra
using Distributions
using Optim
using RCall

#####################
### Sub-functions ###
#####################

"""
# _______________________________________
# Polyalgorithmic matrix inverse function

`inverse(A::Array{Float64,2})`

# Inputs
1. *A* [Array{Float64,2}]: square or rectangular matrix

# Output
1. Inverse or pseudo inverse of A

# Examples
```
S = convert(Array{Float64,2}, reshape(collect(1:9), 3, 3))
LMM_module.inverse(S)
R = convert(Array{Float64,2}, reshape(collect(1:10), 5, 2))
LMM_module.inverse(R)
```
"""
function inverse(A::Array{Float64,2})
	out = try
			LinearAlgebra.inv(A)
		catch
			LinearAlgebra.pinv(A)
		end
	return(out)
end


"""
# ________________________________________________________________
# Linear mixed model parameter optimization and effects estimation

`LMM_optim(parameters::Array{Float64,1}; X::Array{Float64,2}, y::Array{Float64,1}, Z::Array{Float64,2}, METHOD_VAR_EST::String=["ML", "REML"][1], OPTIMIZE::Bool=true, VARFIXEF::Bool=false)`

The mixed linear model is defined as y = Xb + Zu + e, where:
	- X [n,p] is the centered matrix of allele frequencies
	- Z [n,n] is the square symmetric matrix of relatedness
	- y [n,1] is the centered vector of phenotypic values
	- no intercept is explicitly fitted but implicitly set at the mean phenotypic value as a consequence of centering y
	- u ~ N(0, σ²uI)
	- e ~ N(0, σ²eI)
	- y ~ N(0, V); V = (Z (σ²uI) Z') + (σ²eI)
	- variance component (σ²e, σ²u) are estimated via maximum likelihood (ML) or restricted maximum likelihood (REML)
	- fixed effects (b) are estimated by solving: (X' * inverse(V)) * inverse(X * X' * inverse(V)) * y
	- random effects (u) are estimated by solving: (σ²uI) * Z' * inverse(V) * (y - (X*b))

# Inputs
1. parameters [Array{Float64,1}]: the error variance (σ2e) and random effects variance (σ2u)
2. *X* [Array{Float64,2}; xij ∈ [0.00, 1.00]]: matrix of allele frequencies (Pool-seq data) set as the fixed effects
3. *y* [Array{Float64,1}]: vector of phenotypes
4. *Z* [Array{Any,2}]: matrix of relatedness (Fst or XX'/p where p is the number of columns in X) set as the random effects
5. *METHOD_VAR_EST* [String]: method of variance parameter estimation choose "ML" for maximum likelihood or "REML" for restricted maximum likelihood (default="ML")
6. *OPTIMIZE* [Bool]: logical value to perform optimization to find estimate the variance parameters if TRUE or to estimate the fixed and random effects if FALSE (default=true)
7. *VARFIXEF* [Bool]: logical value to estimate the variances of the fixed effects which can be computationally and memory intensive for large datasets (default=false)

# Outputs
1. *If OPTIMIZE==true*: Negative log-likelihood of the parameters (σ2e and σ2u) given the dataset (X, y, and Z)
2. *If OPTIMIZE==false*: Fixed effects (b_hat), random effects (u_hat), and the optional fixed effects variances (σ2b_hat)

# Examples
```
using Distributions
using LinearAlgebra
using Optim
n=5; l=200; Ve=1.00; Vu=1.00;
X = reshape(abs.(rand(n*l)), n, l)
Z = reshape(abs.(rand(n^2)), n, n)
b = zeros(l); b[[5,10,100]] = [1.00, 1.00, 1.00]
u = rand(Distributions.Normal(0.00, Vu), n)
e = rand(Distributions.Normal(0.00, Ve), n)
y = (X*b) .+ (Z*u) .+ e
par = Optim.optimize(par -> LMM_module.LMM_optim(par, X=X, y=y, Z=Z), [0.0,0.0], [Inf,Inf], [1.00,1.00]).minimizer
b_hat, u_hat, Vb_hat = LMM_module.LMM_optim(par, X=X, y=y, Z=Z, OPTIMIZE=false)
```
"""
function LMM_optim(parameters::Array{Float64,1}; X::Array{Float64,2}, y::Array{Float64,1}, Z::Array{Float64,2}, METHOD_VAR_EST::String=["ML", "REML"][1], OPTIMIZE::Bool=true, VARFIXEF::Bool=false)
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
	b_hat = (X' * inv_V) * inverse(X * X' * inv_V) * y
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
			push!(OUT, 0.5 * ( log(abs(det(V_new))) .+ ((y_new - intercept_new)' * inverse(V_new) * (y_new - intercept_new)) .+ (n*log(2*pi)) )[1,1]) ### negative log-likelihood of σ2e and σ2u given X, y, and Z
		else
			println(string("Sorry. ", METHOD_VAR_EST, " is not a valid method of estimating the variances of the random effects effects. Please pick ML or REML."))
			exit()
		end
	else
		################
		### ESTIMATE ###
		################
		### Return the fixed effects (b_hat), random effects (u_hat), and optional fixed effects variance (σ2b_hat) estimates
		OUT = []
		if VARFIXEF
			### Computationally intensive!!!
			Vb_hat = inverse(X' * inverse(V) * X)
			σ2b_hat = diag(Vb_hat)
		else
			σ2b_hat = nothing
		end
		# push!(OUT, (fixef=b_hat, ranef=u_hat, var_fixef=σ2b_hat))
		push!(OUT, [b_hat, u_hat, σ2b_hat])
	end
	return(OUT[1])
end

#####################
### main function ###
#####################
"""
# ______________________
# Linear Mixed Modelling

`LMM(;X::Array{Float64,2}, y::Array{Float64,1}, Z::Array{Float64,2}, METHOD_VAR_EST::String=["ML", "REML"][1])`

The mixed linear model is defined as y = Xb + Zu + e, where:
	- no intercept is fitted
	- X is centered
	- y is centered which means the intercept is automatically the mean of y
	- Z is a square symmetric matrix of relatedness
	- u ~ N(0, σ2uI)
	- e ~ N(0, σ2eI)
	- y ~ (Z (σ2uI) Z') + (σ2eI)

# Inputs
1. *X* [Array{Float64,2}; xij ∈ [0.00, 1.00]]: matrix of allele frequencies (Pool-seq data) set as the fixed effects
2. *y* [Array{Float64,1}]: vector of phenotypes
3. *Z* [Array{Any,2}]: matrix of relatedness (Fst or XX'/p where p is the number of columns in X) set as the random effects
4. *METHOD_VAR_EST* [String]: method of variance parameter estimation choose "ML" for maximum likelihood or "REML" for restricted maximum likelihood (default="ML")

# Outputs
1. Intercept (b0) the mean of the phenotype data
2. Fixed effects estimates of X (b_hat)
3. Random effects estimates of Z (u_hat)

# Examples
```
using Distributions
using LinearAlgebra
n=5; l=200; Ve=1.00; Vu=1.00;
X = reshape(abs.(rand(n*l)), n, l)
Z = reshape(abs.(rand(n^2)), n, n)
b = zeros(l); b[[5,10,100]] = [1.00, 1.00, 1.00]
u = rand(Distributions.Normal(0.00, Vu), n)
e = rand(Distributions.Normal(0.00, Ve), n)
y = (X*b) .+ (Z*u) .+ e
b0, b_hat, u_hat = LMM_module.LMM(X=X, y=y, Z=Z)
b0, b_hat, u_hat = LMM_module.LMM(X=X, y=y, Z=Z, METHOD_VAR_EST="REML")
```
"""
function LMM(;X::Array{Float64,2}, y::Array{Float64,1}, Z::Array{Float64,2}, METHOD_VAR_EST::String=["ML", "REML"][1])
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
	# optim_out = Optim.optimize(VeVu -> LMM_optim(VeVu, X=X, y=y, Z=Z, METHOD_VAR_EST=METHOD_VAR_EST, METHOD_FIX_EST=METHOD_FIX_EST, alfa=alfa, lambda=lambda, OPTIMIZE=true), lower_limit, upper_limit, initial_VeVu)
	### Optimization in Julia 1.0 is too darn slow! Going with R then! Importing Julia functions and variables into R.
	@rput inverse;
	@rput LMM_optim;
	@rput X;
	@rput y;
	@rput Z;
	@rput METHOD_VAR_EST;
	@rput lower_limit;
	@rput upper_limit;
	@rput initial_VeVu;
	R"optim_out = optim(par=initial_VeVu, fn=LMM_optim, X=X, y=y, Z=Z, METHOD_VAR_EST=METHOD_VAR_EST, OPTIMIZE=TRUE, method='L-BFGS-B', lower=lower_limit, upper=upper_limit)"
	@rget optim_out;
	Ve, Vu = optim_out[:par]
	### Fixed effects (b), random effects (u), and fixed effects variance (Vb) estimation (NOTE: No need to back-transform these estimated effects because we simply centered the explanatory variables meaning they're still on the same space)
	b, u, Vb = LMM_optim([Ve, Vu], X=X, y=y, Z=Z, METHOD_VAR_EST=METHOD_VAR_EST, OPTIMIZE=false)
	### Output
	return (b0=uy, b=b, u=u)
end

end #end of LMM module
