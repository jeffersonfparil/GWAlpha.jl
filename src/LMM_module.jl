##################################
###							   ###
### Linear Mixed Model Fitting ###
###							   ###
##################################

module LMM_module

using Statistics
using LinearAlgebra

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


end ### end of LMM_module
