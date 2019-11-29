####################################################
###											     ###
### Estimating heuristic p-values of SNP effects ###
###												 ###
####################################################

module pval_heuristic_module

#####################
### load packages ###
#####################
using DataFrames
using Distributions
using Statistics

############
### test ###
############
# using CSV
# using UnicodePlots
# data = convert(Array{Float64,1}, CSV.read("/data/Lolium/Quantitative_Genetics/GWAS_GP_2018_Inverleigh_Urana/GPAS/UG_pheno-FIXED_LS_Alphas.csv").BETA[2:end])
# data = convert(Array{Float64,1}, CSV.read("/data/Lolium/Quantitative_Genetics/GWAS_GP_2018_Inverleigh_Urana/GPAS/UG_pheno-FIXED_RR_Alphas.csv").BETA[2:end])
# UnicodePlots.histogram(data)

##############################################
### finding the best fitting distributions ###
##############################################
"""
# ____________________________
# Effects distribution fitting

`best_fitting_distribution(data::Array{Float64,1})`

Determine the best fitting distribution to model the data and heuristically estimate the p-values.

# Ditributions
- [Bernoulli](https://en.wikipedia.org/wiki/Bernoulli_distribution)
- [Beta](https://en.wikipedia.org/wiki/Beta_distribution)
- [Binomial](https://en.wikipedia.org/wiki/Binomial_distribution)
- [Categorical](https://en.wikipedia.org/wiki/Categorical_distribution)
- [Dicscrete Uniform](https://en.wikipedia.org/wiki/Discrete_uniform_distribution)
- [Exponential](https://en.wikipedia.org/wiki/Exponential_distribution)
- [Gaussian](https://en.wikipedia.org/wiki/Normal_distribution)
- [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution)
- [Geometric](https://en.wikipedia.org/wiki/Geometric_distribution)
- [Laplace](https://en.wikipedia.org/wiki/Laplace_distribution)
- [Pareto](https://en.wikipedia.org/wiki/Pareto_distribution)
- [Poisson](https://en.wikipedia.org/wiki/Poisson_distribution)
- [Inverse Gaussian or Wald](https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution)
- [Uniform](https://en.wikipedia.org/wiki/Uniform_distribution_(continuous))

# Examples
```
using Distributions
using UnicodePlots
x1 = Distributions.rand(Distributions.Normal(0, 1), 100)
histogram(x1)
pval_heuristic_module.best_fitting_distribution(x1)
x2 = Distributions.rand(Distributions.Laplace(0, 1), 100)
histogram(x2)
pval_heuristic_module.best_fitting_distribution(x2)
```
"""
function best_fitting_distribution(data::Array{Float64,1})
    DIST_NAMES =   [Distributions.Bernoulli, Distributions.Beta, Distributions.Binomial, Distributions.Categorical,
                    Distributions.DiscreteUniform, Distributions.Exponential, Distributions.Normal, Distributions.Gamma,
                    Distributions.Geometric, Distributions.Laplace, Distributions.Pareto, Distributions.Poisson,
                    Distributions.InverseGaussian, Distributions.Uniform]
    DIST_INSTANCES = [try Distributions.fit_mle(D, data); catch nothing; end for D in DIST_NAMES]
    NEG_LOGLIK = [try -sum(Distributions.logpdf.(D, data)); catch nothing; end for D in DIST_INSTANCES]
    # NEG_LOGLIK = [try sum(Distributions.pdf.(D, data)); catch nothing; end for D in DIST_INSTANCES]
    DISTRIBUTIONS_DF = DataFrames.DataFrame(NAME=DIST_NAMES[NEG_LOGLIK .!= nothing],
                                            INSTANCE=DIST_INSTANCES[NEG_LOGLIK .!= nothing],
                                            NEG_LOGLIK=convert(Array{Float64,1}, NEG_LOGLIK[NEG_LOGLIK .!= nothing]))
    D = try
        DISTRIBUTIONS_DF.INSTANCE[argmin(DISTRIBUTIONS_DF.NEG_LOGLIK)]
    catch
        nothing
    end
    return(D)
end

####################################################
### estimate p-values and LOD (-log10(p-values)) ###
####################################################
"""
# ________________________________________________________________________________________
# Heuristic estimates of p-values and -log10(p-values) using the best-fitting distribution

`estimate_PVAL_and_LOD(data)`

# Output
1. p-values
2. -log10(p-values)

# Examples
```
using Distributions
using UnicodePlots
x1 = Distributions.rand(Distributions.Normal(0, 1), 100)
histogram(x1)
PVAL, LOD = pval_heuristic_module.estimate_PVAL_and_LOD(x1)
scatterplot(LOD)
x2 = Distributions.rand(Distributions.Laplace(0, 1), 100)
histogram(x2)
PVAL, LOD = pval_heuristic_module.estimate_PVAL_and_LOD(x2)
scatterplot(LOD)
```
"""
function estimate_PVAL_and_LOD(data)
    D = best_fitting_distribution(data)
    if (D == nothing) | (std(data) == 0)
        PVAL = repeat([1.0], inner=length(data))
        LOD = PVAL .- 1.0
    else
        PVAL = [Distributions.cdf(D, x) <= Distributions.ccdf(D, x) ? 2*Distributions.cdf(D, x) : 2*Distributions.ccdf(D, x) for x in data]
        LOD = -log.(10, PVAL)
        # UnicodePlots.scatterplot(LOD)
    end
    return(PVAL, LOD)
end

end # end of pval_heuristic_module
