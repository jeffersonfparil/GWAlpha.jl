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
