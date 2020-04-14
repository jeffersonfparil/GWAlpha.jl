module significance_testing_module

using DataFrames
using Distributions
using Statistics
using RCall

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
UnicodePlots.histogram(x1)
significance_testing_module.best_fitting_distribution(x1)
x2 = Distributions.rand(Distributions.Laplace(0, 1), 100)
UnicodePlots.histogram(x2)
significance_testing_module.best_fitting_distribution(x2)
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

"""
# ________________________________________________________________________________________
# Heuristic estimates of p-values and -log10(p-values) using the best-fitting distribution

`estimate_pval_lod(data::Array{Float64,1})`

# Output
1. p-values
2. -log10(p-values)

# Examples
```
using Distributions
using UnicodePlots
x1 = Distributions.rand(Distributions.Normal(0, 1), 100)
UnicodePlots.histogram(x1)
PVAL, LOD = significance_testing_module.estimate_pval_lod(x1)
UnicodePlots.scatterplot(LOD)
x2 = Distributions.rand(Distributions.Laplace(0, 1), 100)
UnicodePlots.histogram(x2)
PVAL, LOD = significance_testing_module.estimate_pval_lod(x2)
UnicodePlots.scatterplot(LOD)
```
"""
function estimate_pval_lod(data::Array{Float64,1})
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

"""
# _______________________________________________________________________________________________
# Plot -log10(p-values) with a given Bonferroni threshold at user-defined false positive rate (α)

`plot_manhattan(;gwas_out::DataFrames.DataFrame, fpr::Float64=0.01, png_fname::String=string("Manhattan_plot-", hash(rand()), ".jpeg"))`

# Inputs
1. gwas_out

# Output
1. Manhattan plot in portable network graphics (png) format

# Examples
```
using Distributions
x1 = Distributions.rand(Distributions.Normal(0, 1), 100)
PVAL, LOD = pval_heuristic_module.estimate_PVAL_and_LOD(x1)
DataFrames.DataFrame(LOCUS_IDX=vcat([0], LOCI), CHROM=CHROM, POS=POS, ALLELE=ALLELE_ID, FREQ=ALLELE_FREQ, BETA=vcat([INTERCEPT], EFF), PVALUES=vcat([NaN], P_VALUES), LOD=vcat([NaN], LOD))
plot_manhattan
```
"""
function plot_manhattan(;chrom::Array{Any,1}, pos::Array{Int64,1}, pval::Array{Float64,1}, lod::Array{Float64,1}, fpr::Float64=0.01, png_fname::String=string("Manhattan_plot-", hash(rand()), ".jpeg"))
	### plot in R
    @rput chrom;
    @rput pos;
    @rput pval;
	@rput lod;
	@rput fpr;
	@rput png_fname;
	### remove missing values and convert infinite LOD into the maximum non-infinite LOD
    R"gwas_out = data.frame(CHROM=as.factor(unlist(chrom)), POS=as.numeric(unlist(pos)), PVALUES=as.numeric(unlist(pval)), LOD=as.numeric(unlist(lod)))";
	R"gwas_out = droplevels(gwas_out[gwas_out$CHROM != 'Intercept', ])";
	R"gwas_out$LOD[is.na(gwas_out$LOD)] = 0.0";
	R"gwas_out$LOD[is.infinite(gwas_out$LOD)] = max(gwas_out$LOD[!is.infinite(gwas_out$LOD)])";
    ### arranging chromosomes or scaffold end to end
    R"gwas_out$POS_CONTINUOUS = gwas_out$POS";
    R"contigs = levels(gwas_out$CHROM)";
    R"min_pos = aggregate(POS ~ CHROM, data=gwas_out, FUN=min)"
    R"max_pos = aggregate(POS ~ CHROM, data=gwas_out, FUN=max)"
    R"idx_chrom = max_pos$CHROM[2:nrow(max_pos)]"
    R"add_us = cumsum(max_pos$POS)[1:nrow(max_pos)-1]"
    R"for (i in 1:length(idx_chrom)){
            chr = idx_chrom[i]
            add = add_us[i]
            gwas_out$POS_CONTINUOUS[gwas_out$CHROM == chr] = gwas_out$POS_CONTINUOUS[gwas_out$CHROM == chr] + add
        }"
	### coloring
	R"colours = RColorBrewer::brewer.pal(6, 'Paired')";
	R"colours_points = rep(colours[c(1,3)], times=ceiling(length(contigs)/2))";
	R"colours_labels = rep(colours[c(2,4)], times=ceiling(length(contigs)/2))";
	R"bonferroni_threshold = -log10(fpr/nrow(gwas_out))";
	### plot -log10(pvalue) points
	R"png(png_fname, width=2000, height=800)";
	R"layout(matrix(c(1,1,1,2), nrow=1, byrow=TRUE))";
	R"par(cex=1.5)";
	R"plot(x=gwas_out$POS_CONTINUOUS, y=gwas_out$LOD, col=colours_points[as.numeric(gwas_out$CHROM)], ylim=c(0, max(c(gwas_out$LOD, bonferroni_threshold))), pch=20, main=paste0(' Bonferroni Threshold at ', round(bonferroni_threshold, 2)), xaxt='n', xlab='Chromosome or Scaffold', ylab=expression(log[10]~(pvalue)), las=1)";
	R"xaxis_pos = aggregate(gwas_out$POS_CONTINUOUS ~ gwas_out$CHROM, FUN=median)";
	R"mtext(side=1, at=xaxis_pos[,2], text=xaxis_pos[,1], cex=1.5)";
	R"abline(h=bonferroni_threshold, lty=2, lwd=2, col=colours[5])";
	R"grid()";
	# QQ plot (p-values are assumed to be uniformly distributed)
	# observed_pval = runif(100); observed_pval = observed_pval[order(observed_pval, decreasing=FALSE)] ### TEST
	# sort the observed p-values
	R"observed_pval = gwas_out$PVALUES[order(gwas_out$PVALUES, decreasing=FALSE)]";
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
	return(0)
end

end # end of significance_testing_module
