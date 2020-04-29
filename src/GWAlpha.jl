module GWAlpha

using DelimitedFiles
using Statistics
using RCall

include("sync_processing_module.jl")
include("relatedness_module.jl")
include("LMM_module.jl")
include("GWAlpha_ML_module.jl")
include("significance_testing_module.jl")

"""
# __________________________________________________________________
# Genomic prediction and genome-wide association using Pool-seq data

`PoolGPAS(;filename_sync::String, filename_phen::String, maf::Float64=0.001, depth::Int64=1, model::String=["GWAlpha", "MIXED", "GLMNET"][1], filename_random_covariate=nothing, random_covariate::String=["FST", "RELATEDNESS"][1], varcomp_est::String=["ML", "REML"][1], glmnet_alpha::Float64=collect(range(0.0,1.0,step=0.01,))[1], fpr::Float64=0.01, plot::Bool=false)`

Build genomic prediction models and perform genome-wide association on quantitative traits (GPAS) by inferring additive allelic effects using pool sequencing (Pool-seq) data.

# Inputs
1. *filename_sync* [String]: filename of the genotype data in [synchronized pileup format (.sync)](https://sourceforge.net/p/popoolation2/wiki/Manual/)
2. *filename_phen* [String]: filename of the phenotype data in one of two formats:
- **.py** extension for iterative maximum likelihood estimation i.e. `MODEL="GWAlpha"`, e.g.:
```
	Pheno_name='Phenotype Name';
	sig=0.06724693662723039;		# standard deviation
	MIN=0.0;						# minimum phenotype value
	MAX=0.424591738712776;			# maximum phenotype value
	perc=[0.2,0.4,0.6,0.8];			# cummulative pool sizes percentiles excluding the last pool
	q=[0.16,0.20,0.23,0.27,0.42];	# phenotype values corresponding to each percentile
```
- **.csv** extension for comma-separated headerless pool sizes and corresponding mean phenotypic values, e.g.:
```
	200.0,0.11988952929875112
	200.0,0.18030259365994225
	200.0,0.21548030739673382
	200.0,0.24966378482228616
	200.0,0.31328530259365983
```
3. *maf* [Float64]: minimum allele frequency threshold (default=0.001)
4. *depth* [Int64]: minimum sequencing depth threshold (default=1)
5. *model* [String]: model to use (default="GWAlpha")
- **"GWAlpha"** - iterative maximum likelihood estimation
- **"MIXED"** - linear mixed model
- **"GLMNET"** - [elastic-net penalization](https://web.stanford.edu/~hastie/Papers/glmnet.pdf)
6. *filename_random_covariate* [String]: filename of a precomputed headerless square symmetric matrix of pool relatedness (default=nothing)
7. *random_covariate* [String]: type of relatedness matrix to compute if `filename_random_covariate==nothing` (default="FST")
- **"FST"** - pairwise estimation of fixation indices using Pool-seq data using [Weir and Cockerham, 1984 method](https://www.jstor.org/stable/2408641?seq=1) (additionally [Hivert et al, 2018 method](https://www.biorxiv.org/content/biorxiv/early/2018/03/20/282400.full.pdf) is also available: see `?GWAlpha.relatedness_module.Fst_pairwise`)
- **"RELATEDNESS"** - simple standardized relatedness matrix `XX'/p`, where `X` is the allele frequency matrix (Pool-seq data) and `p` is the number of columns of `X`
8. *glmnet_alpha* [Float64]: elastic-net penalty (default=0.00 or ridge regression penalty)
9. *varcomp_est* [String]: variance component estimation method
- **"ML"** - maximum likelihood estimation
- **"REML"** - restricted maximum likelihood estimation
10. *glmnet_alpha* [Float64]: elastic-net penalty (from 0.00 for ridge regression up to 1.00 for LASSO)
11. *fpr* [Float64]: false positive rate for computing the Bonferroni threshold in significance testing (default=0.01)
12. *plot* [Bool]: generate Manhattan and quantile-quantile (QQ) plots, and save in a portable network graphics (.png) format (default=false)

# Outputs
1. Additive allelic effects array (header: CHROM, POS, ALLELE, FREQ, ALPHA, PVALUES, LOD) written into a comma-separated (.csv) file
- "GWAlpha": `string(join(split(filename_sync, ".")[1:(end-1)], '.'), "-GWAlpha-OUTPUT.csv")`
- "MIXED": `string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "-", model, varcomp_est, "_", random_covariate, "-OUTPUT.csv")`
- "GLMNET": `string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "-", model, "_ALPHA", glmnet_alpha, "-OUTPUT.csv")`
2. Random covariate effects vector (headerless: RANDOM_EFFECTS) written into a comma-separated (.csv) file
- "GWAlpha": nothing
- "MIXED": `string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "-", model, "_", random_covariate, "-RANEF-OUTPUT.csv")`
- "GLMNET": nothing
3. Manhattan and QQ plots in .png format
- `string(join(split(filename_output_csv, ".")[1:(end-1)], '.'), ".png")`
4. Parsing, filtering, and relatedness matrix outputs:
- Parsed sync data into a .csv file of allele frequency data (headerless: CHROM, POS, ALLELE, FREQ_POOL1, ..., FREQ_POOLn): `string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv")`
- Filtered sync data into a sync file: filename_sync_filtered = `string(join(split(filename_sync, ".")[1:(end-1)], "."), "_MAF", maf, "_DEPTH", depth, ".sync")`
- Pairwise pool relatedness (square, symmetric and headerless): `string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "_COVARIATE_", random_covariate, ".csv")`

# Examples
```
### Input files:
filename_sync = "test/test.sync"
filename_phen_py = "test/test.py"
filename_phen_csv = "test/test.csv"
### Single-threaded execution:
using GWAlpha
@time OUT_GWAS = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_py, maf=0.001, depth=10, model="GWAlpha", fpr=0.01, plot=true)
@time OUT_MIXEDML_FST = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, maf=0.001, depth=10, model="MIXED", random_covariate="FST", fpr=0.01, plot=true)
@time OUT_MIXEDREML_RELATEDNESS = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, maf=0.001, depth=10, model="MIXED", random_covariate="RELATEDNESS", varcomp_est="REML", fpr=0.01, plot=true)
@time OUT_GLMNET_RIDGE = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, maf=0.001, depth=10, model="GLMNET", glmnet_alpha=0.00, fpr=0.01, plot=true)
@time OUT_GLMNET_LASSO = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, maf=0.001, depth=10, model="GLMNET", glmnet_alpha=1.00, fpr=0.01, plot=true)
### Multi-threaded execution (parallel execution only applicable to model=="GWAlpha"):
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
@everywhere using GWAlpha
@time OUT_GWAS = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_py, maf=0.001, depth=10, model="GWAlpha", fpr=0.01, plot=true)
```
# Details
The GWAlpha model is defined as α = W(μₐₗₗₑₗₑ-μₐₗₜₑᵣₙₐₜᵢᵥₑ)/σᵧ, where:
- μ is the mean of the beta distribution Beta(θ) where θ={θ₁,θ₂}
- θ is estimated via maximum likelihood L(θ|Q) ∝ πᵢ₌₁₋ₖf(qᵢ|θ)
- Q = {q₁,...,qₖ} is the cumulative sum of allele frequencies across increasing-phenotypic-value-sorted pools where k is the number of pools
- E(allele|θ) = Beta_cdf(yᵢ',θ) - Beta_cdf(yᵢ₋₁',θ), where yᵢ' ∈ Y'
- Y' is the inverse quantile-normalized into phenotype data such that Y' ∈ [0,1]
- W = 2√{E(allele)*(1-E(allele))} is the penalization for low allele frequency

The linear mixed model is defined as y = Xb + Zu + e, where:
- X [n,p] is the centered matrix of allele frequencies
- Z [n,n] is the square symmetric matrix of relatedness
- y [n,1] is the centered vector of phenotypic values
- no intercept is explicitly fitted but implicitly set at the mean phenotypic value as a consequence of centering y
- u ~ N(0, σ²uI)
- e ~ N(0, σ²eI)
- y ~ N(0, V); V = (Z (σ²uI) Z') + (σ²eI)
- variance component (σ²e, σ²u) are estimated via maximum likelihood (ML) or restricted maximum likelihood (REML)
- fixed effects (b) are estimated via least squares: (X' * inverse(V)) * inverse(X * X' * inverse(V)) * y
- random effects (u) are estimated by solving: (σ²uI) * Z' * inverse(V) * (y - (X*b))

Empirical p-values were calculated by modelling the additive allelic effects (α) using a normal distribution with mean and variance estimated using maximum likelihood.

"""
function PoolGPAS(;filename_sync::String, filename_phen::String, maf::Float64=0.001, depth::Int64=1, model::String=["GWAlpha", "MIXED", "GLMNET"][1], filename_random_covariate=nothing, random_covariate::String=["FST", "RELATEDNESS"][1], varcomp_est::String=["ML", "REML"][1], glmnet_alpha::Float64=collect(range(0.0,1.0,step=0.01,))[1], fpr::Float64=0.01, plot::Bool=false)
	# ### TEST
	# filename_sync = "test/test.sync"
	# filename_phen = ["test/test.py", "test/test.csv"][2]
	# maf = 0.001
	# depth = 10
	# model = ["GWAlpha", "MIXED", "GLMNET"][2]
	# filename_random_covariate = nothing
	# random_covariate = ["FST", "RELATEDNESS"][1]
	# varcomp_est = "REML"
	# glmnet_alpha = 0.50
	# fpr = 0.05
	if model == "GWAlpha"
		println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
		println("PERFORMING ITERATIVE GWALPHA")
		println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
		OUT = GWAlpha_ML_module.GWAlpha_ML(filename_sync=filename_sync, filename_phen_py=filename_phen, MAF=maf)
		filename_output_csv = string(join(split(filename_sync, ".")[1:(end-1)], '.'), "-GWAlpha-OUTPUT.csv")
		writedlm(filename_output_csv, reshape(vcat(string.(keys(OUT))...), 1, length(OUT)), ',') ### write header
		io = open(filename_output_csv, "a"); writedlm(io, hcat(OUT...), ','); close(io) ### append the column-binded GWAlpha output
	else
		println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
		println("PERFORMING NON-ITERATIVE GWALPHA")
		println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
		### filter the sync file by maf and depth if the input file has not yet been filtered
		split_fname = split(filename_sync, "_")
		bool_maf = try 
						maf == parse(Float64, split(split_fname[.!isnothing.(match.(r"MAF", split_fname))][1], "MAF")[2])
					catch
						false
					end ### filtered by minimum allele frequency?
		bool_depth = try
						depth == parse(Int64, split(split(split_fname[.!isnothing.(match.(r"DEPTH", split_fname))][1], "DEPTH")[2], ".")[1])
					catch
						false
					end ### filtered by minimum sequencing depth?
		if (bool_maf & bool_depth) == false
			println(string("Filtering the genotype data (sync format): \"", filename_sync, "\" by minimum allele frequency: ", maf, " and minimum depth: ", depth, "."))
			sync_processing_module.sync_filter(filename_sync=filename_sync, MAF=maf, DEPTH=depth);
			filename_sync_filtered = string(join(split(filename_sync, ".")[1:(end-1)], "."), "_MAF", maf, "_DEPTH", depth, ".sync")
			
		else
			filename_sync_filtered = filename_sync
		end
		### parse the filtered sync file if the input file has not yet been parsed
		filename_sync_parsed = string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv")
		if !isfile(filename_sync_parsed)
			println(string("Parsing the filtered genotype data (sync format): \"", filename_sync, "\""))
			sync_processing_module.sync_parse(filename_sync_filtered)
			filename_sync_parsed = string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv")
		else
			filename_sync_parsed = filename_sync_parsed
		end
		println(string("Load the filtered and parsed genotype data (csv format): \"", filename_sync_parsed, "\"."))
		GENO = DelimitedFiles.readdlm(filename_sync_parsed, ',')
		idx_non_zero = sum(GENO[:,4:end],dims=2)[:,1] .!= 0.0
		X = convert(Array{Float64,2}, GENO[idx_non_zero, 4:end]')
		println(string("Loading phenotype data (csv format): ", filename_phen))
		PHENO = DelimitedFiles.readdlm(filename_phen, ',')
		n_pools = size(PHENO,1)
		pool_sizes = convert(Array{Int64,1},PHENO[:,1])
		y = convert(Array{Float64,1},PHENO[:,2])
		println(string("Extracting loci information from \"", filename_sync_parsed, "\"."))
		CHROM = convert(Array{Any,1}, GENO[idx_non_zero, 1])
		POS = convert(Array{Int64,1}, GENO[idx_non_zero, 2])
		ALLELE = convert(Array{String,1}, GENO[idx_non_zero, 3])
		FREQ = convert(Array{Float64,1}, sum(X .* (pool_sizes/sum(pool_sizes)'), dims=1)[1,:])
		if model == "MIXED"
			println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
			println("Performing Mixed Modelling")
			println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
			### load random covariate of relatedness (Fst matrix or standardized relatedness matrix) computed from the filtered sync file
			if isnothing(filename_random_covariate)
				filename_random_covariate = string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "_COVARIATE_", random_covariate, ".csv")
				println(string("Building relatedness matrix to be used as the random covariate: ", filename_random_covariate))
				if random_covariate == "FST"
					relatedness_module.Fst_pairwise(filename_sync=filename_sync_filtered, window_size=100000, pool_sizes=pool_sizes) ### default="WeirCock" [Weir and Cockerham, 1984 method](https://www.jstor.org/stable/2408641?seq=1)
				elseif random_covariate == "RELATEDNESS"
					relatedness_module.standardized_relatedness(filename_sync_filtered)
				end
			else
				random_covariate = "USER_DEFINED_COVARIATE"
			end
			println(string("Loading relatedness matrix to be used as the random covariate: ", filename_random_covariate))
			Z = DelimitedFiles.readdlm(filename_random_covariate, ',')
			println(string("Identifying the mixed model to use: ", model))
			METHOD_VAR_EST = string(split(model, "_")[1])
			println(string("Fitting a mixed model using ", METHOD_VAR_EST, " variance estimation with ", random_covariate, " as the random covariate."))
			b0, b_hat, u_hat = LMM_module.LMM(X=X, y=y, Z=Z, METHOD_VAR_EST=varcomp_est)
			filename_output_csv = string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "-", model, varcomp_est, "_", random_covariate, "-OUTPUT.csv")
			filename_ranef_csv = string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "-", model, varcomp_est, "_", random_covariate, "-RANEF-OUTPUT.csv")
			writedlm(filename_ranef_csv, u_hat, ',') ### write headerless random effects
		else
			println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
			println("Performing Elastic-Net Regression")
			println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
			@rput X
			@rput y
			@rput glmnet_alpha
			R"library(glmnet)";
			R"lambda = glmnet::cv.glmnet(x=X, y=y, alpha=glmnet_alpha, intercept=TRUE)$lambda.min";
			R"glmnet_out = glmnet::glmnet(x=X, y=y, alpha=glmnet_alpha, intercept=TRUE, lambda=lambda)";
			R"b0 = as.numeric(glmnet_out$a0)";
			R"b_hat = as.numeric(glmnet_out$beta)";
			@rget b0
			@rget b_hat
			filename_output_csv = string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "-", model, "_ALPHA", glmnet_alpha, "-OUTPUT.csv")
		end
		### write-out the additive allelic effects
		println("Significance testing")
		PVAL, LOD = significance_testing_module.estimate_pval_lod(b_hat)
		println("Summarizing")
		OUT = (CHROM=vcat(["Intercept"], CHROM), POS=vcat([0], POS), ALLELE=vcat(["Intercept"], ALLELE), FREQ=vcat([0.0], FREQ), ALPHA=vcat([b0], b_hat), PVALUES=vcat([0.0], PVAL), LOD=vcat([0.0], LOD))
		writedlm(filename_output_csv, reshape(vcat(string.(keys(OUT))...), 1, length(OUT)), ',') ### write header
		io = open(filename_output_csv, "a"); writedlm(io, hcat(OUT...), ','); close(io) ### append the column-binded GWAlpha output
	end
	### Manhattan plotting
	if plot
		filename_output_png = string(join(split(filename_output_csv, ".")[1:(end-1)], '.'), ".png")
		significance_testing_module.plot_manhattan(chrom=OUT.CHROM, pos=OUT.POS, pval=OUT.PVALUES, lod=OUT.LOD, fpr=fpr, png_fname=filename_output_png)
	end
	### Output reporting
	println("===============================================================")
	println("Everything went well. Please check the output files:")
	println(filename_output_csv)
	if model == "MIXED"
		println(filename_ranef_csv)
	end
	if plot
		println(filename_output_png)
	end
	println("===============================================================")
	return(0)
end

end
