module GWAlpha

using DelimitedFiles
using DataFrames
using CSV
using Statistics

include("sync_processing_module.jl")
include("relatedness_module.jl")
include("LMM_module.jl")
include("GWAlpha_ML_module.jl")
include("significance_testing_module.jl")

"""
# __________________________________________________________________
# Genomic prediction and genome-wide association using Pool-seq data

`PoolGPAS(filename_sync::String, filename_phen::String, MAF::Float64, DEPTH::Int64; MODEL="FIXED_GWAlpha", COVARIATE=nothing)`

Build genomic prediction models and perform genome-wide association (GPAS) on quantitative traits by inferring additive genetic effects
using Pool sequencing (Pool-seq) data.

# Input
1. [synchronized pileup filename](https://sourceforge.net/p/popoolation2/wiki/Manual/)
2. phenotype data filename
- **.py** extension for iterative maximum likelihood estimation i.e. `MODEL="FIXED_GWAlpha"`, e.g.:
```
	Pheno_name='Phenotype Name';
	sig=0.06724693662723039;		# standard deviation
	MIN=0.0;						# minimum phenotype value
	MAX=0.424591738712776;			# maximum phenotype value
	perc=[0.2,0.4,0.6,0.8];			# cummulative pool sizes percentiles excluding the last pool
	q=[0.16,0.20,0.23,0.27,0.42];	# phenotype values corresponding to each percentile
```
- **.csv** extension for comma-separated headerless poolsizes and corresponding mean phenotype values, e.g.:
```
	200.0,0.11988952929875112
	200.0,0.18030259365994225
	200.0,0.21548030739673382
	200.0,0.24966378482228616
	200.0,0.31328530259365983
```
3. *MAF*: minimum allele frequency threshold (default=0.01)
4. *DEPTH*: minimum sequencing depth threshold (default=10)
5. *MODEL*: GPAS model to use (default="FIXED_GWAlpha")
	- FIXED_GWAlpha
	- FIXED_LS
	- FIXED_RR (alpha=0.0)
	- FIXED_GLMNET (alpha=0.5)
	- FIXED_LASSO (alpha=1.0)
	- MIXED_RR (alpha=0.0)
	- MIXED_GLMNET (alpha=0.5)
	- MIXED_LASSO (alpha=1.0)
6. *COVARIATE*: array of covariate/s to use (default=nothing; currently not applicable for FIXED_GWAlpha model)
7. *FPR*: False positive rate or the significance level to use to define the Bonferroni threshold (default=0.01)
8. **PARALLEL**: Parallel execution of the *FIXED_GWAlpha* model (default=false; NOTE: Load GWAlpha using @everywhere for this to work truly parallel.)

# Output
1. DataFrames.DataFrame of additive allele effects with the corresponding identification (CHROM, POS, ALLELE, FREQ)
2. Array of covariate effects
3. Additive allele effects file (csv format): `string(dir, replace(filename, ".py" => string("-", MODEL, "_Alphas.csv")))` or `string(dir, replace(filename, ".csv" => string("-", MODEL, "_Alphas.csv")))`
4. Manhattan plot (png format): `string(dir, replace(filename, ".py" => string("-", MODEL, "_Manhattan.png")))` or `string(dir, replace(filename, ".csv" => string("-", MODEL, "_Manhattan.png")))`

# Examples
```
### Genome-wide association study
filename_sync = "test/test.sync"
filename_phen_py = "test/test.py"
filename_phen_csv = "test/test.csv"

filename_sync = "/home/student.unimelb.edu.au/jparil/Downloads/test_GWAlpha.jl/INPUT/POP_03.sync"
filename_phen_py = "/home/student.unimelb.edu.au/jparil/Downloads/test_GWAlpha.jl/INPUT/POP_03.py"
filename_phen_csv = "/home/student.unimelb.edu.au/jparil/Downloads/test_GWAlpha.jl/INPUT/POP_03.csv"

@time OUT_GWAS = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_py, maf=0.001, depth=10, model="GWAlpha", fpr=0.01, plot=true)
@time OUT_ML_LS_FST = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, maf=0.001, depth=10, model="ML_LS", random_covariate="FST", fpr=0.01, plot=true)
@time OUT_ML_LS_RELATEDNESS = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, maf=0.001, depth=10, model="ML_LS", random_covariate="RELATEDNESS", fpr=0.01, plot=true)
@time OUT_ML_GLMNET_FST = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, maf=0.001, depth=10, model="ML_GLMNET", random_covariate="RELATEDNESS", glmnet_alpha=0.50, fpr=0.01, plot=true)
@time OUT_ML_GLMNET_RELATEDNESS = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, maf=0.001, depth=10, model="ML_GLMNET", random_covariate="RELATEDNESS", glmnet_alpha=0.50, fpr=0.01, plot=true)
@time OUT_REML_LS_FST = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, maf=0.001, depth=10, model="REML_LS", random_covariate="FST", fpr=0.01, plot=true)
@time OUT_REML_LS_RELATEDNESS = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, maf=0.001, depth=10, model="REML_LS", random_covariate="RELATEDNESS", fpr=0.01, plot=true)
@time OUT_REML_GLMNET_FST = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, maf=0.001, depth=10, model="REML_GLMNET", random_covariate="RELATEDNESS", glmnet_alpha=0.50, fpr=0.01, plot=true)
@time OUT_REML_GLMNET_RELATEDNESS = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, maf=0.001, depth=10, model="REML_GLMNET", random_covariate="RELATEDNESS", glmnet_alpha=0.50, fpr=0.01, plot=true)

### For parallel execution of FIXED_GWAlpha model load Distributed package, set the number of threads, and (re)load GWAlpha with @everywhere
using Distributed
Distributed.addprocs(length(Sys.cpu_info()))
@everywhere using GWAlpha
@time OUT_GWAS = ...
```
"""
function PoolGPAS(;filename_sync::String, filename_phen::String, maf::Float64=0.001, depth::Int64=1, model::String=["GWAlpha", "ML_LS", "ML_GLMNET", "REML_LS", "REML_GLMNET"][1], filename_random_covariate=nothing, random_covariate::String=["FST", "RELATEDNESS"][1], glmnet_alpha::Float64=1.00, fpr::Float64=0.01, plot::Bool=false)
	# ### TEST
	# filename_sync = "test/test.sync"
	# filename_phen = ["test/test.py", "test/test.csv"][2]
	# maf = 0.001
	# depth = 10
	# model = ["GWAlpha", "ML_LS", "ML_GLMNET", "REML_LS", "REML_GLMNET"][2]
	# filename_random_covariate = nothing
	# random_covariate = ["FST", "RELATEDNESS"][1]
	# glmnet_alpha = 1.00
	# fpr = 0.05
	if model == "GWAlpha"
		println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
		println("Performing GWAlpha: SNP-wise MLE")
		println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
		OUT = GWAlpha_ML_module.GWAlpha_ML(filename_sync=filename_sync, filename_phen_py=filename_phen, MAF=maf)
		filename_output_csv = string(join(split(filename_sync, ".")[1:(end-1)], '.'), "-GWAlpha-OUTPUT.csv")
		writedlm(filename_output_csv, reshape(vcat(string.(keys(OUT))...), 1, length(OUT)), ',') ### write header
		io = open(filename_output_csv, "a"); writedlm(io, hcat(OUT...), ','); close(io) ### append the column-binded GWAlpha output
	else
		println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
		println("Performing Mixed Modelling")
		println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
		println(string("Filtering the genotype data (sync format): \"", filename_sync, "\" by minimum allele frequency: ", maf, " and minimum depth: ", depth, "."))
		idx = sync_processing_module.sync_filter(filename_sync=filename_sync, MAF=maf, DEPTH=depth);
		p = sum(idx) ### number of predictors i.e. the total number of alleles across loci after filtering by MAF and minimum depth
		filename_sync_filtered = string(join(split(filename_sync, ".")[1:(end-1)], "."), "_MAF", maf, "_DEPTH", depth, ".sync")
		println(string("Parsing the unfiltered genotype data (sync format): \"", filename_sync, "\" parse, load, and filter using the filtering indices output."))
		sync_processing_module.sync_parse(filename_sync) ### use the unfiltered sync file and then filter with the indices (idx)
		filename_sync_parsed = string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv")
		println(string("Load the filtered and parsed genotype data (csv format): \"", filename_sync_parsed, "\"."))
		GENO = DelimitedFiles.readdlm(filename_sync_parsed, ',')[idx, :]
		X = convert(Array{Float64,2}, GENO[:,4:end]')
		println(string("Loading phenotype data (csv format): ", filename_phen))
		PHENO = DelimitedFiles.readdlm(filename_phen, ',')
		n_pools = size(PHENO,1)
		pool_sizes = convert(Array{Int64,1},PHENO[:,1])
		y = convert(Array{Float64,1},PHENO[:,2])
		println(string("Extracting loci information from \"", filename_sync_parsed, "\"."))
		CHROM = convert(Array{Any,1}, GENO[:,1])
		POS = convert(Array{Int64,1}, GENO[:,2])
		ALLELE = convert(Array{String,1}, GENO[:,3])
		FREQ = convert(Array{Float64,1}, sum(X .* (pool_sizes/sum(pool_sizes)), dims=1)[1,:])
		### load random covariate of relatedness (Fst matrix or standardized relatedness matrix) computed from the filtered sync file
		if isnothing(filename_random_covariate)
			filename_random_covariate = string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "_COVARIATE_", random_covariate, ".csv")
			println(string("Building relatedness matrix to be used as the random covariate: ", filename_random_covariate))
			if random_covariate == "FST"
				relatedness_module.Fst_pairwise(sync_fname=filename_sync_filtered, window_size=100000, pool_sizes=pool_sizes) ### default="WeirCock" [Weir and Cockerham, 1984 method](https://www.jstor.org/stable/2408641?seq=1)
			elseif random_covariate == "RELATEDNESS"
				relatedness_module.standardized_relatedness(filename_sync_filtered)
			end
		end
		println(string("Loading relatedness matrix to be used as the random covariate: ", filename_random_covariate))
		Z = DelimitedFiles.readdlm(filename_random_covariate, ',')
		println(string("Identifying the mixed model to use: ", model))
		METHOD_VAR_EST = string(split(model, "_")[1])
		METHOD_FIX_EST = string(split(model, "_")[2])
		println(string("Fitting a mixed model using ", METHOD_VAR_EST, " variance estimation and ", METHOD_FIX_EST, " fixed effects estimation with ", random_covariate, " as the random covariate."))
		b0, b_hat, u_hat = LMM_module.LMM(X=X, y=y, Z=Z, METHOD_VAR_EST=METHOD_VAR_EST, METHOD_FIX_EST=METHOD_FIX_EST, alfa=fpr)
		println("Significance testing")
		PVAL, LOD = significance_testing_module.estimate_pval_lod(b_hat)
		println("Summarizing")
		OUT = (CHROM=vcat(["Intercept"], CHROM), POS=vcat([0], POS), ALLELE=vcat(["Intercept"], ALLELE), FREQ=vcat([0.0], FREQ), ALPHA=vcat([b0], b_hat), PVALUES=vcat([0.0], PVAL), LOD=vcat([0.0], LOD))
		filename_output_csv = string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "-", model, "_", random_covariate, "-OUTPUT.csv")
		writedlm(filename_output_csv, reshape(vcat(string.(keys(OUT))...), 1, length(OUT)), ',') ### write header
		io = open(filename_output_csv, "a"); writedlm(io, hcat(OUT...), ','); close(io) ### append the column-binded GWAlpha output
		filename_ranef_csv = string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "-", model, "_", random_covariate, "-RANEF-OUTPUT.csv")
		writedlm(filename_ranef_csv, u_hat, ',') ### write headerless random effects
	end
	println("===============================================================")
	println("Everything went well. Please check the output files:")
	println(filename_output_csv)
	if model != "GWAlpha"
		println(filename_ranef_csv)
	end
	if plot
		filename_output_png = string(join(split(filename_output_csv, ".")[1:(end-1)], '.'), ".png")
		significance_testing_module.plot_manhattan(chrom=OUT.CHROM, pos=OUT.POS, pval=OUT.PVALUES, lod=OUT.LOD, fpr=fpr, png_fname=filename_output_png)
	end
	println("===============================================================")
	return(0)
end

end
