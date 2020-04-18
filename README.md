# GWAlpha.jl

|                                                          **Lab Website**                                                          |                                                            **Build Status**                                                             |                                                                             **Documentation**                                                                             |
|:---------------------------------------------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| <a href="https://adaptive-evolution.biosciences.unimelb.edu.au/"><img src="misc/Adaptive Evolution Logo mod.png" width="150"></a> | [![Build Status](https://travis-ci.com/jeffersonfparil/GWAlpha.jl.svg?branch=master)](https://travis-ci.com/jeffersonfparil/GWAlpha.jl) | <a href="https://github.com/jeffersonfparil/GWAlpha.jl/wiki" target="_blank"><img src="https://img.shields.io/badge/docs-latest-blue.svg" alt="Latest documentation"></a> |

<!--- [![CircleCI](https://circleci.com/gh/jeffersonfparil/GWAlpha.svg?style=shield)](https://circleci.com/gh/jeffersonfparil/GWAlpha) --->

A [Julia](https://julialang.org/downloads/) package for building genomic prediction models and performing genome-wide association on quantitative traits by inferring additive allelic effects using pool sequencing data (Pool-seq; i.e. allele frequencies).

## Installation
Install dependencies (see [.travis.yml](https://github.com/jeffersonfparil/GWAlpha.jl/tree/master/.travis.yml)) for Ubuntu18.04-specific installation:
- [Julia](https://julialang.org/downloads/)
- [R](https://www.r-project.org/)
- R packages:
  - [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html)
	- [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)

Installation in Julia:
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/jeffersonfparil/GWAlpha.jl.git", rev="master"))
using GWAlpha
```

## Usage
```julia
GWAlpha.PoolGPAS(;filename_sync::String, filename_phen::String, maf::Float64=0.001, depth::Int64=1, model::String=["GWAlpha", "ML_LS", "ML_GLMNET", "REML_LS", "REML_GLMNET"][1], filename_random_covariate=nothing, random_covariate::String=["FST", "RELATEDNESS"][1], glmnet_alpha::Float64=collect(range(0.0,1.0,step=0.01,))[1], fpr::Float64=0.01, plot::Bool=false)
```

## Inputs
1. *filename_sync* [String]: filename of the genotype data in [synchronized pileup format (.sync)](https://sourceforge.net/p/popoolation2/wiki/Manual/)
2. *filename_phen* [String]: filename of the phenotype data in one of two formats:
- **.py** extension for iterative maximum likelihood estimation i.e. `MODEL="GWAlpha"`, e.g.:
```julia
	Pheno_name='Phenotype Name';
	sig=0.06724693662723039;		# standard deviation
	MIN=0.0;				# minimum phenotype value
	MAX=0.424591738712776;			# maximum phenotype value
	perc=[0.2,0.4,0.6,0.8];			# cummulative pool sizes percentiles excluding the last pool
	q=[0.16,0.20,0.23,0.27,0.42];		# phenotype values corresponding to each percentile
```
- **.csv** extension for comma-separated headerless pool sizes and corresponding mean phenotypic values, e.g.:
```julia
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
	- **"ML_LS"** - linear mixed model using maximum likelihood (ML) estimation of variances and least squares (LS) estimation of fixed effects (additive allelic effects)
	- **"ML_GLMNET"** - linear mixed model using ML and the elastic-net penalization (GLMNET) to estimate the additive allelic effects
	- **"REML_LS"** - linear mixed model using restricted maximum likelihood (REML) and LS
	- **"REML_GLMNET"** - linear mixed model using REML and GLMNET
6. *filename_random_covariate* [String]: filename of a precomputed headerless square symmetric matrix of pool relatedness (default=nothing)
7. *random_covariate* [String]: type of relatedness matrix to compute if `filename_random_covariate==nothing` (default="FST")
	- **"FST"** - pairwise estimation of fixation indices using Pool-seq data using [Weir and Cockerham, 1984 method](https://www.jstor.org/stable/2408641?seq=1) (additionally [Hivert et al, 2018 method](https://www.biorxiv.org/content/biorxiv/early/2018/03/20/282400.full.pdf) is also available: see `?GWAlpha.relatedness_module.Fst_pairwise`)
	- **"RELATEDNESS"** - simple standardized relatedness matrix `XX'/p`, where `X` is the allele frequency matrix (Pool-seq data) and `p` is the number of columns of `X`
8. *glmnet_alpha* [Float64]: elastic-net penalty (default=0.00 or ridge regression penalty)
9. *fpr* [Float64]: false positive rate for computing the Bonferroni threshold in significance testing (default=0.01)
9. *plot* [Bool]: generate Manhattan and quantile-quantile (QQ) plots, and save in a portable network graphics (.png) format (default=false)

## Outputs
1. Additive allelic effects array (header: CHROM, POS, ALLELE, FREQ, ALPHA, PVALUES, LOD) written into a comma-separated (.csv) file
	- "GWAlpha": `string(join(split(filename_sync, ".")[1:(end-1)], '.'), "-GWAlpha-OUTPUT.csv")`
	- ["ML". "REML"]_["LS", "GLMNET"]: `string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "-", model, "_", random_covariate, "-OUTPUT.csv")`
2. Random covariate effects vector (headerless: RANDOM_EFFECTS) written into a comma-separated (.csv) file
	- "GWAlpha": nothing
	- ["ML". "REML"]_["LS", "GLMNET"]: `string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "-", model, "_", random_covariate, "-RANEF-OUTPUT.csv")`
3. Manhattan and QQ plots in .png format
	- `string(join(split(filename_output_csv, ".")[1:(end-1)], '.'), ".png")`
4. Parsing, filtering, and relatedness matrix outputs:
	- Parsed sync data into a .csv file of allele frequency data (headerless: CHROM, POS, ALLELE, FREQ_POOL1, ..., FREQ_POOLn): `string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv")`
	- Filtered sync data into a sync file: filename_sync_filtered = `string(join(split(filename_sync, ".")[1:(end-1)], "."), "_MAF", maf, "_DEPTH", depth, ".sync")`
	- Boolean indices of sync filtering across loci*6 alleles (A:T:C:G:N:DEL) (headerless: IS_INCLUDED): `string(join(split(filename_sync, ".")[1:(end-1)], "."), "_MAF", maf, "_DEPTH", depth, "_IDX_OUT.txt")`
	- Pairwise pool relatedness (square, symmetric and headerless): `string(join(split(filename_sync_filtered, ".")[1:(end-1)], '.'), "_COVARIATE_", random_covariate, ".csv")`

## Examples
```julia
### Input files:
filename_sync = "test/test.sync"
filename_phen_py = "test/test.py"
filename_phen_csv = "test/test.csv"
### Single-threaded execution:
using GWAlpha
OUT_GWAS = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_py)
OUT_REML_GLMNET_RELATEDNESS = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_csv, maf=0.001, depth=10, model="REML_GLMNET", random_covariate="RELATEDNESS", glmnet_alpha=0.50, fpr=0.01, plot=true)
### Multi-threaded execution (parallel execution only applicable to model="GWAlpha"):
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
@everywhere using GWAlpha
OUT_GWAS = GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen_py, plot=true)
```

## Details
The GWAlpha model is defined as **α = W(μₐₗₗₑₗₑ-μₐₗₜₑᵣₙₐₜᵢᵥₑ)/σᵧ**, where:
- **μ** is the mean of the beta distribution, **Beta(θ)** where **θ={θ₁,θ₂}**
- **θ** is estimated via maximum likelihood **L(θ|Q) ∝ πᵢ₌₁₋ₖf(qᵢ|θ)**
- **Q = {q₁,...,qₖ}** is the cumulative sum of allele frequencies across increasing-phenotypic-value-sorted pools where **k** is the number of pools
- **E(allele|θ) = Beta_cdf(yᵢ',θ) - Beta_cdf(yᵢ₋₁',θ)**, where **yᵢ' ∈ Y'**
- **Y'** is the inverse quantile-normalized into phenotype data such that **Y' ∈ [0,1]**
- **W = 2√{E(allele)*(1-E(allele))}** is the penalization for low allele frequency

Empirical p-values were calculated by modelling the additive allelic effects (α) using a normal distribution with mean and variance estimated using maximum likelihood.

The linear mixed model is defined as **y = Xb + Zu + e**, where:
- **X** [n,p] is the centered matrix of allele frequencies
- **Z** [n,n] is the square symmetric matrix of relatedness
- **y** [n,1] is the centered vector of phenotypic values
- no intercept is explicitly fitted but implicitly set at the mean phenotypic value as a consequence of centering **y**
- **u ~ N(0, σ²uI)**
- **e ~ N(0, σ²eI)**
- **y ~ N(0, V)**
- **V = (Z (σ²uI) Z') + (σ²eI)**
- variance components (**σ²e**, **σ²u**) are estimated via maximum likelihood (ML) or restricted maximum likelihood (REML)
- fixed effects (**b**) are estimated via least squares (LS) or elastic-net penalization (GLMNET*; default: α=0.00 which is ridge regression)
- random effects (**u**) are estimated by solving: **(σ²uI) * Z' * V⁻¹ * (y - (X*b))**

GLMNET cross-validation to find the optimum tuning parameter (**λ**) was performed once for the fixed model: **y = Xb + e** to expedite variance components estimation vial ML or REML. The tuning parameter which minimized the mean squared error is selected.

## Citation
Fournier-Level A, Robin C, Balding DJ (2016). [GWAlpha: Genome-Wide estimation of additive effects (Alpha) based on trait quantile distribution from pool-sequencing experiments.](https://doi.org/10.1093/bioinformatics/btw805)
