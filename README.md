# GWAlpha.jl

|                                                          **Lab Website**                                                          |                                                            **Build Status**                                                             |                                                                             **Documentation**                                                                             |
|:---------------------------------------------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| <a href="https://adaptive-evolution.biosciences.unimelb.edu.au/"><img src="misc/Adaptive Evolution Logo mod.png" width="150"></a> | [![Build Status](https://travis-ci.com/jeffersonfparil/GWAlpha.jl.svg?branch=master)](https://travis-ci.com/jeffersonfparil/GWAlpha.jl) | <a href="https://github.com/jeffersonfparil/GWAlpha.jl/wiki" target="_blank"><img src="https://img.shields.io/badge/docs-latest-blue.svg" alt="Latest documentation"></a> |

<!--- [![CircleCI](https://circleci.com/gh/jeffersonfparil/GWAlpha.svg?style=shield)](https://circleci.com/gh/jeffersonfparil/GWAlpha) --->

A [Julia](https://julialang.org/downloads/) package for genome-wide association (GWAS) and genomic prediction (GP) of quantitative traits from pool sequencing (Pool-seq) data.

This repository include a suite of scripts for genomic landscape simulation using [quantinemo2](https://github.com/jgx65/quantinemo) to test individual-based and Pool-seq-based GWAS and GP algorithms and models, which makes use of [npstat](https://github.com/lucaferretti/npstat) for extracting population genetics summary statistics as covariates for some models.

A mirror repository is found in [gitlab](https://gitlab.com/jeffersonfparil/genomic_prediction).

## Installation
Install dependencies:
```shell
sudo apt install at-spi2-core libgtk-3-dev xauth xvfb
pip3 install --upgrade pip==19.3.1 --user
pip install setuptools --user
pip install numpy --user
pip install matplotlib --user
```
Installation in Julia:
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/jeffersonfparil/GWAlpha.jl.git", rev="master"))
Pkg.build()
```

## Testing
```julia
using Pkg
Pkg.update("GWAlpha")
Pkg.test("GWAlpha")
```

## Inputs

1. [synchronized pileup filename](https://sourceforge.net/p/popoolation2/wiki/Manual/)
2. phenotype data filename
- **.py** extension for iterative maximum likelihood estimation i.e. `MODEL="FIXED_GWAlpha"`, e.g.:
```julia
Pheno_name='Phenotype Name';
sig=0.06724693662723039;	# standard deviation
MIN=0.0;			# minimum phenotype value
MAX=0.424591738712776;		# maximum phenotype value
perc=[0.2,0.4,0.6,0.8];		# cummulative pool sizes percentiles excluding the last pool
q=[0.16,0.20,0.23,0.27,0.42];	# phenotype values corresponding to each percentile
```
- **.csv** extension for comma-separated headerless poolsizes and corresponding mean phenotype values, e.g.:
```julia
200.0,0.11988952929875112
200.0,0.18030259365994225
200.0,0.21548030739673382
200.0,0.24966378482228616
200.0,0.31328530259365983
```
3. minimum allele frequency threshold
4. minimum sequencing depth threshold
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

## Outputs

1. DataFrames.DataFrame of additive allele effects with the corresponding identification (CHROM, POS, ALLELE, FREQ)
2. Array of covariate effects
3. Additive allele effects csv file:
```julia
string(dir, replace(filename, ".py" => string("-", MODEL, "_Alphas.csv"))) ### or...
string(dir, replace(filename, ".csv" => string("-", MODEL, "_Alphas.csv")))
```
4. Manhattan plot png format:
```julia
string(dir, replace(filename, ".py" => string("-", MODEL, "_Manhattan.png"))) ### or...
string(dir, replace(filename, ".csv" => string("-", MODEL, "_Manhattan.png")))
```

## More details

Open julia, load the GWAlpha library,
```julia
using GWAlpha
?GWAlpha.PoolGPAS
```

## Contents

- original GWAlpha implemented in python in the [legacy directory](https://github.com/jeffersonfparil/GWAlpha.jl/tree/master/legacy)
- Julia, shell, and R scripts are located in the [src directory](https://github.com/jeffersonfparil/GWAlpha.jl/tree/master/src)
- testing scripts are found in the [test directory](https://github.com/jeffersonfparil/GWAlpha.jl/tree/master/test)

## Citations

Fournier-Level A, Robin C, Balding DJ (2016). [GWAlpha: Genome-Wide estimation of additive effects (Alpha) based on trait quantile distribution from pool-sequencing experiments.](https://doi.org/10.1093/bioinformatics/btw805)
