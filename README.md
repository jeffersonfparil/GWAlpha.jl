# GWAlpha.jl

<a href="https://adaptive-evolution.biosciences.unimelb.edu.au/"><img src="misc/Adaptive Evolution Logo mod.png" width="150"></a>
[![Build Status](https://travis-ci.com/jeffersonfparil/GWAlpha.jl.svg?branch=master)](https://travis-ci.com/jeffersonfparil/GWAlpha.jl)
<a href="https://github.com/jeffersonfparil/GWAlpha.jl/wiki" target="_blank"><img src="https://img.shields.io/badge/docs-latest-blue.svg" alt="Latest documentation"></a>

<!--- [![CircleCI](https://circleci.com/gh/jeffersonfparil/GWAlpha.svg?style=shield)](https://circleci.com/gh/jeffersonfparil/GWAlpha) --->

A [Julia](https://julialang.org/downloads/) package for genome-wide association (GWAS) and genomic prediction (GP) of quantitative traits from pool sequencing (Pool-seq) data.

This repository include a suite of scripts for genomic landscape simulation using [quantinemo2](https://github.com/jgx65/quantinemo) to test individual-based and Pool-seq-based GWAS and GP algorithms and models, which makes use of [npstat](https://github.com/lucaferretti/npstat) for extracting population genetics summary statistics as covariates for some models.

A mirror repository is found in [gitlab](https://gitlab.com/jeffersonfparil/genomic_prediction).

## Installation
```
using Pkg
Pkg.add(PackageSpec(url="https://github.com/jeffersonfparil/GWAlpha.jl.git", rev="master"))
Pkg.build()
```

## Testing
```
using Pkg
Pkg.update("GWAlpha")
Pkg.test("GWAlpha")
```

## Contents

- original GWAlpha implemented in python in the [legacy directory](https://github.com/jeffersonfparil/GWAlpha.jl/tree/master/legacy)
- Julia, shell, and R scripts are located in the [src directory](https://github.com/jeffersonfparil/GWAlpha.jl/tree/master/src)
- testing scripts are found in the [test directory](https://github.com/jeffersonfparil/GWAlpha.jl/tree/master/test)

## Citations

Fournier-Level A, Robin C, Balding DJ (2016). GWAlpha: Genome-Wide estimation of additive effects (Alpha) based on trait quantile distribution from pool-sequencing experiments. Submitted to Bioinformatics
