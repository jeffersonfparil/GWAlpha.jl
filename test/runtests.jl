using Test
using GWAlpha

function runGWAlpha(filename_sync::String, filename_phen::String, MAF::Float64, DEPTH::Int64; MODEL="FIXED_GWAlpha", COVARIATE=nothing)
    PoolGPAS(filename_sync::String, filename_phen::String, MAF::Float64, DEPTH::Int64; MODEL="FIXED_GWAlpha", COVARIATE=nothing)
    return(0)
end

### GWAlpha for GWAS
Test.@test runGWAlpha("UG_MAPQ20_BASQ20_MAF0.001_DEPTH10.sync",
                      "UG_pheno.py",
                      0.01,
                      10,
                      MODEL="FIXED_GWAlpha",
                      COVARIATE=nothing) == 0
### GWAlpha for GP
Test.@test runGWAlpha("UG_MAPQ20_BASQ20_MAF0.001_DEPTH10.sync",
                      "UG_pheno.csv",
                      0.01,
                      10,
                      MODEL="FIXED_RR",
                      COVARIATE="UG_MAPQ20_BASQ20_MAF0.001_DEPTH10_COVARIATE_FST.csv") == 0
