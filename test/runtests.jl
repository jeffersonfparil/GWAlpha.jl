using Test
using GWAlpha

function runGWAlpha(filename_sync::String, filename_phen::String, MAF::Float64, DEPTH::Int64; MODEL="FIXED_GWAlpha", COVARIATE=nothing)
    GWAlpha.PoolGPAS(filename_sync::String, filename_phen::String, MAF::Float64, DEPTH::Int64; MODEL="FIXED_GWAlpha", COVARIATE=nothing)
    return(0)
end

### GWAlpha for GWAS
Test.@test runGWAlpha("LOLIUM_1rep_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_g500_p16_POOLS_GENO.sync",
                      "LOLIUM_1rep_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_g500_p16_POOLS_PHENO.py",
                      0.01,
                      10,
                      MODEL="FIXED_GWAlpha",
                      COVARIATE=nothing) == 0
### GWAlpha for GP
Test.@test runGWAlpha("LOLIUM_1rep_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_g500_p16_POOLS_GENO.sync",
                      "LOLIUM_1rep_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_g500_p16_POOLS_PHENO.csv",
                      0.01,
                      10,
                      MODEL="FIXED_RR",
                      COVARIATE="LOLIUM_1rep_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_g500_p16_POOLS_GENO_COVARIATE_FST.csv") == 0
