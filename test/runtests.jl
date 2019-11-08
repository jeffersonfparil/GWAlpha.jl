using Test
using DelimitedFiles
using GWAlpha

# DIR = replace(dirname(pathof(GWAlpha)), "/src" => "/test")
# geno_sync_fname = string(DIR, "/LOLIUM_1rep_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_g500_p16_POOLS_GENO.sync")
# pheno_py_fname = string(DIR, "/LOLIUM_1rep_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_g500_p16_POOLS_PHENO.py")
# pheno_csv_fname = string(DIR, "/LOLIUM_1rep_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_g500_p16_POOLS_PHENO.csv")
# covariate_fname = string(DIR, "/LOLIUM_1rep_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_g500_p16_POOLS_GENO_COVARIATE_FST.csv")

geno_sync_fname = "LOLIUM_1rep_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_g500_p16_POOLS_GENO.sync"
pheno_py_fname = "LOLIUM_1rep_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_g500_p16_POOLS_PHENO.py"
pheno_csv_fname = "LOLIUM_1rep_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_g500_p16_POOLS_PHENO.csv"
covariate_fname = "LOLIUM_1rep_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_g500_p16_POOLS_GENO_COVARIATE_FST.csv"

function runGWAlpha(filename_sync::String, filename_phen::String, MAF::Float64, DEPTH::Int64; MODEL="FIXED_GWAlpha"::String, COVARIATE=nothing)
    GWAlpha.PoolGPAS(filename_sync, filename_phen, MAF, DEPTH, MODEL=MODEL, COVARIATE=COVARIATE)
    return(0)
end

### GWAlpha for GWAS
Test.@test runGWAlpha(geno_sync_fname,
                      pheno_py_fname,
                      0.01,
                      10,
                      MODEL="FIXED_GWAlpha",
                      COVARIATE=nothing) == 0
### GWAlpha for GP
Test.@test runGWAlpha(geno_sync_fname,
                      pheno_csv_fname,
                      0.01,
                      10,
                      MODEL="FIXED_RR",
                      COVARIATE=DelimitedFiles.readdlm(covariate_fname, ',')) == 0
