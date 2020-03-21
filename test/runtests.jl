using Test
using DelimitedFiles
using GWAlpha

geno_sync_fname = "test.sync"
pheno_py_fname = "test.py"
pheno_csv_fname = "test.csv"

function runGWAlpha(filename_sync::String, filename_phen::String, MAF::Float64, DEPTH::Int64, MODEL::String, COVARIATE::Any, FPR::Float64)
    GWAlpha.PoolGPAS(filename_sync, filename_phen, MAF=MAF, DEPTH=DEPTH, MODEL=MODEL, COVARIATE=COVARIATE, FPR=FPR)
    return(0)
end

### GWAlpha for GWAS
Test.@test runGWAlpha(geno_sync_fname,
                      pheno_py_fname,
                      0.01,
                      10,
                      "FIXED_GWAlpha",
                      nothing,
                      0.01) == 0
### GWAlpha for GP
pool_sizes = convert(Array{Int}, DelimitedFiles.readdlm(pheno_csv_fname, ',')[:,1])
Test.@test GWAlpha.poolFST_module.Fst_pairwise(sync_fname=geno_sync_fname, window_size=100000, pool_sizes=pool_sizes, METHOD="Hivert") == 0
covariate_fname = "test_COVARIATE_FST.csv"
Test.@test runGWAlpha(geno_sync_fname,
                      pheno_csv_fname,
                      0.01,
                      10,
                      "FIXED_RR",
                      DelimitedFiles.readdlm(covariate_fname, ','),
                      0.01) == 0
