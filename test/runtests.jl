using Test
using GWAlpha

geno_sync_fname = "test.sync"
pheno_py_fname = "test.py"
pheno_csv_fname = "test.csv"

function runGWAlpha(filename_sync::String, filename_phen::String, MAF::Float64, DEPTH::Int64, MODEL::String)
    GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen, maf=MAF, depth=DEPTH, model=MODEL, plot=true)
    return(0)
end

### GWAlpha for GWAS (non-parallel)
Test.@test runGWAlpha(geno_sync_fname,
                      pheno_py_fname,
                      0.01,
                      10,
                      "GWAlpha") == 0

### GWAlpha for GP
Test.@test runGWAlpha(geno_sync_fname,
                      pheno_csv_fname,
                      0.01,
                      10,
                      "REML_GLMNET") == 0

### GWAlpha for GWAS (parallel)
using Distributed
Distributed.addprocs(length(Sys.cpu_info()))
@everywhere using GWAlpha
function runGWAlpha(filename_sync::String, filename_phen::String, MAF::Float64, DEPTH::Int64, MODEL::String)
    GWAlpha.PoolGPAS(filename_sync=filename_sync, filename_phen=filename_phen, maf=MAF, depth=DEPTH, model=MODEL, plot=true)
    return(0)
end
Test.@test runGWAlpha(geno_sync_fname,
                      pheno_py_fname,
                      0.01,
                      10,
                      "GWAlpha") == 0
