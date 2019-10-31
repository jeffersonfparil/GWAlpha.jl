### INPUTS
SRC_DIR = ARGS[1]
filename_sync = ARGS[2]
filename_phen_py = ARGS[3]
filename_phen_csv = ARGS[4]
filename_covariate_csv = ARGS[5]
MAF = parse(Float64, ARGS[6])
DEPTH = parse(Int64, ARGS[7])
MODEL = ARGS[8]
# ### TESTS
# cd("/data/Lolium/Quantitative_Genetics/GWAS_GP_2018_Inverleigh_Urana")
# SRC_DIR = "/data/Lolium/Softwares/genomic_prediction/src"
# filename_sync = "VCF/UG_MAPQ20_BASQ20_MAF0.001_DEPTH10.sync"
# filename_phen_py = "UG_pheno.py"
# filename_phen_csv = "UG_pheno.csv"
# filename_covariate_csv = 'VCF/UG_MAPQ20_BASQ20_MAF0.001_DEPTH10_COVARIATE_FST.csv'
# MAF = 0.001
# DEPTH = 10
# MODEL = "FIXED_GWAlpha" ### FIXED_GWAlpha, FIXED_LS, FIXED_RR, FIXED_GLMNET, FIXED_LASSO, MIXED_RR, MIXED_GLMNET, and MIXED_LASSO
### LOAD LIBRARIES
using DelimitedFiles
using DataFrames
using CSV
push!(LOAD_PATH, SRC_DIR)
using PoolGPAS_module
### LOAD COVARIATE
if MODEL == "FIXED_GWAlpha"
COVARIATE = nothing
else
COVARIATE = DelimitedFiles.readdlm(filename_covariate_csv, ',')
  end
### EXECUTE
println("#######################################################")
println(MODEL)
println("#######################################################")
if MODEL == "FIXED_GWAlpha"
  OUT, COVAR_EFF = PoolGPAS_module.PoolGPAS(filename_sync, filename_phen_py, MAF, DEPTH, MODEL=MODEL, COVARIATE=COVARIATE)
else
  OUT, COVAR_EFF = PoolGPAS_module.PoolGPAS(filename_sync, filename_phen_csv, MAF, DEPTH, MODEL=MODEL, COVARIATE=COVARIATE)
end
### OUTPUTS:
### (1) Allelic effects file (*-_Alphas.csv)
### (2) Manhattan plot (*-_Manhattan.png)
