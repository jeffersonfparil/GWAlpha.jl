##########################################################
###                                                    ### @Weir1984
### Fst across all populations and all loci per window ### Estimation F-statistics for the analysis of population structure. 1984. Evolution. Vol 38. No 6)
###       using Pool-seq data (i.e. sync format)       ### (Weir and Cockerham's method;
##########################################################
# sample execution
# julia ${GEN_PRED_SRC_DIR}/GPASim_03_Fst.jl \
# QUANTI_g140_POOLS_GENO.sync \
# 1000000 \
# 20,20,20,20,20 \
# WeirCock

##############
### inputs ###
##############
sync_fname = ARGS[1]
window_size = parse(Int64, ARGS[2])
pool_sizes = parse.(Int64, split(ARGS[3], ",")) #poolsizes delimited by comma e.g. 10,20,30,20,10
METHOD = ARGS[4] ### WeirCock or Hivert
# ################
# ### test
# window_size = parse(Int64, "1000000") #non-overlapping 10 Mb windows
# # sync_fname = "test_1kloci_ALLPOP_GENO.sync"
# sync_fname = "Lolium_2019_60pop.sync"
# # # sync_fname = "p09_p13_p12_p08_MERGED_POOLS_GENO.sync"
# # # pool_sizes = [46,46,46,46]
# # pool_sizes = repeat([20], inner=16)
# pool_sizes = repeat([42], inner=60)
# # sync_fname = "test_1kloci_g1000_p01_POOLS_GENO.sync"
# # pool_sizes = repeat([20], inner=5)
# # sync_fname = "ZZZ_SYNC_DIVERGED.sync"
# # sync_fname = "ZZZ_SYNC_MERGED.sync"
# # pool_sizes = repeat([20], inner=2)
# METHOD = "WeirCock"
# ################

######################
### load libraries ###
######################
using CSV
using DataFrames
using Statistics
using DelimitedFiles
JULIA_SCRIPT_HOME = @__DIR__
# JULIA_SCRIPT_HOME = "/data/Lolium/Softwares/genomic_prediction/src"
push!(LOAD_PATH, JULIA_SCRIPT_HOME)
using poolFST_module

poolFST_module.Fst_across_pools(sync_fname=sync_fname, window_size=window_size, pool_sizes=pool_sizes, METHOD=METHOD)
poolFST_module.Fst_pairwise(sync_fname=sync_fname, window_size=window_size, pool_sizes=pool_sizes, METHOD=METHOD)
