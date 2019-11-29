##############################
### ANALYSIS PER LANDSCAPE ### MODULE 2 of 3:
############################## ABC simulations via parallel processing

### simulate summary stats across this parameter space for NREP replications
### In Julia
# julia -p12
### inputs
DIR = ARGS[1]   ### GPASIM output directory
nQTL = parse(Int64, ARGS[2]) ### number of QTL
nREP = parse(Int64, ARGS[3]) ### number of replications
# ### TEST:
# DIR = "/data/Lolium/Quantitative_Genetics/LOLSIM_2019_TEST/LOLSIM_1rep_5qtl_0.001mr_0.25fgs_0.00bgs_0grad"
# nQTL = 5
# nREP = 10
### load libraries and modules
using Distributed
Distributed.addprocs(length(Sys.cpu_info()))
@everywhere using Distributed
@everywhere using DataFrames
@everywhere using CSV
@everywhere using DelimitedFiles
@everywhere using Statistics
@everywhere using ProgressMeter
@everywhere using SharedArrays
using ProgressMeter
### set working directory
cd(DIR)
### load the streamlined GPAS cross-validation output
@everywhere ACROSS_POOL = CSV.read("ACROSS_POOL.csv")
@everywhere ACROSS_INDI = CSV.read("ACROSS_INDI.csv")
@everywhere WITHIN_POOL = CSV.read("WITHIN_POOL.csv")
@everywhere WITHIN_INDI = CSV.read("WITHIN_INDI.csv")
### a very inelegant fix to putting nQTL on all threads
@everywhere nQTL = convert(Int64, DelimitedFiles.readdlm("nQTL.temp", ',')[1,1])
### generate the parameter space
@everywhere across_pool = []
@everywhere across_indi = []
@everywhere within_pool = []
@everywhere within_indi = []
epsilon = 0.01
pb = ProgressMeter.Progress(length(collect(0.0:epsilon:1.0))^4, dt=1, desc="Generating parameter space: ",  barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow); counter=[1]
for i in collect(0.0:epsilon:1.0)
    for j in collect(0.0:epsilon:1.0)
        for k in collect(0.0:epsilon:1.0)
            for l in collect(0.0:epsilon:1.0)
                if (i+j+k+l == 1.00)
                    push!(across_pool, i)
                    push!(across_indi, j)
                    push!(within_pool, k)
                    push!(within_indi, l)
                end
                ProgressMeter.update!(pb, counter[1]); counter[1] = counter[1] + 1
            end
        end
    end
end
@everywhere PAR = DataFrames.DataFrame(ACROSS_POOL=across_pool, ACROSS_INDI=across_indi, WITHIN_POOL=within_pool, WITHIN_INDI=within_indi)
### define the sampling_summstats_function
@everywhere function SAMPLE_ABC(;par::Array{Float64,1}, nLib::Int64=10, nQTL::Int64=nQTL, ACROSS_POOL::DataFrames.DataFrame=ACROSS_POOL, ACROSS_INDI::DataFrames.DataFrame=ACROSS_INDI, WITHIN_POOL::DataFrames.DataFrame=WITHIN_POOL, WITHIN_INDI::DataFrames.DataFrame=WITHIN_INDI)
  ### test
  # par = [0.25, 0.25, 0.25, 0.25]
  # nLib = 10
  # nQTL = 5
  ### random sampling and summary statistics i.e CORRELATION, LOG10_RMSD, TRUE_POSITIVE_RATE, and FALSE_POSITIVE_RATE
  DATA_LEVELS = [ACROSS_POOL, ACROSS_INDI, WITHIN_POOL, WITHIN_INDI]
  CORRELATION = []
  LOG10_RMSD = []
  TRUE_POSITIVE_ID = []
  FALSE_POSITIVE_ID = []
  for i in 1:length(DATA_LEVELS)
    # i = 1
    LEVEL = DATA_LEVELS[i]
    n = convert(Int64, round(par[i]*nLib))
    if (n > 0)
      idx = rand(collect(1:size(LEVEL,1)), n)
      # println(idx)
      push!(CORRELATION, mean(dropmissing(LEVEL[idx,:]).CORRELATION))
      push!(LOG10_RMSD, mean(dropmissing(LEVEL[idx,:]).LOG10_RMSD))
      push!(TRUE_POSITIVE_ID, split(join(convert(Array{String}, dropmissing(LEVEL[idx,:]).TRUE_POSITIVE_ID), ";"), ";"))
      push!(FALSE_POSITIVE_ID, split(join(convert(Array{String}, dropmissing(LEVEL[idx,:]).FALSE_POSITIVE_ID), ";"), ";"))
    end
  end
  if (sum(isnan.(CORRELATION))>0) | (sum(ismissing.(CORRELATION))>0)
    SAMPLE_ABC(par=par, nLib=nLib, nQTL=nQTL)
  else
    OUT = [mean(CORRELATION[.!isnan.(CORRELATION)]),
           mean(LOG10_RMSD[.!isnan.(LOG10_RMSD)]),
           length(unique(TRUE_POSITIVE_ID))/nQTL,
           length(unique(FALSE_POSITIVE_ID))/(length(unique(TRUE_POSITIVE_ID))+length(unique(FALSE_POSITIVE_ID)))]
    return(OUT)
  end
end
### parallel computation using @distributed and prefixed with
### @sync to wait for the whole parallel computation to finish before runing the next line
OUT = SharedArrays.SharedArray{Float64,2}(size(PAR,1)*nREP, 4)
@time x = @sync @distributed for i = 1:size(OUT,1)
  idx_PAR = convert(Int64, i - (floor((i-1)/size(PAR,1))*size(PAR,1)))
  OUT[i,:] = SAMPLE_ABC(par=convert(Array{Float64,1},PAR[idx_PAR,:]), nLib=10)
end
### merge the parameters and the corresponding summary statistics from the simulated (sampled) data
MERGED = DataFrames.DataFrame(REP=repeat(collect(1:nREP), inner=size(PAR,1)),
                              ACROSS_POOL=repeat(PAR.ACROSS_POOL, outer=nREP),
                              ACROSS_INDI=repeat(PAR.ACROSS_INDI, outer=nREP),
                              WITHIN_POOL=repeat(PAR.WITHIN_POOL, outer=nREP),
                              WITHIN_INDI=repeat(PAR.WITHIN_INDI, outer=nREP),
                              CORRELATION=OUT[:,1],
                              LOG10_RMSD=OUT[:,2],
                              TRUE_POSITIVE_RATE=OUT[:,3],
                              FALSE_POSITIVE_RATE=OUT[:,4])
CSV.write(string("ABC_SAMPLING_OUTPUT_", nREP, "REPS.csv"), MERGED)
