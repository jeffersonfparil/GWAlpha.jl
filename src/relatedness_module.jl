module relatedness_module

using Statistics
using DelimitedFiles
include("sync_processing_module.jl")

#####################
### Sub-functions ###
#####################
### [Weir and Cockerham, 1984 method](https://www.jstor.org/stable/2408641?seq=1)
function Fst_weir1984_func(sync::Array{Any,2}, pool_sizes::Array{Int64,1}, window_size::Int64)
    CHROM=[]; WINDOW=[]; FST=[]
    for chrom in unique(sync[:,1])
        # test
        # chrom = 1
        subsync = sync[sync[:,1] .== chrom, :]
        nWindows = convert(Int, ceil(maximum(subsync[:,2]) / window_size))
        for window in 1:nWindows
            # window = 1
            # println(chrom)
            # println(window)
            # println("#######################")
            _start_ = window_size*(window-1) + 1
            _end_ = window_size*(window)
            windowsync = subsync[(subsync[:,2] .>= _start_) .& (subsync[:,2] .<= _end_), :]
            if (size(windowsync,1) != 0)
                l = size(windowsync,1)     #number of loci in this window
                a = 6                               #number of alleles which is 6 for the sync formal
                r = size(windowsync,2) - 3 #number of pools or sub-populations or populations in the sync file
                p_3d = Array{Float64}(undef, l, a, r)
                for pop in 1:r
                    for locus in 1:l
                        # #test:
                        # pop=1; locus=1
                        COUNTS = parse.(Int64, split(windowsync[locus, 3+pop], ":"))
                        FREQS = COUNTS ./ sum(COUNTS)
                        p_3d[locus, :, pop] = FREQS
                    end
                end
                # # remove fixed loci
                # idx_fixed_loci = Statistics.var(Statistics.maximum(p_3d, dims=2), dims=3) .> 1e-5
                # p_3d = p_3d[idx_fixed_loci[:],:,:]
                # # redefine the number of loci
                # l = size(p_3d)[1]
                # start computations
                n_bar = sum(pool_sizes) / r
                n_c = ( r*n_bar - ((sum(pool_sizes.^2))/(r*n_bar)) ) / (r-1)
                n_3d = reshape(repeat(pool_sizes, inner=l*a), (l,a,r))
                p_bar_2d = sum(n_3d .* p_3d, dims=3) ./ (r*n_bar)
                s_sq_2d = sum(n_3d .* ((p_3d .- p_bar_2d) .^ 2), dims=3) ./ (n_bar*(r-1))
                # building the HWE-based P(heterozygotes)
                h_3d = 2 .* p_3d .* (1 .- p_3d)
                h_bar_2d = sum(n_3d .* h_3d, dims=3) ./ (r * n_bar)
                # calculate the components of Fst
                a_2d = (n_bar/n_c) * ( (s_sq_2d) - ((1/(n_bar-1)) .* ( (p_bar_2d .* (1 .- p_bar_2d)) .- (((r-1)/r) .* s_sq_2d) .- (h_bar_2d ./ 4) )) )
                b_2d = (n_bar/(n_bar-1)) .* ( (p_bar_2d.* (1 .- p_bar_2d)) .- (((r-1)/r) .* s_sq_2d) .- ((((2*n_bar)-1)/(4*n_bar)) .* h_bar_2d) )
                c_2d = h_bar_2d ./ 2
                # calculate Fst and append to output vectors
                if sum(a_2d) <= 0.0 #for Fst values less than or equal to zero are essentially zero, i.e. Fst = (H_between - H_within) / H_between; --> H_between </= H_within
                    Fst = 0.0
                    append!(FST, Fst)
                else
                    Fst = sum(a_2d) / sum( a_2d .+ b_2d .+ c_2d )
                    append!(FST, Fst)
                end
                push!(CHROM, chrom)
                append!(WINDOW, window)
            end #non-zero window size
        end #window
    end #chrom
    OUT_Fst = try
                hcat((chr=convert(Array{String},CHROM), window=convert(Array{Int64},WINDOW), Fst=convert(Array{Float64},FST))...)
            catch
                hcat((chr=convert(Array{String},string.(CHROM)), window=convert(Array{Int64},WINDOW), Fst=convert(Array{Float64},FST))...)
            end
    return(OUT_Fst)
end

### [Hivert et al, 2018 method](https://www.biorxiv.org/content/biorxiv/early/2018/03/20/282400.full.pdf)
function Fst_hivert2018_func(sync::Array{Any,2}, pool_sizes::Array{Int64,1}, window_size::Int64)
    CHROM=[]; WINDOW=[]; FST=[]
    for chrom in unique(sync[:,1])
        # # test:
        # chrom = 1
        subsync = sync[sync[:,1] .== chrom, :]
        nWindows = convert(Int, ceil(maximum(subsync[:,2]) / window_size))
        for window in 1:nWindows
            # # test:
            # window = 2
            _start_ = window_size*(window-1) + 1
            _end_ = window_size*(window)
            windowsync = subsync[(subsync[:,2] .>= _start_) .& (subsync[:,2] .<= _end_), :]
            if (size(windowsync,1) != 0)
                l = size(windowsync,1)     #number of loci in this window
                a = 6                               #number of alleles which is 6 for the sync formal
                r = size(windowsync,2) - 3 #number of pools or sub-populations or populations in the sync file
                c_3d = Array{Int64}(undef, l, a, r)
                p_3d = Array{Float64}(undef, l, a, r)
                for pop in 1:r
                    for locus in 1:l
                        COUNTS = parse.(Int64, split(windowsync[locus, 3+pop], ":"))
                        FREQS = COUNTS ./ sum(COUNTS)
                        c_3d[locus, :, pop] = COUNTS
                        p_3d[locus, :, pop] = FREQS
                    end
                end
                C1 = sum(c_3d)
                C1i = reshape(sum(c_3d, dims=[1, 2]), r)
                D2 = sum( (C1i .+ (pool_sizes .- 1)) ./ pool_sizes )
                Pik = reshape(mean(p_3d, dims=1), a, r)'
                Pk = reshape(mean(p_3d, dims=[1,3]), a)'
                # ### filter-out alleles with zero frequencies
                # Pik = Pik[:, (Pk .!= 0.0)[:]]
                # Pk = Pk[:, (Pk .!= 0.0)[:]]
                D2a = sum(C1i .* (C1i .+ (pool_sizes .- 1)) ./ pool_sizes) / C1
                C2 = sum(C1i.^2)
                nc = (C1 - (C2/C1)) / (D2 - 1)
                MSP_k = (C1 - D2) .* sum( C1i .* ((Pik .- Pk).^2), dims=1 )
                MSI_k = (D2 - D2a) .* sum( C1i .* Pik .* (1.0 .- Pik), dims=1 )
                Fst = sum( MSP_k .- MSI_k ) / sum( MSP_k .+ ((nc-1) .* MSI_k) )
                if sum( MSP_k .+ ((nc-1) .* MSI_k) ) == 0.0 ### zero denominator indicating ZERO DIFFERENTIATON!
                    Fst = 0.0
                elseif (Fst < 0.0) ### negative Fst interpreted as zero!
                    Fst = 0.0
                end
                append!(FST, Fst)
                push!(CHROM, chrom)
                append!(WINDOW, window)
            end #non-zero window size
        end #window
    end #chrom
    OUT_Fst = try
                hcat((chr=convert(Array{String},CHROM), window=convert(Array{Int64},WINDOW), Fst=convert(Array{Float64},FST))...)
            catch
                hcat((chr=convert(Array{String},string.(CHROM)), window=convert(Array{Int64},WINDOW), Fst=convert(Array{Float64},FST))...)
            end
    return(OUT_Fst)
end

######################
### Main functions ###
######################
"""
# ___________________________________________________________________
# Fixation index (Fst) estimation from Pool-seq data across all pools

`Fst_across_pools(;sync_fname::String, window_size::Int64, pool_sizes::Array{Int64,1}, METHOD::String)`

Compute **across pools Fst** using [Weir and Cockerham, 1984 method](https://www.jstor.org/stable/2408641?seq=1)
or using [Hivert et al, 2018 method](https://www.biorxiv.org/content/biorxiv/early/2018/03/20/282400.full.pdf)

# Input
1. *sync_fname*: synchronized pileup filename
2. *window_size*: size of non-overlapping sliding window in base-pairs
3. *pool-sizes*: array of pool sizes corresponding to each genotype column in the sync file
4. *METHOD*: Fst calculation method to use:
- "WeirCock"
- "Hivert"

# Output
1. Fst estimates per window with the filename: `string(join(split(sync_fname, ".")[1:(end-1)], '.'), "_window_", window_size, "bp_Fst_data.csv")`
2. Mean Fst estimate with the filename: `string(join(split(sync_fname, ".")[1:(end-1)], '.'), "_window_", window_size, "bp_Fst_sumstats.csv")`

# Examples
```
relatedness_module.Fst_across_pools(sync_fname="test/test.sync", window_size=100000, pool_sizes=[20,20,20,20,20], METHOD="WeirCock")
relatedness_module.Fst_across_pools(sync_fname="test/test.sync", window_size=100000, pool_sizes=[20,20,20,20,20], METHOD="Hivert")
```
"""
function Fst_across_pools(;sync_fname::String, window_size::Int64, pool_sizes::Array{Int64,1}, METHOD::String)
    println("#################################################################################################")
    println("Fst estimation across pools:")
    ### load the sync allele frequency data
    sync = DelimitedFiles.readdlm(sync_fname)
    ### input parameters
    nPools = size(sync,2)  - 3
    nLoci = size(sync,1)
    println(string("Synchronized pileup file: ", sync_fname))
    println(string("Non-overlapping window size: ", window_size, " bp"))
    println(string("Number of pools: ", nPools))
    for i in 1:length(pool_sizes)
        println(string("Pool ", i, " size : ", pool_sizes[i]))
    end
    println(string("Number of loci: ", nLoci))
    ### Across pools Fst estimation
    if METHOD == "WeirCock"
        OUT_Fst = Fst_weir1984_func(sync, pool_sizes, window_size)
    elseif METHOD == "Hivert"
        OUT_Fst = Fst_hivert2018_func(sync, pool_sizes, window_size)
    else
        println(string("AwWwwW! SowWwWyyyy! We have not implemented the Pool-seq Fst calculation method: ", METHOD, " yet."))
        println("¯\\_(๑❛ᴗ❛๑)_/¯ ʚ(´◡`)ɞ")
    end
    println("Writing the output into the files:")
    Fst_data_out_fname = string(join(split(sync_fname, ".")[1:(end-1)], '.'), "_window_", window_size, "bp_Fst_data.csv")
    Fst_sumstats_out_fname = string(join(split(sync_fname, ".")[1:(end-1)], '.'), "_window_", window_size, "bp_Fst_sumstats.csv")
    println(string("(1/2) ", Fst_data_out_fname, " for the Fst per chromosome per locus."))
    println(string("(2/2) ", Fst_sumstats_out_fname, " for the summary statistics of Fst across chromosomes and loci (i.e. mean, min, max and var)."))
    println("#################################################################################################")
    ### performing a Shapiro-Wilk's test for normality shows that the Fst across windows is normally distributed!
    SUMSTATS_Fst = (mean_Fst=Statistics.mean(OUT_Fst[.!isnan.(OUT_Fst[:,3]),3]), min_Fst=Statistics.minimum(OUT_Fst[.!isnan.(OUT_Fst[:,3]),3]), max_Fst=Statistics.maximum(OUT_Fst[.!isnan.(OUT_Fst[:,3]),3]), var_Fst=Statistics.var(OUT_Fst[.!isnan.(OUT_Fst[:,3]),3]))
    ### write-out
    writedlm(Fst_data_out_fname, reshape(["CHROM", "POS", "FST"],1,3), ',') ### write header
    io = open(Fst_data_out_fname, "a"); writedlm(io, OUT_Fst, ','); close(io) ### append the column-binded GWAlpha output
    writedlm(Fst_sumstats_out_fname, reshape(vcat(string.(keys(SUMSTATS_Fst))...), 1, length(SUMSTATS_Fst)), ',') ### write header
    io = open(Fst_sumstats_out_fname, "a"); writedlm(io, hcat(SUMSTATS_Fst...), ','); close(io) ### append the column-binded GWAlpha output
    return(0)
end

"""
# _______________________________________________________________
# Pairwise estimation of fixation indices (Fst) from Pool-seq data

`Fst_pairwise(;sync_fname::String, window_size::Int64, pool_sizes::Array{Int64,1}, METHOD::String)`

Compute **pairwise Fst** using [Weir and Cockerham, 1984 method](https://www.jstor.org/stable/2408641?seq=1)
or using [Hivert et al, 2018 method](https://www.biorxiv.org/content/biorxiv/early/2018/03/20/282400.full.pdf)

# Input
1. *sync_fname*: synchronized pileup filename
2. *window_size*: size of non-overlapping sliding window in base-pairs
3. *pool-sizes*: array of pool sizes corresponding to each genotype column in the sync file
4. *METHOD*: Fst calculation method to use:
    - "WeirCock"
    - "Hivert"

# Output
Pairwise Fst estimates with the filename: `string(join(split(sync_fname, ".")[1:(end-1)], '.'), "_COVARIATE_FST.csv")`

# Examples
```
relatedness_module.Fst_pairwise(sync_fname="test/test.sync", window_size=100000, pool_sizes=[20,20,20,20,20], METHOD="WeirCock")
relatedness_module.Fst_pairwise(sync_fname="test/test.sync", window_size=100000, pool_sizes=[20,20,20,20,20], METHOD="Hivert")
```
"""
function Fst_pairwise(;sync_fname::String, window_size::Int64, pool_sizes::Array{Int64,1}, METHOD::String=["WeirCock", "Hivert"][1])
    println("#################################################################################################")
    println("Pairwise Fst estimation:")
    ### load the sync allele frequency data
    sync = DelimitedFiles.readdlm(sync_fname)
    ### input parameters
    nPools = size(sync,2)  - 3
    nLoci = size(sync,1)
    println(string("Synchronized pileup file: ", sync_fname))
    println(string("Non-overlapping window size: ", window_size, " bp"))
    println(string("Number of pools: ", nPools))
    for i in 1:length(pool_sizes)
        println(string("Pool ", i, " size : ", pool_sizes[i]))
    end
    println(string("Number of loci: ", nLoci))
    ### Pairwise Fst estimation
    nPairs = convert(Int64, (nPools * (nPools - 1)) / 2)
    PAIRWISE_Fst = zeros(nPools, nPools)
    for i in 1:(nPools-1)
        for j in (i+1):nPools
            sub = sync[:, [1,2,3,(3+i), (3+j)]]
            if METHOD == "WeirCock"
                Fst = Fst_weir1984_func(sub, pool_sizes[[i,j]], window_size)[:,3]
                PAIRWISE_Fst[i, j] = mean(Fst[.!isnan.(Fst)])
            elseif METHOD == "Hivert"
                Fst = Fst_hivert2018_func(sub, pool_sizes[[i,j]], window_size)[:,3]
                PAIRWISE_Fst[i, j] = mean(Fst[.!isnan.(Fst)])
            else
                println(string("AwWwwW! SowWwWyyyy! We have not implemented the Pool-seq Fst calculation method: ", model, " yet."))
                println("¯\\_(๑❛ᴗ❛๑)_/¯ ʚ(´◡`)ɞ")
            end
            PAIRWISE_Fst[j, i] = PAIRWISE_Fst[i, j]
        end
    end
    Fst_data_out_fname = string(join(split(sync_fname, ".")[1:(end-1)], '.'), "_COVARIATE_FST.csv")
    DelimitedFiles.writedlm(Fst_data_out_fname, PAIRWISE_Fst, ',')
    println("#################################################################################################")
    return(0)
end

"""
# _______________________________
# Standardized relatedness matrix

`standardized_relatedness(sync_fname::String)`

Compute a simple standardized relatedness matrix XX'/p, where X is the allele frequency matrix (Pool-seq data) and p is the number of columns of X

# Input
1. *sync_fname*: synchronized pileup filename

# Output
Standardized relatedness matrix with the filename: `string(join(split(sync_fname, ".")[1:(end-1)], '.'), "_COVARIATE_RELATEDNESS.csv")`

# Examples
```
relatedness_module.standardized_relatedness("test/test.sync")
```
"""
function standardized_relatedness(sync_fname::String)
    println("#################################################################################################")
    println("Relatedness estimation:")
    ### check if the parsed sync file i.e. the allele frequency csv file exists
    ### if it does then use it, otherwise parse the sync file
    filename_parsed_sync = string(join(split(sync_fname, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv")
	if isfile(filename_parsed_sync)
		# rm(filename_out)
		println(string("Using the existing parsed sync file in csv format: ", filename_parsed_sync))
    else
        sync_processing_module.sync_parse(sync_fname)
	end
    X = convert(Array{Float64,2},DelimitedFiles.readdlm(filename_parsed_sync, ',')[:,4:end])'
    K = X*X' ./ size(X,2)
    relatedness_fname = string(join(split(sync_fname, ".")[1:(end-1)], '.'), "_COVARIATE_RELATEDNESS.csv")
    DelimitedFiles.writedlm(relatedness_fname, K, ',')
    println("#################################################################################################")
    return(0)
end


end ### end relatedness_module
