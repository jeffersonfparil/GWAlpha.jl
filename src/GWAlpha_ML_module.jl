module GWAlpha_ML_module

using Distributed
using SharedArrays
using DelimitedFiles
using Distributions
using LinearAlgebra
using Optim
using ProgressMeter
include("significance_testing_module.jl")

"""
# _____________________________________________________________________________________________________________________________________________
# Negative log-likelihood of the Beta distribution parameters (shape1 and shape2) given the distribution of the allele frequencies across pools

`neg_log_likelihood_cdfbeta(beta::Array{Float64,1}, data_A::Array{Float64,1}, data_B::Array{Float64,1})`

# Inputs
1. beta [Array{Float64,1}]: 4-element vector of the beta distribution shape parameters, 2 for each of the data distributions
2. data_A [Array{Float64,1}; xi ∉ 0]: allele frequency distribution across pools
3. data_B [Array{Float64,1}; xi ∉ 0]: additive inverse of the allele frequency distribution across pools

# Output
1. Negative log-likelihood of the Beta distribution parameters (shape1 and shape2) given the distribution of the allele frequencies across pools

# Examples
```
beta = [1.0, 1.0, 0.5, 1.0]
f = rand(5); bins = repeat([0.2], inner=5); p = 0.5
data_A = cumsum( (f .* bins ./ p) )
data_B = cumsum( ((1 .- f) .* bins ./ (1 .- p)) )
GWAlpha_ML_module.neg_log_likelihood_cdfbeta(beta, data_A, data_B)
```
"""
function neg_log_likelihood_cdfbeta(beta::Array{Float64,1}, data_A::Array{Float64,1}, data_B::Array{Float64,1})
	-sum(
		 log.(10,
		 	 (Distributions.cdf.(Distributions.Beta(beta[1], beta[2]), data_A)) .-
		 	 (Distributions.cdf.(Distributions.Beta(beta[1], beta[2]), append!(zeros(1), data_A[1:(length(data_A)-1)])))
		     )
		 )     -sum(
		 log.(10,
		 	 (Distributions.cdf.(Distributions.Beta(beta[3], beta[4]), data_B)) .-
		 	 (Distributions.cdf.(Distributions.Beta(beta[3], beta[4]), append!(zeros(1), data_B[1:(length(data_B)-1)])))
		     )
		 )
end

"""
# _____________________________________________________________________________________________________________________________________________
# Iterating function per locus for genome-wide additive allelic effects (alpha) association using Pool-seq data

`GWAlpha_ML_iterator(COUNTS::SharedArray{Int64,3}, snp::Int64, BINS::Array{Float64,1}, MAF::Float64, MIN::Float64, MAX::Float64, SD::Float64)`

# Inputs
1. COUNTS [SharedArray{Int64,3}]: 3-dimensional array of allele counts for each locus for each of the 6 SNP alleles (A:T:C:G:N:DEL) for each pool
2. snp [Int64]: index number of the SNP
3. BINS [Array{Float64,1}]: vector of pool sizes
4. MAF [Float64]: minor allele frequency threshold to filter out alleles (MAF <= f > 1-MAF)
5. MIN [Float64]: minimum phenotype value
6. MAX [Float64]: maximum phenotype value
7. SD [Float64]: standard deviation of phenotype value

# Outputs
1. OUT_ALPHA: vector of additive allelic effects of alleles in the locus
2. OUT_allele: vector of alleles passing minor allele frequency (MAF) threshold filtering
3. OUT_snp: vector of repeated locus id for each of the allelic effects or allele
4, OUT_1: vector of ones and zeros indicating the indices of the alleles passing minor allele frequency (MAF) threshold filtering among the 6 SNP alleles per locus (A:T:C:G:N:DEL)
5, OUT_pA: vector of allele frequencies per locus

# Examples
```
NPOOLS=5; NLOCI=20; BINS=repeat([1.0/NPOOLS],inner=NPOOLS); MAF=0.001; MIN=0.00; MAX=10.00; SD=1.00;
COUNTS = SharedArrays.SharedArray{Int64,3}(NPOOLS, 6, NLOCI)
COUNTS[:,:,:] = Int.(round.(rand(NPOOLS, 6, NLOCI)*100))
for snp in 1:NLOCI
	out = GWAlpha_ML_module.GWAlpha_ML_iterator(COUNTS, snp, BINS, MAF, MIN, MAX, SD)
	println(out)
end
```
"""
function GWAlpha_ML_iterator(COUNTS::SharedArray{Int64,3}, snp::Int64, BINS::Array{Float64,1}, MAF::Float64, MIN::Float64, MAX::Float64, SD::Float64)
	OUT_ALPHA = convert(Array{Float64}, zeros(6))
	OUT_allele = zeros(6)
	OUT_snp = zeros(6)
	OUT_1 = zeros(6)
	OUT_pA = convert(Array{Float64}, zeros(6))
	# #parse allele counts from the sync file
	# COUNTS = zeros(Int64, NPOOLS, 6)
	# for i in 1:NPOOLS
	# 	COUNTS[i,:] = [parse(Int64, x) for x in split.(SYNC[snp, 4:(NPOOLS+3)], [':'])[i]]
	# end
	#convert to frequencies per pool
	counts = COUNTS[:,:,snp]
	FREQS = counts ./ ( sum(counts, dims=2) .+ 1e-10 ) #added 1e-10 to the denominator to avoid NAs in pools with no allele counts (zero depth; which should actually have been filtered out after mpileup using awk)
	allele_freqs = sum(FREQS .* BINS, dims=1)
	#iterate across alleles while filtering by MAF
	if (sum(counts) != 0.0)
		if (minimum(allele_freqs[allele_freqs .!= 0.0]) >= MAF) & (maximum(allele_freqs) < (1.0 - MAF)) #locus filtering by mean MAF
			for allele in 1:6
				if (allele_freqs[allele] > 0.0) & (maximum(FREQS[:,allele]) < 0.999999)  #allele filtering remove alleles with no counts and that the number of pools with allele frequency close to one should not occur even once!
				# if (sum(FREQS[:,allele] .== 0.0) < NPOOLS) & (sum(FREQS[:,allele] .> 0.999999) < 1) #filter-out alleles with at least 1 pool fixed for that allele because it causes a failure in the optimization
					freqA = FREQS[:, allele]
					pA = sum(freqA .* BINS)
					pB = 1 - pA
					BINA = (freqA .* BINS) ./ pA
					BINB = ( (1 .- freqA) .* BINS ) ./ (1-pA)
					percA = cumsum(BINA)
					percB = cumsum(BINB)

					### optimize (minimize) -log-likelihood of these major allele frequencies modelled as a beta distribution
					# using Nelder-Mead optimization or Box minimisation (try-catch if one or the other fails with preference to Nelder-Mead)
					lower_limits = [1e-20, 1e-20, 1e-20, 1e-20]
					upper_limits = [1.0, 1.0, 1.0, 1.0]
					initial_values = [0.1, 0.1, 0.1, 0.1]
					BETA = try
						Optim.optimize(beta->neg_log_likelihood_cdfbeta(beta, percA, percB), initial_values, NelderMead())
					catch
						try
							Optim.optimize(beta->neg_log_likelihood_cdfbeta(beta, percA, percB), lower_limits, upper_limits, initial_values)
						catch ### lower limits of 1e-20 to 1e-6 causes beta dist parameter values to shrink to zero somehow - so we'r'e setting lower limits to 1e-5 instead
							lower_limits = [1e-5, 1e-5, 1e-5, 1e-5]
							Optim.optimize(beta->neg_log_likelihood_cdfbeta(beta, percA, percB), lower_limits, upper_limits, initial_values)
						end
					end
					MU_A = MIN + ((MAX-MIN)*BETA.minimizer[1]/(BETA.minimizer[1]+BETA.minimizer[2]))
					MU_B = MIN + ((MAX-MIN)*BETA.minimizer[3]/(BETA.minimizer[3]+BETA.minimizer[4]))

					### compute alpha
					W_PENAL = 2*sqrt(pA*pB)
					ALPHA = W_PENAL*(MU_A - MU_B) / SD
					OUT_ALPHA[allele] =  ALPHA
					OUT_allele[allele] =  allele
					OUT_snp[allele] =  snp
					OUT_1[allele] =  1
					OUT_pA[allele] =  pA
				else
					## for explicit zero effects for null (low to none) frequency alleles
					OUT_ALPHA[allele] =  0.0
					OUT_allele[allele] =  allele
					OUT_snp[allele] =  snp
					OUT_1[allele] =  0
					OUT_pA[allele] =  allele_freqs[allele]
				end
			end
		end
	end
	return([OUT_ALPHA, OUT_allele, OUT_snp, OUT_1, OUT_pA])
end

"""
# _____________________________________________________________________________________________________________________________________________
# Genome-wide additive allelic effects (alpha) association using Pool-seq data

`GWAlpha_ML(filename_sync::String, filename_phen_py::String, MAF::Float64)`

# Inputs
1. *filename_sync* [String]: filename of the synchronized pileup file, i.e. the genotype file from pool sequencing (Pool-seq)
2. *filename_phen_py* [String]; filename of the phenotype specifications (with ".py" extension name for cross-compatibility with [GWAlpha.py](https://github.com/aflevel/GWAlpha))
3. *MAF* [Float64]: minor allele frequency threshold to filter out alleles (MAF <= f > 1-MAF)

# Output
1. Tuple of vectors of the additive allelic affects across loci with the corresponding locus identification and significance test statistics
	- CHROM
	- POS
	- ALLELE
	- FREQ
	- ALPHA
	- PVALUES
	- LOD

# Example
`OUT = GWAlpha_ML_module.GWAlpha_ML(filename_sync="test/test.sync", filename_phen_py="test/test.py", MAF=0.001)`
"""
function GWAlpha_ML(;filename_sync::String, filename_phen_py::String, MAF::Float64)
	### load the sync and phenotype files
	SYNC = DelimitedFiles.readdlm(filename_sync, '\t')
	phen = DelimitedFiles.readdlm(filename_phen_py)

	### gather phenotype specifications
	NPOOLS = length(split(phen[5], ['=', ',', '[', ']', ';'])) - 3 #less the first leading and trailing elemenets
	if length(split(phen[1], ['=', '\"'])) < 3
		global NAME = split(phen[1], ['=', '\''])[3]
	else
		global NAME = split(phen[1], ['=', '\"'])[3]
	end
	SD = parse.(Float64, split(phen[2], ['=',';'])[2])
	MIN = parse.(Float64, split(phen[3], ['=',';'])[2])
	MAX = parse.(Float64, split(phen[4], ['=',';'])[2])
	PERC = parse.(Float64, split(phen[5], ['=', ',', '[', ']', ';'])[3:(NPOOLS+1)])
	QUAN = parse.(Float64, split(phen[6], ['=', ',', '[', ']', ';'])[3:(NPOOLS+1)])
	BINS = append!([x for x in PERC], 1) - append!(zeros(1), PERC)

	### gather genotype (allele frequency) specificications
	NSNP = size(SYNC)[1]
	n_pools_sync = size(SYNC)[2] - 3
	if NPOOLS != n_pools_sync
		println("The number of pools with phenotype data does not match the number of pools with allele frequency data!")
		println("Please check you input files :-)")
		println("Remove leading and intervening whitespaces in the phenotype file.")
		exit()
	else
		n_pools_sync = nothing #clear out contents of this redundant n_pools variable
	end

	### creating an allele counts SharedArray from SYNC to minimize memory usage
	### input shared array
	COUNTS = SharedArrays.SharedArray{Int64,3}(NPOOLS, 6, size(SYNC)[1])
	progress_bar = ProgressMeter.Progress(NPOOLS, dt=1, desc="Converting the sync file into a SharedArray of allele counts: ",  barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
	for i in 1:NPOOLS
		COUNTS[i,:,:] = parse.(Int64, hcat(split.(SYNC[:, 4:(NPOOLS+3)], [':'])[:, i]...))
		ProgressMeter.update!(progress_bar, i)
	end
	### ouput shared arrays
	alpha_out = SharedArrays.SharedArray{Float64,1}(NSNP*6)
	allele_id = SharedArrays.SharedArray{Int64,1}(NSNP*6)
	locus_id = SharedArrays.SharedArray{Int64,1}(NSNP*6)
	locus_w_eff = SharedArrays.SharedArray{Int64,1}(NSNP*6)
	allele_freq = SharedArrays.SharedArray{Float64,1}(NSNP*6)

	# p = ProgressMeter.Progress(NSNP, 1, "GWAlpha_ML_iterator...", 50)
	println("Performing GWAlpha_ML in parallel...")
 	@time x = @sync @distributed for snp in 1:NSNP
		# println(snp)
		OUT = GWAlpha_ML_iterator(COUNTS, snp, BINS, MAF, MIN, MAX, SD)
		idx = collect(((6*snp)-5):(6*snp))
		alpha_out[idx] = OUT[1]
		allele_id[idx] = OUT[2]
		locus_id[idx] = OUT[3]
		locus_w_eff[idx] = OUT[4]
		allele_freq[idx] = OUT[5]
		# ProgressMeter.next!(p)
	end
	ALPHA_OUT = alpha_out[locus_w_eff .== 1]
	ALLELE_ID_INT = allele_id[locus_w_eff .== 1]
	LOCUS_ID = locus_id[locus_w_eff .== 1]
	LOCUS_W_EFF = locus_w_eff[locus_w_eff .== 1]
	ALLELE_FREQ = allele_freq[locus_w_eff .== 1]

	### estimate heuristic p-values
	P_VALUES, LOD = significance_testing_module.estimate_pval_lod(convert(Array{Float64,1}, ALPHA_OUT))
	### output
	ALLELE_ID = repeat(["N"], inner=length(ALLELE_ID_INT))
	for i in 1:length(ALLELE_ID) #convert int allele ID into corresponting A, T, C, G, N, DEL
		if ALLELE_ID_INT[i] == 1; ALLELE_ID[i] = "A"
		elseif ALLELE_ID_INT[i] == 2; ALLELE_ID[i] = "T"
		elseif ALLELE_ID_INT[i] == 3; ALLELE_ID[i] = "C"
		elseif ALLELE_ID_INT[i] == 4; ALLELE_ID[i] = "G"
		elseif ALLELE_ID_INT[i] == 5; ALLELE_ID[i] = "N"
		elseif ALLELE_ID_INT[i] == 6; ALLELE_ID[i] = "DEL"
		end
	end
	OUT = (CHROM=convert(Array{Any,1},SYNC[LOCUS_ID,1]), POS=convert(Array{Int64,1},SYNC[LOCUS_ID,2]), ALLELE=convert(Array{Any,1},ALLELE_ID), FREQ=convert(Array{Float64,1},ALLELE_FREQ), ALPHA=convert(Array{Float64},ALPHA_OUT), PVALUES=convert(Array{Float64},P_VALUES), LOD=convert(Array{Float64},LOD))
	return(OUT)
end

end # end of GWAlpha_ML_module
