using Distributed
@everywhere using DelimitedFiles
@everywhere using Distributions
@everywhere using LinearAlgebra
@everywhere using Optim
@everywhere using ProgressMeter
@everywhere using DataFrames

@everywhere function neg_log_likelihood_cdfbeta(beta::Array{Float64,1}, data_A::Array{Float64,1}, data_B::Array{Float64,1})
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

@everywhere function GWAlpha_ML_iterator(sync::Array{Any,2}, snp::Int64, NPOOLS::Int64, BINS::Array{Float64,1}, MAF::Float64, MIN::Float64, MAX::Float64, SD::Float64)
	OUT_ALPHA = []
	OUT_allele = []
	OUT_snp = []
	OUT_1 = []
	OUT_pA = []
	#parse allele counts from the sync file
	COUNTS = zeros(Int64, NPOOLS, 6)
	for i in 1:NPOOLS
		COUNTS[i,:] = [parse(Int64, x) for x in split.(sync[snp, 4:(NPOOLS+3)], [':'])[i]]
	end
	#convert to frequencies per pool
	FREQS = COUNTS ./ ( sum(COUNTS, dims=2) .+ 1e-10 ) #added 1e-10 to the denominator to avoid NAs in pools with no allele counts (zero depth; which should actually have been filtered out after mpileup using awk)
	allele_freqs = sum(FREQS .* BINS, dims=1)
	#iterate across alleles while filtering by MAF
	if (sum(COUNTS) != 0.0)
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
					append!(OUT_ALPHA, ALPHA)
					append!(OUT_allele, allele)
					append!(OUT_snp, snp)
					append!(OUT_1, 1)
					append!(OUT_pA, pA)
				else
					### for explicit zero effects for null (low to none) frequency alleles
					# append!(ALPHA_OUT, 0.0)
					# append!(ALLELE_ID, allele)
					# append!(LOCUS_ID, snp)
					# append!(LOCUS_W_EFF, 0)
					# append!(ALLELE_FREQ, allele_freqs[allele])
				end
			end
		end
	end
	return([OUT_ALPHA, OUT_allele, OUT_snp, OUT_1, OUT_pA])
end

@everywhere function GWAlpha_ML_parallel(filename_sync::String, filename_phen_py::String, MAF::Float64)
	###test:
	# filename_sync = "test.sync"
	# filename_phen_py = "test.py"
	# MAF = 0.01
	### load the sync and phenotype files
	sync = DelimitedFiles.readdlm(filename_sync, '\t')
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
	NSNP = size(sync)[1]
	n_pools_sync = size(sync)[2] - 3
	if NPOOLS != n_pools_sync
		println("The number of pools with phenotype data does not match the number of pools with allele frequency data!")
		println("Please check you input files :-)")
		println("Remove leading and intervening whitespaces in the phenotype file.")
		exit()
	else
		n_pools_sync = nothing #clear out contents of this redundant n_pools variable
	end

	### iterate across SNPs
	alpha_out = SharedArrays.SharedArray{Float64,1}(NSNP*6)
	allele_id = SharedArrays.SharedArray{Int64,1}(NSNP*6)
	locus_id = SharedArrays.SharedArray{Int64,1}(NSNP*6)
	locus_w_eff = SharedArrays.SharedArray{Int64,1}(NSNP*6)
	allele_freq = SharedArrays.SharedArray{Float64,1}(NSNP*6)

	# p = ProgressMeter.Progress(NSNP, 1, "GWAlpha_ML_iterator...", 50)
 	@time x = @sync @distributed for snp in 1:NSNP
		# println(snp)
		OUT = GWAlpha_ML_iterator(sync, snp, NPOOLS, BINS, MAF, MIN, MAX, SD)
		idx = collect(snp:(snp+length(OUT[1])-1))
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
	# alpha_mean = mean(ALPHA_OUT[LOCUS_W_EFF .== 1])
	# alpha_sd = std(ALPHA_OUT[LOCUS_W_EFF .== 1])
	# P_VALUES = ones(length(LOCUS_ID))
	# P_VALUES[LOCUS_W_EFF .== 1] .= [pval_Normal(x, alpha_mean, alpha_sd) for x in ALPHA_OUT[LOCUS_W_EFF .== 1]]
	# LOD = -Distributions.log.(10, P_VALUES)
	P_VALUES, LOD = pval_heuristic_module.estimate_PVAL_and_LOD(convert(Array{Float64,1}, ALPHA_OUT))
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
	OUT = DataFrames.DataFrame(LOCUS_ID=LOCUS_ID, CHROM=sync[LOCUS_ID,1], POS=sync[LOCUS_ID,2], ALLELE=ALLELE_ID, FREQ=ALLELE_FREQ, ALPHA=ALPHA_OUT, PVALUES=P_VALUES, LOD=LOD)
	return(OUT)
end
