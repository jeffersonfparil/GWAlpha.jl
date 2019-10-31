###########################
### Building covariates ###
###########################

#####################
###				  ###
### load packages ###
###				  ###
#####################
using DelimitedFiles
using Statistics
using LinearAlgebra
using ProgressMeter
using DataFrames

#####################
###				  ###
### main function ###
###				  ### #inputs are the genoype matrix of allele dose (0,1,2), the covariate you wish to generate: PC or K, and the filename of the output file to store the covariate
##################### #outputs are the covariate matrix as julia array and as a csv file
function build_covariates(;X_raw, MAF=0.001, covariate_out="PC")
	if covariate_out == "PC"
		### principal components
		n = size(X_raw, 1)
		l = size(X_raw, 2)
		### MAF filtering
		LOC = collect(1:l)[ ((mean(X_raw, dims=1) ./ 2) .> MAF)[:] .& ((mean(X_raw, dims=1) ./ 2) .< (1.0-MAF))[:] ]
		X = X_raw[:, LOC]
		X_std = (X .- mean(X, dims=1)) ./ std(X, dims=1)
		X_SVD = LinearAlgebra.svd(X)
		COVARIATE = X_SVD.U * LinearAlgebra.Diagonal(X_SVD.S)
	elseif covariate_out == "K"
		### kinship matrix (equation 1 of Goudet, Kay & Weir, 2018)
		idx_biallelic = [] # using only biallelic loci
		for i in 1:convert(Int, round(size(X_raw,2)/5))
			allele_sums = sum( X_raw[:, (((i-1)*5)+1):(i*5)], dims=1 )
			to_retain_1_allele = collect(1:5) .* (allele_sums .== maximum(allele_sums))[1,:]
			idx = to_retain_1_allele .== maximum(to_retain_1_allele)
			append!(idx_biallelic, idx)
		end
		X_biallelic = X_raw[:, convert(Array{Bool}, idx_biallelic)]
		### MAF filtering
		l = size(X_biallelic, 2)
		LOC = collect(1:l)[ ((mean(X_biallelic, dims=1) ./ 2) .> MAF)[:] .& ((mean(X_biallelic, dims=1) ./ 2) .< (1.0-MAF))[:] ]
		X_biallelic_filtered = X_biallelic[:, LOC]
		n = size(X_biallelic_filtered, 1)
		M = Array{Float64}(undef, n, n)
		# pb = ProgressMeter.Progress(n; dt=0.1, desc="Building Kinship matrix: ", barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
		for i in 1:n
			for j in 1:n
				M[i,j] = sum( (1 .+ ( (X_biallelic_filtered[i,:] .- 1) .* (X_biallelic_filtered[j,:] .- 1) )) ./ 2 )
			end
			# ProgressMeter.update!(pb, i)
		end
		M_nodiag = copy(M)
		M_nodiag[diagind(M_nodiag)] .= 0
		Ms = sum(M_nodiag) / (n*(n-1))
		COVARIATE = (M .- Ms) ./ (1 .- Ms)
	else
		println("Please specify a covariate_out.")
		println("Choose 'PC' for principal componenets of eigenvectors or 'K' for kinship matrix")
		COVARIATE = nothing
		return(0)
	end
	# return(COVARIATE)
	return(COVARIATE)
end

FNAME=ARGS[1]
MAF = parse(Float64, ARGS[2])
COVARIATE = ARGS[3]
# FNAME="test_1kloci_g1000_p01_GENO.csv"
# MAF = 0.001
# COVARIATE = "K"
# COVARIATE = "PC"
X_raw = convert(Array{Float64}, DelimitedFiles.readdlm(FNAME, ',')[:, 4:end]')
Z = build_covariates(X_raw=X_raw, MAF=MAF, covariate_out=COVARIATE)
fname_out = string( join(split(basename(FNAME), '.')[1:(end-1)], '.'), "_COVARIATE_", COVARIATE, ".csv" )
DelimitedFiles.writedlm(fname_out, Z, ',')
