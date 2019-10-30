###############################################
###											###
### Parsing synchronised pileup (SYNC) file ###
###									        ### NOTE!!!! I APPEND!!!!!! 20190925
############################################### same module name and filename

module sync_parsing_module

#####################
###				  ###
### load packages ###
###				  ###
#####################
using DelimitedFiles
using ProgressMeter
using DataFrames
using CSV

#####################
###				  ###
### main function ###
###				  ###
##################### #input is the genotype sync filename
function sync_parse(filename_sync::String)
	### load the sync file
	sync = DelimitedFiles.readdlm(filename_sync, '\t')

	### gather genotype (allele frequency) specificications
	NSNP = size(sync)[1]
	NPOOLS = size(sync)[2] - 3

	### iterate across SNPs
	COUNTS = zeros(Int64, NPOOLS, 6) #nrow=n_pools and ncol=A,T,C,G,N,DEL
	# global OUT = DataFrames.DataFrame()
	progress_bar = ProgressMeter.Progress(NSNP, dt=1, desc="Progress: ",  barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow) #progress bar
	for snp in 1:NSNP
		# snp = 65
		# println(snp)
		#parse allele counts from the sync file
		for i in 1:NPOOLS
			COUNTS[i,:] = [parse(Int64, x) for x in split.(sync[snp, 4:(NPOOLS+3)], [':'])[i]]
		end
		#convert to frequencies per pool (FREQS: npools x (nalleles=6))
		FREQS = COUNTS ./ ( sum(COUNTS, dims=2) .+ 1e-10 ) #added 1e-10 to the denominator to avoid NAs in pools with no allele counts (zero depth; which should actually have been filtered out after mpileup using awk)
		# FREQS = COUNTS ./ sum(COUNTS, dims=2)
		allele_freq_across_pools = sum(FREQS, dims=1) ./ NPOOLS
		# #filter by MAF
		# non_zero_index = convert(Array{Bool}, allele_freq_across_pools .!= 0.0)
		# allele_freq_across_pools_non_zero = allele_freq_across_pools[non_zero_index]
		# if (isempty(allele_freq_across_pools_non_zero) == false) #for locui without any allele counts due to popoolations' filtering by read quality?! not sure but they exist!
		# 	if (minimum(allele_freq_across_pools_non_zero) > MAF) & (maximum(allele_freq_across_pools_non_zero) < (1.0 - MAF)) # filter by MAF and with allele freq ~1.0 that slipped
				# non_zero_nonref_index = convert(Array{Bool}, allele_freq_across_pools_non_zero .!= maximum(allele_freq_across_pools_non_zero))
				# OUT_FREQS = FREQS[:, non_zero_index[:]][:, non_zero_nonref_index[:]]' #transpose so we can append each locus (can have multiple alleles per locus) to the next line
				# OUT_FREQS = FREQS[:, non_zero_index[:]]' #transpose so we can append each locus including the reference allele!!! (can have multiple alleles per locus) to the next line
				OUT_FREQS = FREQS'
				# idx_remove_major_allele = allele_freq_across_pools_non_zero .!= maximum(allele_freq_across_pools_non_zero) #remove major allele to allocate 1 degree of freedom to the mean thereby reducing the dimesionality of the predictions substantially
				# OUT_FREQS = OUT_FREQS[idx_remove_major_allele, :]  #remove major allele to allocate 1 degree of freedom to the mean thereby reducing the dimesionality of the predictions substantially
				OUT_CHR = sync[snp, 1]
				OUT_POS = sync[snp, 2]
				# OUT_ALLELE = ["A", "T", "C", "G", "N", "DEL"][non_zero_index[:]]
				OUT_ALLELE = ["A", "T", "C", "G", "N", "DEL"]
				# OUT_ALLELE = ["A", "T", "C", "G", "N", "DEL"][non_zero_index[:]][idx_remove_major_allele]  #remove major allele to allocate 1 degree of freedom to the mean thereby reducing the dimesionality of the predictions substantially
				# OUT_REFALLELE = ["A", "T", "C", "G", "N", "DEL"][convert(Array{Bool}, allele_freq_across_pools .== maximum(allele_freq_across_pools))[:]]
				for a in 1:size(OUT_FREQS)[1] #iterate across alleles
					# global OUT = DataFrames.DataFrame(CHROM=OUT_CHR, POS=OUT_POS, REF=OUT_REFALLELE, ALLELE=OUT_ALLELE[a])
					global OUT = DataFrames.DataFrame(CHROM=OUT_CHR, POS=OUT_POS, ALLELE=OUT_ALLELE[a])
					# global OUT = DataFrames.DataFrame(CHROM=OUT_CHR, POS=OUT_POS, REF=OUT_REFALLELE, ALLELE=OUT_ALLELE[a], POOL1=OUT_FREQS[a,1])
					# insert!(OUT, OUT_CHR, :CHROM)
					# insert!(OUT, OUT_POS, :POS)
					# insert!(OUT, OUT_REFALLELE, :REF)
					# insert!(OUT, OUT_ALLELE[a], :ALLELE)
					# for p in 2:NPOOLS #add allele frequency per pool iteratively into the dataframe
					for p in 1:NPOOLS #add allele frequency per pool iteratively into the dataframe
						insert!(OUT, size(OUT)[2]+1, OUT_FREQS[a,p], Symbol("POOL", p))
					end
					CSV.write(string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv"), OUT, append=true)
				end
		# 	end
		# end
		ProgressMeter.update!(progress_bar, snp)
	end
	println("===============================================================")
 	return 0
end

end #end of GWAlpha module

### Sample usage:
# include("src/sync_parsing_module.jl")
# # DIR="/data/Lolium/Quantitative_Genetics/Genomic_prediction_2018Lolium_50s"
# # DIR="/data/Lolium/Population_Genetics/2018November_Ryegrass_Collection/VCF"
# DIR="/data/Lolium/Quantitative_Genetics/Genomic_prediction_2018Lolium_50s"
# ALL_FILES=readdir(DIR)
# SYNC_FILES=ALL_FILES[occursin.(".sync", ALL_FILES)]
# # for sync in string(DIR, "/") .* SYNC_FILES
# # 	println(sync)
# # 	INPUT = [sync, 0.01]
# # 	sync_parse(INPUT)
# # end
# individual sync file execution:
# INPUT = ["IF.sync", 1/500.00]
# sync_parse(INPUT)
# # testing parallele execution
# using Distributed
# function multiThreading_test(DIR, SYNC_FILES)
# 	@Distributed.sync for sync in string(DIR, "/") .* SYNC_FILES
# 		println(sync)
# 		INPUT = [sync, 0.01]
# 		@Distributed.async sync_parse(INPUT)
# 	end
# end
# multiThreading_test(DIR, SYNC_FILES)
# #tests
# INPUT = ["/data/Lolium/Quantitative_Genetics/Genomic_prediction_2018Lolium_50s/IF.sync", 0.001]
# sync_parse(INPUT)
# INPUT = ["/data/Lolium/Population_Genetics/2018November_Ryegrass_Collection/VCF/Lolium2019_100X_30PHRED_filtered.sync", 0.01]
# sync_parse(INPUT)
# # giving up on julia in-house parallelization solutions and just joing with GNU-parallel
# # copy the library loading and function definitions into a new file: "_parallelize_me_SYNC_PARSE_.jl"
# # then add these to the bottom:
# sync_fname = ARGS[1]
# MAF = parse(Float64, ARGS[2])
# INPUT = [sync_fname, MAF]
# sync_parse(INPUT)
# # then execute in bash via
# parallel julia _parallelize_me_SYNC_PARSE_.jl {} 0.002 ::: $(ls *.sync)
