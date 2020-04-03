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
#####################
"""
# ___________________________________________________________
# Parse synchronize pileup file into an allele frequency file

`sync_parse(filename_sync::String)`

Parse the synchronized pileup file line by line and appends the output into the output file.

# Output
Allele frequency csv file with the filname: `string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv")`

# Example
```
sync_parsing_module.sync_parse("test/test.sync")
```

# Note
Make sure the output file does not exist before executing this function.
**This function appends into the existing output file** with the filename: `string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv")`
"""
function sync_parse(filename_sync::String)
	### load the sync file
	sync = DelimitedFiles.readdlm(filename_sync, '\t')

	### gather genotype (allele frequency) specificications
	NSNP = size(sync)[1]
	NPOOLS = size(sync)[2] - 3

	### prepare output csv file of allele frequencies
	filename_out = string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv")
	if isfile(filename_out)
		rm(filename_out)
	end
	io = open(filename_out, "a")
	### iterate across SNPs
	COUNTS = zeros(Int64, NPOOLS, 6) #nrow=n_pools and ncol=A,T,C,G,N,DEL
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
		allele_freq_across_pools = sum(FREQS, dims=1) ./ NPOOLS
		OUT_FREQS = convert(Array{Any,2}, FREQS')
		OUT_CHR = repeat([sync[snp, 1]], inner=size(OUT_FREQS,1))
		OUT_POS = repeat([sync[snp, 2]], inner=size(OUT_FREQS,1))
		OUT_ALLELE = ["A", "T", "C", "G", "N", "DEL"]
		OUT = hcat(OUT_CHR, OUT_POS,OUT_ALLELE, OUT_FREQS)
		DelimitedFiles.writedlm(io, OUT, ',')
		ProgressMeter.update!(progress_bar, snp)
	end
	close(io)
	println("===============================================================")
	println(string("Sync file parsed into a csv file: ", filename_out))
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
