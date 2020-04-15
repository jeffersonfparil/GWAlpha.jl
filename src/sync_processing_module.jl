module sync_processing_module

using DelimitedFiles
using Statistics
using ProgressMeter

"""
# ________________________________________
# Filter synchronized pileup genotype data

`sync_filter(;filename_sync::String, MAF::Float64, DEPTH::Int64=0)`

Filter a synchronized pileup file to include only sites
with user defined minimum allele frequency (MAF) and minimum sequencing depth (DEPTH).

# Input
1. *filename_sync* [String]: [synchronized pileup filename](https://sourceforge.net/p/popoolation2/wiki/Manual/)
2. *MAF* [Float64]: minimum allele frequency threshold (default=0.001)
3. *DEPTH* [Int64]: minimum depth threshold (default=0)

# Output:
1. filtered sync file with the suffix: `string("_MAF", MAF, "_DEPTH", DEPTH, ".sync")`
2. boolean indices of sync loci * 6 alleles passing the MAF and DEPTH filtering (can be used to filter parsed sync genotype data [sync_parsing_module.sync_parse()])

# Example
`idx = sync_processing_module.sync_filter(filename_sync="test/test.sync", MAF=0.001, DEPTH=10);`
"""
function sync_filter(;filename_sync::String, MAF::Float64=0.001, DEPTH::Int64=0)
	### test if an output file exists and use that instead of re-filtering the sync file
	filename_out = string(join(split(filename_sync, ".")[1:(end-1)], "."), "_MAF", MAF, "_DEPTH", DEPTH, ".sync")
	filename_out_idx = string(join(split(filename_sync, ".")[1:(end-1)], "."), "_MAF", MAF, "_DEPTH", DEPTH, "_IDX_OUT.txt")
	if (isfile(filename_out)) & (isfile(filename_out_idx))
		OUT = convert(Array{Bool,1}, DelimitedFiles.readdlm(filename_out_idx)[:,1])
		println(string("Using the existing filtered sync file: ", filename_out))
		println(string("And the corresponding existing filtered sync file indices: ", filename_out_idx))
	else
		### load the sync and phenotype files
		sync = DelimitedFiles.readdlm(filename_sync, '\t')

		### gather genotype (allele frequency) specifications
		NSNP = size(sync)[1]
		NPOOLS = size(sync)[2] - 3

		### iterate across SNPs
		OUT = repeat([false], inner=(NSNP*6)) ### boolean vector of loci that passed MAF and DEPTH filtering
		COUNTS = zeros(Int64, NPOOLS, 6) #nrow=n_pools and ncol=A,T,C,G,N,DEL
		progress_bar = ProgressMeter.Progress(NSNP, dt=1, desc="Filter sync by MAF Progress: ",  barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow) #progress bar
		for snp in 1:NSNP
	        # snp = 1
			# println(snp)
			#parse allele counts from the sync file
			for i in 1:NPOOLS
				COUNTS[i,:] = [parse(Int64, x) for x in split.(sync[snp, 4:(NPOOLS+3)], [':'])[i]]
			end
	        # test if all pools have the desired depth or higher
	        if sum(sum(COUNTS, dims=2) .>= DEPTH) == NPOOLS
	    		#convert to frequencies per pool
	    		FREQS = COUNTS ./ ( sum(COUNTS, dims=2) .+ 1e-10 ) #added 1e-10 to the denominator to avoid NAs in pools with no allele counts (zero depth; which should actually have been filtered out after mpileup using awk)
	    		allele_freqs = mean(FREQS, dims=1)
	    		#iterate across alleles while filtering by MAF
	    		if length(allele_freqs[allele_freqs .!= 0.0]) != 0
	    			if (minimum(allele_freqs[allele_freqs .!= 0.0]) >= MAF) & (maximum(allele_freqs) < (1.0 - MAF)) #locus filtering by mean MAF
	    				for allele in 1:6
	    					if (allele_freqs[allele] > 0.0) & (maximum(FREQS[:,allele]) < 0.999999)  #allele filtering remove alleles with no counts and that the number of pools with allele frequency close to one should not occur even once!
	    						OUT[((snp-1)*6) + allele] = true
	    					end
	    				end
	    			end
	    		end
	    	end
	        # println("$snp")
	        ProgressMeter.update!(progress_bar, snp)
		end
		#write out filtered sync
		idx_in_sync = maximum(reshape(OUT, 6, NSNP), dims=1)'[:,1]
		DelimitedFiles.writedlm(filename_out, sync[idx_in_sync, :])
		DelimitedFiles.writedlm(filename_out_idx, OUT)
	end
	return(OUT)
end

"""
# ____________________________________________________________
# Parse synchronized pileup file into an allele frequency file

`sync_parse(filename_sync::String)`

Parse the synchronized pileup file line by line and appends the output into a comma-separated (csv) file.

# Output
Allele frequency csv file with the filename: `string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv")`

# Example
`sync_processing_module.sync_parse("test/test.sync")`
"""
function sync_parse(filename_sync::String)
	### load the sync file
	sync = DelimitedFiles.readdlm(filename_sync, '\t')

	### gather genotype (allele frequency) specifications
	NSNP = size(sync)[1]
	NPOOLS = size(sync)[2] - 3

	### prepare output csv file of allele frequencies
	filename_out = string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv")
	if isfile(filename_out)
		# rm(filename_out)
		println(string("Parsed sync file exists: ", filename_out))
		println("You can use the existing file or you can remove it and re-parse the sync file.")
		return(0)
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


end #end of sync_processing_module
