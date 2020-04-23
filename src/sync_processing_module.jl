module sync_processing_module

using DelimitedFiles
using Statistics

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

# Example
`sync_processing_module.sync_filter(filename_sync="test/test.sync", MAF=0.001, DEPTH=10)`
"""
function sync_filter(;filename_sync::String, MAF::Float64=0.001, DEPTH::Int64=0)
	### test if the output file exists and use that instead of re-filtering the sync file
	filename_out = string(join(split(filename_sync, ".")[1:(end-1)], "."), "_MAF", MAF, "_DEPTH", DEPTH, ".sync")
	if (isfile(filename_out))
		println(string("Filtered sync file exists: ", filename_out))
		return(1)
	else
		### open ouput file for appending
		out = open(filename_out, "a")
		### iterate across snp
		# x = []
		@time open(filename_sync, "r") do f
			for snp in eachline(f)
				### parse counts per pool (6 alleles x n pools)
				counts = parse.(Int64, hcat(split.(split(snp, '\t')[4:end], ":")...))
				# push!(x, counts); break; end; end; counts = x[1]
				# mean depth per pool
				depth = Statistics.mean(sum(counts, dims=1))
				# minimum allele frequency
				freqs = counts ./ sum(counts, dims=1)
				freqs = freqs[.!ismissing.(freqs) .& .!isnan.(freqs) .& .!isnothing.(freqs) .& (freqs .!= 0.0)]
				if length(freqs) == 0
					min_freq = 0
					max_freq = 0
				else
					min_freq = minimum(freqs)
					max_freq = maximum(freqs)
				end
				if (min_freq >= MAF) & (max_freq < 1-MAF) & (depth >= DEPTH)
					DelimitedFiles.writedlm(out, hcat(split(snp, '\t')...), '\t')
				end
			end
		end ### faster than just using eachline() with array output
		### close output filtered sync file
		close(out)
	end
	println("===============================================================")
	println(string("Sync file filtered: ", filename_out))
	println("===============================================================")
	return(0)
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
	### test if the output file exists and use that instead of re-parsing the sync file
	filename_out = string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv")
	if (isfile(filename_out))
		println(string("Parsed sync file exists: ", filename_out))
		return(1)
	else
	### open ouput file for appending
		out = open(filename_out, "a")
		### iterate across snp
		# x = []
		allele = ["A", "T", "C", "G", "N", "DEL"]
		@time open(filename_sync, "r") do f
			for snp in eachline(f)
				### parse counts per pool (6 alleles x n pools)
				counts = parse.(Int64, hcat(split.(split(snp, '\t')[4:end], ":")...))
				### calculate allele frequencies
				freqs = counts ./ sum(counts, dims=1)
				freqs[ismissing.(freqs) .| isnan.(freqs) .| isnothing.(freqs)] .= 0.0
				# push!(x, counts); break; end; end; counts = x[1]
				chrom = repeat([split(snp, '\t')[1]], 6)
				pos = repeat([split(snp, '\t')[2]], 6)
				DelimitedFiles.writedlm(out, hcat(chrom, pos, allele, freqs), ',')
			end
		end ### faster than just using eachline() with array output
		### close output filtered sync file
		close(out)
	end
	println("===============================================================")
	println(string("Sync file parsed into a csv file: ", filename_out))
	println("===============================================================")
 	return 0
end

end #end of sync_processing_module
