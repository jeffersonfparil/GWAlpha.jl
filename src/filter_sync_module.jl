module filter_sync_module

using DelimitedFiles
using Statistics
using ProgressMeter

### LEGACY MAF FILTERING
function filter_sync_by_MAF(filename_sync, MAF)
	### load the sync and phenotype files
	sync = DelimitedFiles.readdlm(filename_sync, '\t')

	### gather genotype (allele frequency) specificications
	NSNP = size(sync)[1]
	NPOOLS = size(sync)[2] - 3

	### iterate across SNPs
	OUT = repeat([false], inner=(NSNP*6))
	COUNTS = zeros(Int64, NPOOLS, 6) #nrow=n_pools and ncol=A,T,C,G,N,DEL
	progress_bar = ProgressMeter.Progress(NSNP, dt=1, desc="Filter sync by MAF Progress: ",  barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow) #progress bar
	for snp in 1:NSNP
		# println(snp)
		#parse allele counts from the sync file
		for i in 1:NPOOLS
			COUNTS[i,:] = [parse(Int64, x) for x in split.(sync[snp, 4:(NPOOLS+3)], [':'])[i]]
		end
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
		# println("$snp")
		ProgressMeter.update!(progress_bar, snp)
	end
	return(OUT)
end

### UPDATED MAF AND DEPTH FILTERING FUNCTION
function filter_sync(;filename_sync::String, MAF::Float64, DEPTH::Int64=0)
	### load the sync and phenotype files
	sync = DelimitedFiles.readdlm(filename_sync, '\t')

	### gather genotype (allele frequency) specificications
	NSNP = size(sync)[1]
	NPOOLS = size(sync)[2] - 3

	### iterate across SNPs
	OUT = repeat([false], inner=(NSNP*6))
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
	filename_out = string(join(split(filename_sync, ".")[1:(end-1)], "."), "_MAF", MAF, "_DEPTH", DEPTH, ".sync")
	idx_in_sync = maximum(reshape(OUT, 6, NSNP), dims=1)'[:,1]
	DelimitedFiles.writedlm(filename_out, sync[idx_in_sync, :])
	return(OUT)
end

end #end of filter_sync_module

### INPUTS:
### (1). synchronized mpileup filename
### (2). minimum allele frequency threshold (MAF)
### (3). minimum depth threshold
### OUTPUT:
### (1) filtered sync file ((_MAF${MAF}_DEPTH${DEPTH}.sync))
### (2) indices of sync loci passing the MAF-DEPTH filtering