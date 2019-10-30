#####################################
###								  ###
### GENERATE POPULATION GROUPINGS ###
###								  ###
#####################################

using DelimitedFiles
using DataFrames
using CSV
using Random
using Statistics

POP_FNAME_INDI_DF = CSV.read(ARGS[1])
ALL_POP_GENO_SYNC = DelimitedFiles.readdlm(ARGS[2], '\t')
ALL_POP_PHENO_CSV = CSV.read(ARGS[3])
N_LIB = parse(Int64, ARGS[4])
### test:
# POP_FNAME_INDI_DF = CSV.read("POP_FNAME_WITHIN_INDI_DF.csv")
# ALL_POP_GENO_SYNC = DelimitedFiles.readdlm("LOLIUM_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_ALLPOP_GENO.sync", '\t')
# ALL_POP_PHENO_CSV = CSV.read("LOLIUM_10QTL_0.001mr_0.25fgs_0.00bgs_1grad_ALLPOP_PHENO.pool")
# N_LIB = parse(Int64, "400") ### squarable too

### merge randomly sampled individuals per population as defined in the GROUP_MAT grouping scheme at the "idx_pos_train" column\
function merge_individuals_across_pops_in_column_func(;POP_FNAME_INDI_DF, GROUP_MAT, idx_pops, pop_sample_size)
	### INPUTS:
	# POP_FNAME_INDI_DF	= dataframe (npop x 3): col1 (String, population ID), col2 (String, geno filename csv), col3 (String, pheno filename csv)
	# GROUP_MAT 		= array (s x m) of groups per columns under the current grouping scheme
	# idx_pops			= integer as the index of the current groups or column in the GROUP_MAT array
	# N_LIB					= total number of individuals across all poplations merged
	### extract filenames within the current group under the current grouping scheme
	fnames = Array{String}(undef, size(GROUP_MAT, 1), size(POP_FNAME_INDI_DF, 2))
	for p in 1:size(GROUP_MAT, 1)
		fnames[p, :] = convert(Matrix{String}, POP_FNAME_INDI_DF[POP_FNAME_INDI_DF.POP .== GROUP_MAT[p, idx_pops], :])
	end
	### merge data of all the populations per group
	### QUESTION!!!!! HOW DO WE DECIDE THE NUMBER OF INDIVIDUALS TO SAMPLE PER POPULATION ::: DOES IT SCALE WITH THE NUMBER OF POPULATIONS PER GROUP PER GROUPING SCHEME
	# pop_sample_size = [convert(Int, round(N_LIB / s[idx_pops]))] # sample size per population given the current sampling scheme
	for p in 1:size(fnames, 1)
		# p = 1 #test
		X_p = DelimitedFiles.readdlm(fnames[p, 2], ',')
		y_p = DelimitedFiles.readdlm(fnames[p, 3], ','; skipstart=1)[:, 3] # column 3 is the transformed phenotype data
		idx = sort!(rand(1:size(y_p, 1), pop_sample_size)) # random sampling of individuals per population
		if (p==1)
			global loci_merge = X_p[:, 1:3] #extract loci information
			global X_merge = convert(Array{Float64}, X_p[:, 4:end]'[idx, :]) #remove the 3 loci info columns --> transpose --> get random samples
			global y_merge = y_p[idx, 1]
		else
			X_merge = vcat(X_merge, convert(Array{Float64}, X_p[:, 4:end]'[idx, :]))
			y_merge = vcat(y_merge, y_p[idx, 1])
		end
	end
	pop_id_merge = repeat(fnames[:,1], inner=pop_sample_size)
	return(X_merge, y_merge, loci_merge, pop_id_merge)
end

### build input phenotype data in python script format
function build_py_phenotype_input_func(PHENO, fname)
	### INPUTS:
	# PHENO = array of phenotype data (col1:pop_id, col2:y_untransformed, col3:y)
	# fname = filename of the resulting py phenotype file to be written out
	npools = size(PHENO, 1)
	Pheno_name = string("Pheno_name='pheno'", ";")
	sig = string("sig=", sqrt(var(PHENO[:,3])), ";")
	MIN = string("MIN=", minimum(PHENO[:,3]), ";")
	MAX = string("MAX=", maximum(PHENO[:,3]), ";")
	perc = string("perc=[", join(cumsum(repeat([1.0/npools], npools-1)), ","), "];")
	q = string("q=[", join(PHENO[:,3], ","), "];")
	write_me_out = [Pheno_name, sig, MIN, MAX, perc, q]
	DelimitedFiles.writedlm(fname, write_me_out, '\t')
	return(0)
end

### random uniform sampling across the landscape
function random_uniform(;POP_FNAME_INDI_DF::DataFrames.DataFrame, ALL_POP_GENO_SYNC::Array{Any,2}, ALL_POP_PHENO_CSV::DataFrames.DataFrame, N_LIB::Int64)
	### generate a s x m_max matrix of randomised population names
	POP_NAMES = convert(Array{String}, POP_FNAME_INDI_DF.POP)
	k = length(POP_NAMES)												### total number of populations (NOTE: MUST HAVE A NATURAL NUMBER AS SQUARE-ROOT!!!)
	m_max = convert(Int, sqrt(k))										### maximum number of population-groups must be a multiple of k
	s = convert(Array{Int}, cumsum(repeat([(1/m_max)*k], inner=m_max)))	### set of number of populations per group as the number of population groups decreases
	idx_rand = Random.randperm(k)										### randomize the sequence of populations to simulate a random uniform landscape sampling of populations
	SQUARE_GROUP_MAT = reshape(POP_NAMES[idx_rand], s[1], m_max)
	### arrays of population id and corresponding genotype and phenotype filenames
	POP_ID_INDI = []
	FNAME_GENO_INDI_CSV = []
	FNAME_PHENO_INDI_CSV = []
	FNAME_GENO_COVAR_PC_CSV = []
	FNAME_GENO_COVAR_K_CSV = []
	POP_ID_POOL = []
	FNAME_GENO_POOL_SYNC = []
	FNAME_GENO_POOL_CSV = []
	FNAME_PHENO_POOL_PY = []
	FNAME_PHENO_POOL_CSV = []
	FNAME_GENO_POOL_COVAR_NPSTAT_CSV = []
	FNAME_GENO_POOL_COVAR_FST_CSV = []
	### individual genotype and phenotype data
	for i in 1:m_max
		# i = m_max
		n_groups = convert(Int64, floor(m_max*m_max/s[i]) )
		GROUP_MAT = reshape(SQUARE_GROUP_MAT[:, 1:(end - convert(Int, (m_max*m_max%s[i])/m_max) ) ], convert(Int, i*m_max), n_groups)
		pop_sample_size = convert(Int, round(N_LIB / size(GROUP_MAT, 1))) # sample size per population given the current sampling scheme
		for j in 1:size(GROUP_MAT, 2)
			# j = 1
			X_merge, y_merge, loci_merge, pop_id_merge = merge_individuals_across_pops_in_column_func(POP_FNAME_INDI_DF=POP_FNAME_INDI_DF, GROUP_MAT=GROUP_MAT, idx_pops=j, pop_sample_size=pop_sample_size)
			# fname_geno_out = string(join(GROUP_MAT[:, j], '_'), "_MERGED_GENO.csv")
			# fname_pheno_out = string(join(GROUP_MAT[:, j], '_'), "_MERGED_PHENO.csv")
			fname_geno_out = string(hash(join(GROUP_MAT[:, j])), "_MERGED_GENO.csv")
			fname_pheno_out = string(hash(join(GROUP_MAT[:, j])), "_MERGED_PHENO.csv")
			DelimitedFiles.writedlm(fname_geno_out, hcat(loci_merge, X_merge'), ',')
			# DelimitedFiles.writedlm(fname_pheno_out, hcat(pop_id_merge, y_merge), ',')
			y_df = DataFrames.DataFrame(pop=pop_id_merge, y=y_merge)
			CSV.write(fname_pheno_out, y_df)
			push!(POP_ID_INDI, join(GROUP_MAT[:, j], ';'))
			push!(FNAME_GENO_INDI_CSV, fname_geno_out)
			push!(FNAME_PHENO_INDI_CSV, fname_pheno_out)
			### covariates: to be built with by the bash  script that called this julia script
			# fname_geno_covariate_PC = string(join(GROUP_MAT[:, j], '_'), "_MERGED_GENO_COVARIATE_PC.csv")
			# fname_geno_covariate_K = string(join(GROUP_MAT[:, j], '_'), "_MERGED_GENO_COVARIATE_K.csv")
			fname_geno_covariate_PC = string(hash(join(GROUP_MAT[:, j])), "_MERGED_GENO_COVARIATE_PC.csv")
			fname_geno_covariate_K = string(hash(join(GROUP_MAT[:, j])), "_MERGED_GENO_COVARIATE_K.csv")
			push!(FNAME_GENO_COVAR_PC_CSV, fname_geno_covariate_PC)
			push!(FNAME_GENO_COVAR_K_CSV, fname_geno_covariate_K)
		end
	end
	### Pool-based genotype and phenotype data
	for i in 1:m_max
		# i = m_max
		n_groups = convert(Int64, floor(m_max*m_max/s[i]) )
		GROUP_MAT = reshape(SQUARE_GROUP_MAT[:, 1:(end - convert(Int, (m_max*m_max%s[i])/m_max) ) ], convert(Int, i*m_max), n_groups)
		for j in 1:size(GROUP_MAT, 2)
			# j = 1
			idx = [parse(Int64, x[2]) for x in split.(GROUP_MAT[:,j], 'p') ]
			Y_POOLS_MAT = hcat(repeat([round(N_LIB / length(idx))], inner=length(idx)), ALL_POP_PHENO_CSV[idx, :])
			idx_sorter = sortperm(Y_POOLS_MAT[:,3])
			Y_POOLS_MAT = Y_POOLS_MAT[idx_sorter, :]
			X_POOLS_SYNC = hcat(ALL_POP_GENO_SYNC[:,1:3], ALL_POP_GENO_SYNC[:, (idx[idx_sorter] .+ 3)])
			# fname_geno_sync_out = string(join(Y_POOLS_MAT.pop_id, '_'), "_MERGED_POOLS_GENO.sync")
			# fname_geno_csv_out = string(join(Y_POOLS_MAT.pop_id, '_'), "_MERGED_POOLS_GENO_ALLELEFREQ.csv") ### to be generated in the cross-validation julia script
			# fname_pheno_py_out = string(join(Y_POOLS_MAT.pop_id, '_'), "_MERGED_POOLS_PHENO.py")
			# fname_pheno_csv_out = string(join(Y_POOLS_MAT.pop_id, '_'), "_MERGED_POOLS_PHENO.csv")
			fname_geno_sync_out = string(hash(join(Y_POOLS_MAT.pop_id)), "_MERGED_POOLS_GENO.sync")
			fname_geno_csv_out = string(hash(join(Y_POOLS_MAT.pop_id)), "_MERGED_POOLS_GENO_ALLELEFREQ.csv") ### to be generated in the cross-validation julia script
			fname_pheno_py_out = string(hash(join(Y_POOLS_MAT.pop_id)), "_MERGED_POOLS_PHENO.py")
			fname_pheno_csv_out = string(hash(join(Y_POOLS_MAT.pop_id)), "_MERGED_POOLS_PHENO.csv")
			DelimitedFiles.writedlm(fname_geno_sync_out, X_POOLS_SYNC, '\t')
			build_py_phenotype_input_func(Y_POOLS_MAT, fname_pheno_py_out)
			DelimitedFiles.writedlm(fname_pheno_csv_out, convert(Matrix{Float64}, Y_POOLS_MAT[:,[1,4]]), ',')
			push!(POP_ID_POOL, join(Y_POOLS_MAT.pop_id, ';'))
			push!(FNAME_GENO_POOL_SYNC, fname_geno_sync_out)
			push!(FNAME_GENO_POOL_CSV, fname_geno_csv_out)
			push!(FNAME_PHENO_POOL_PY, fname_pheno_py_out)
			push!(FNAME_PHENO_POOL_CSV, fname_pheno_csv_out)
			### covariates: to be built with by the bash  script that called this julia script
			# fname_geno_covariate_PC = string(join(Y_POOLS_MAT.pop_id, '_'), "_MERGED_POOLS_GENO_COVARIATE_NPSTAT.csv")
			# fname_geno_covariate_K = string(join(Y_POOLS_MAT.pop_id, '_'), "_MERGED_POOLS_GENO_COVARIATE_FST.csv")
			fname_geno_covariate_PC = string(hash(join(Y_POOLS_MAT.pop_id)), "_MERGED_POOLS_GENO_COVARIATE_NPSTAT.csv")
			fname_geno_covariate_K = string(hash(join(Y_POOLS_MAT.pop_id)), "_MERGED_POOLS_GENO_COVARIATE_FST.csv")
			push!(FNAME_GENO_POOL_COVAR_NPSTAT_CSV, fname_geno_covariate_PC)
			push!(FNAME_GENO_POOL_COVAR_FST_CSV, fname_geno_covariate_K)
		end
	end
	### ouput dataframes
	POP_FNAME_ACROSS_INDI_DF = DataFrames.DataFrame(POP=POP_ID_INDI,
													GENO=FNAME_GENO_INDI_CSV,
													PHENO=FNAME_PHENO_INDI_CSV,
													COVARIATE_PC=FNAME_GENO_COVAR_PC_CSV,
													COVARIATE_K=FNAME_GENO_COVAR_K_CSV)
	POP_FNAME_ACROSS_POOL_DF = DataFrames.DataFrame(POP=POP_ID_POOL,
													GENO_SYNC=FNAME_GENO_POOL_SYNC,
													GENO_CSV=FNAME_GENO_POOL_CSV,
													PHENO_PY=FNAME_PHENO_POOL_PY,
													PHENO_CSV=FNAME_PHENO_POOL_CSV,
													COVARIATE_NPSTAT=FNAME_GENO_POOL_COVAR_NPSTAT_CSV,
													COVARIATE_FST=FNAME_GENO_POOL_COVAR_FST_CSV)
	return(POP_FNAME_ACROSS_INDI_DF, POP_FNAME_ACROSS_POOL_DF)
end

### pseudo-optimal square quadrants sampling
function pseudo_optimal(;POP_FNAME_INDI_DF::DataFrames.DataFrame, ALL_POP_GENO_SYNC::Array{Any,2}, ALL_POP_PHENO_CSV::DataFrames.DataFrame, N_LIB::Int64)
	### build the l x l landscape matrix of population names
	n = length(POP_FNAME_INDI_DF.POP)			### total number of populations in the landscape
	l_max = convert(Int64, floor(sqrt(n)))		### maximum number of horizontally or vertically adjacent populations in a square-group or the number of horizontally or vertically adjacent populations in each of the 4 sides of the whole square landscape
	k_max = minimum([N_LIB, n])					### maximum number of populations within a square-group which is ultimately restricted by the number of genotyping libraries (N_LIB) set
	x = [k_max]; while x[end] > 4 push!(x, convert(Int64, (sqrt(x[end])-2)^2)); end; push!(x, 1) ### all the permutations of the number of non-overalapping square-groups that can fit within the landscape
	squares_array = reverse(x) ### start with 1 square array (square-group) comprising of all n populations, followed by 4 or 9 then 16 or 25, and so on...
	LANDSCAPE_MATRIX = permutedims(convert(Array{String,2},reshape(POP_FNAME_INDI_DF.POP, l_max, l_max))) ### the arrangement of the populations across the entire square landscape
	### arrays of population id and corresponding genotype and phenotype filenames for outputting
	POP_ID_INDI = []
	FNAME_GENO_INDI_CSV = []
	FNAME_PHENO_INDI_CSV = []
	FNAME_GENO_COVAR_PC_CSV = []
	FNAME_GENO_COVAR_K_CSV = []
	POP_ID_POOL = []
	FNAME_GENO_POOL_SYNC = []
	FNAME_GENO_POOL_CSV = []
	FNAME_PHENO_POOL_PY = []
	FNAME_PHENO_POOL_CSV = []
	FNAME_GENO_POOL_COVAR_NPSTAT_CSV = []
	FNAME_GENO_POOL_COVAR_FST_CSV = []
	### iterate across all possible number of square-groups defined in the squares array
	### i.e. [1, 4, 16, 36, ..., n] for mod(n,2)==0 and [1, 9, 25, 49, ..., n] for mod(n, 2)==1
	for k in squares_array
		# k = 2
		# println(k)
		l = convert(Int64, floor(sqrt(n/k))) ### number of horizontally or verically adjacent populations per square-group
		merge_pop_names = [] ### array population names  initialization
		for quad_row in 1:convert(Int64, round(l_max/l)) ### iterate across square-groups per row
			# quad_row = 1
			row_ini = (l * (quad_row - 1)) ### first member of the square-group
			row_fin = (l * quad_row) ### last member of the square-group
			row = convert(Int64, rand([floor((row_fin-row_ini+1)/2), ceil((row_fin-row_ini+1)/2)], 1)[1]) + row_ini ### the verically middle population where if l (the number of vertically adjacent populations) is even then the "middle population is randomly chosen betweem floor(l/2) and ceil(l/2) "
			for quad_col in 1:convert(Int64, round(l_max/l))
				# quad_col = 1
				col_ini = (l * (quad_col -1))
				col_fin = (l * quad_col)
				col = convert(Int64, rand([floor((col_fin-col_ini+1)/2), ceil((col_fin-col_ini+1)/2)], 1)[1]) + col_ini ### adding 1 to account for the wobbliness of even number of populations
				push!(merge_pop_names, LANDSCAPE_MATRIX[row, col])
			end
		end
		GROUP_MAT = convert(Array{String,2}, reshape(merge_pop_names, length(merge_pop_names), 1))
		###############################
		### INDIVIDUAL GENTYPE DATA ###
		###############################
		pop_sample_size = convert(Int64, round(N_LIB/k))
		### merge
		X_merge, y_merge, loci_merge, pop_id_merge = merge_individuals_across_pops_in_column_func(POP_FNAME_INDI_DF=POP_FNAME_INDI_DF, GROUP_MAT=GROUP_MAT, idx_pops=1, pop_sample_size=pop_sample_size)
		### genotype filenames and write-out
		fname_geno_out = string(hash(join(GROUP_MAT[:, 1])), "_MERGED_GENO.csv")
		fname_pheno_out = string(hash(join(GROUP_MAT[:, 1])), "_MERGED_PHENO.csv")
		DelimitedFiles.writedlm(fname_geno_out, hcat(loci_merge, X_merge'), ',')
		### phenotype filenames and write-out
		y_df = DataFrames.DataFrame(pop=pop_id_merge, y=y_merge)
		CSV.write(fname_pheno_out, y_df)
		### covariates: to be built with by the bash  script that called this julia script
		fname_geno_covariate_PC = string(hash(join(GROUP_MAT[:, 1])), "_MERGED_GENO_COVARIATE_PC.csv")
		fname_geno_covariate_K = string(hash(join(GROUP_MAT[:, 1])), "_MERGED_GENO_COVARIATE_K.csv")
		### push id and names into the evetual output arrays
		push!(POP_ID_INDI, join(GROUP_MAT[:, 1], ';'))
		push!(FNAME_GENO_INDI_CSV, fname_geno_out)
		push!(FNAME_PHENO_INDI_CSV, fname_pheno_out)
		push!(FNAME_GENO_COVAR_PC_CSV, fname_geno_covariate_PC)
		push!(FNAME_GENO_COVAR_K_CSV, fname_geno_covariate_K)
		##############################
		### POOLED POPULATION DATA ###
		##############################
		if k != 1 ### exclude the grouping with just a single population
			### extract population phenotype data sort and extract sync genotype data
			idx = [parse(Int64, x[2]) for x in split.(GROUP_MAT[:,1], 'p') ]
			Y_POOLS_MAT = hcat(repeat([round(N_LIB / length(idx))], inner=length(idx)), ALL_POP_PHENO_CSV[idx, :])
			idx_sorter = sortperm(Y_POOLS_MAT[:,3])
			Y_POOLS_MAT = Y_POOLS_MAT[idx_sorter, :]
			X_POOLS_SYNC = hcat(ALL_POP_GENO_SYNC[:,1:3], ALL_POP_GENO_SYNC[:, (idx[idx_sorter] .+ 3)])
			### genotype and phenotype filenames, build phenotype files and write-out
			fname_geno_sync_out = string(hash(join(Y_POOLS_MAT.pop_id)), "_MERGED_POOLS_GENO.sync")
			fname_geno_csv_out = string(hash(join(Y_POOLS_MAT.pop_id)), "_MERGED_POOLS_GENO_ALLELEFREQ.csv") ### to be generated in the cross-validation julia script
			fname_pheno_py_out = string(hash(join(Y_POOLS_MAT.pop_id)), "_MERGED_POOLS_PHENO.py")
			fname_pheno_csv_out = string(hash(join(Y_POOLS_MAT.pop_id)), "_MERGED_POOLS_PHENO.csv")
			DelimitedFiles.writedlm(fname_geno_sync_out, X_POOLS_SYNC, '\t')
			build_py_phenotype_input_func(Y_POOLS_MAT, fname_pheno_py_out)
			DelimitedFiles.writedlm(fname_pheno_csv_out, convert(Matrix{Float64}, Y_POOLS_MAT[:,[1,4]]), ',')
			### covariates: to be built with by the bash  script that called this julia script
			fname_geno_covariate_PC = string(hash(join(Y_POOLS_MAT.pop_id)), "_MERGED_POOLS_GENO_COVARIATE_NPSTAT.csv")
			fname_geno_covariate_K = string(hash(join(Y_POOLS_MAT.pop_id)), "_MERGED_POOLS_GENO_COVARIATE_FST.csv")
			### push id and names into the evetual output arrays
			push!(POP_ID_POOL, join(Y_POOLS_MAT.pop_id, ';'))
			push!(FNAME_GENO_POOL_SYNC, fname_geno_sync_out)
			push!(FNAME_GENO_POOL_CSV, fname_geno_csv_out)
			push!(FNAME_PHENO_POOL_PY, fname_pheno_py_out)
			push!(FNAME_PHENO_POOL_CSV, fname_pheno_csv_out)
			push!(FNAME_GENO_POOL_COVAR_NPSTAT_CSV, fname_geno_covariate_PC)
			push!(FNAME_GENO_POOL_COVAR_FST_CSV, fname_geno_covariate_K)
		end
	end
	### ouput dataframes
	POP_FNAME_ACROSS_INDI_DF = DataFrames.DataFrame(POP=POP_ID_INDI,
													GENO=FNAME_GENO_INDI_CSV,
													PHENO=FNAME_PHENO_INDI_CSV,
													COVARIATE_PC=FNAME_GENO_COVAR_PC_CSV,
													COVARIATE_K=FNAME_GENO_COVAR_K_CSV)
	POP_FNAME_ACROSS_POOL_DF = DataFrames.DataFrame(POP=POP_ID_POOL,
													GENO_SYNC=FNAME_GENO_POOL_SYNC,
													GENO_CSV=FNAME_GENO_POOL_CSV,
													PHENO_PY=FNAME_PHENO_POOL_PY,
													PHENO_CSV=FNAME_PHENO_POOL_CSV,
													COVARIATE_NPSTAT=FNAME_GENO_POOL_COVAR_NPSTAT_CSV,
													COVARIATE_FST=FNAME_GENO_POOL_COVAR_FST_CSV)
	return(POP_FNAME_ACROSS_INDI_DF, POP_FNAME_ACROSS_POOL_DF)
end

### EXECUTION
@time POP_FNAME_ACROSS_INDI_DF, POP_FNAME_ACROSS_POOL_DF = random_uniform(POP_FNAME_INDI_DF=POP_FNAME_INDI_DF, ALL_POP_GENO_SYNC=ALL_POP_GENO_SYNC, ALL_POP_PHENO_CSV=ALL_POP_PHENO_CSV, N_LIB=N_LIB)
CSV.write("POP_FNAME_ACROSS_INDI_RANDUNIF_DF.csv", POP_FNAME_ACROSS_INDI_DF)
CSV.write("POP_FNAME_ACROSS_POOL_RANDUNIF_DF.csv", POP_FNAME_ACROSS_POOL_DF)

@time POP_FNAME_ACROSS_INDI_DF, POP_FNAME_ACROSS_POOL_DF = pseudo_optimal(POP_FNAME_INDI_DF=POP_FNAME_INDI_DF, ALL_POP_GENO_SYNC=ALL_POP_GENO_SYNC, ALL_POP_PHENO_CSV=ALL_POP_PHENO_CSV, N_LIB=N_LIB)
CSV.write("POP_FNAME_ACROSS_INDI_PSEUDOPT_DF.csv", POP_FNAME_ACROSS_INDI_DF)
CSV.write("POP_FNAME_ACROSS_POOL_PSEUDOPT_DF.csv", POP_FNAME_ACROSS_POOL_DF)
