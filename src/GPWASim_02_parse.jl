############################################################
###                                                      ###
### PARSING QuantiNemo2 SIMULATION PARAMETERS AND OUTPUT ###
###                                                      ###
############################################################

### NOTE: FOR MOST FUNCTIONS I HAVE HARD-CODED NUMBER OF ALLELES TO 5 --> NEED TO MAKE THIS MORE FLEXIBLE!!!!

################
#### inputs ####
################
#### dat_fname =
#### genome_loci_spec_fname = loci specifications (HEADER: CHROM, r, POS)


### load libraries
using CSV
using DataFrames
using ProgressMeter
using DelimitedFiles
using Statistics
using StatsBase

### Input:
dat_fname = ARGS[1]                 ### quantinemo2 genotype FSTAT file output
genome_loci_spec_fname = ARGS[2]    ### loci specifications (HEADER: CHROM, r, POS)
QTL_phen_min_max_fname = ARGS[3]    ### minimum and maximum possible genomic breeding values
npools = parse(Int, ARGS[4])        ### number of pools per population
read_count = parse(Int, ARGS[5])    ### read depth to simulate for building pileup files from the FSTAT genotype output file
nAlleles = parse(Int64, split(open(readlines, dat_fname)[1], " ")[3])
### test
# cd("/data/Lolium/Quantitative_Genetics/genomic_prediction_simulations/QUANTINEMO2_SIM")
# outdir = readdir()[end]
# dat_fname = string(outdir, "/test_1kloci_g140_p01.dat")
# genome_loci_spec_fname = "GENOME_SPEC.csv"
# QTL_phen_min_max_fname = "QTL_min_max_GEBV.spec"
# npools = 5
# read_count = 100
# nAlleles = parse(Int64, split(open(readlines, dat_fname)[1], " ")[3])
#### Output:
#### (1) Genotype matrix (nind x nloci*5) (NO HEADER)
#### (2) Phenotype data (nind x 3) (HEADER: individual ID, untransformed phenotypic value, transformed phenotyp value to range from 0-1 based on possible extremes based on allele effects)
#### (3) Pooling percentiles (npools+1 x 1) (NO HEADER)
#### (4) Mean phenotypic values per pool (npools x 1) (NO HEADER)
#### (5) Allele frequency matrix (npools x nloci*5) (NO HEADER)
#### (6) Synchronized pileup file for the whole population (nloci x npools+3) (NO HEADER)

#################################################################################################################################
############################
### parse the FSTAT file ###
############################
println("Parsing the genotype data.")
### SUBFUNCTION: convert the 2-digit numbers into strings, separate them and assign the corresponding SNP allele
function num_2digit_2_bin_5allele(num_2digit)
    # #test
    # num_2digit = 23
    #initialize the ouput array of 5 elements corrensponding to the 5 possible SNP alleles: A,T,C,G and deletion
    out = [0, 0, 0, 0, 0]
    #convert the 2-digit number allele form into a string
    str_2digit = string(num_2digit)
    str_arr = split(str_2digit, "")
    for str in str_arr
        if str == "1"
            out = out .+ [1, 0, 0, 0, 0] #A
        elseif str == "2"
            out = out .+ [0, 1, 0, 0, 0] #T
        elseif str == "3"
            out = out .+ [0, 0, 1, 0, 0] #C
        elseif str == "4"
            out = out .+ [0, 0, 0, 1, 0] #G
        elseif str == "5"
            out = out .+ [0, 0, 0, 0, 1] #DEL
        end
    end
    return(out)
end
function FUNC_PARSE_DAT(dat_fname)
    ### load the data in FSTAT format and extract the number of loci, number of alleles per locus and the number of individuals
    row1 = open(readlines, dat_fname)[1] # read the first row
    nloci = parse(Int64, split(row1, " ")[2])
    nalleles = parse(Int64, split(row1, " ")[3])
    nlines = countlines(dat_fname)
    header = ["IND"]
    append!(header, ["loc_"] .* string.(collect(1:nloci)))
    append!(header, ["DELETE_ME"]) #for the extra column weirdly included when reading with space as delimiter
    # @time dat = CSV.read(dat_fname, delim=" ", header=header, datarow=nloci + 2)
    @time dat = DelimitedFiles.readdlm(dat_fname, ' '; skipstart=nloci+1)
    dat = dat[:,1:(nloci+1)] #remove the extra column
    nind = size(dat)[1]

    ### rename the individual column with consecutive numbers
    # dat.IND = collect(1:nind)
    dat[:,1] = collect(1:nind)

    ###### iterate across loci: convert the 2-digit numbers into strings, separate them and assign the corresponding SNP allele
    GENO = Array{Int32}(undef, nind, nloci*5) #initialize the genotype matrix output
    progress_bar = ProgressMeter.Progress(nloci, dt=1, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow) #progress bar
    for i in 1:nloci
        # i = 1
        geno = dat[:, i+1]
        for j in 1:length(geno) #iterate across individuals within the locus at hand
            # j = 1
            _start_ = (5*(i-1)) + 1
            _end_ = 5*i
            GENO[j, _start_:_end_] = num_2digit_2_bin_5allele(geno[j])
        end
        ProgressMeter.update!(progress_bar, i)
    end
    ### output the genotype matrix
    return(GENO)
end
### INPUT:
###### (1) FSTAT (*.dat) genotype data (filename)
### OUTPUT:
###### (1) Genotype matrix (nind x nloci*5) (NO HEADER)

### execute the fulnction and save the genotype matrix: nind x nloci*5
# dat_fname = "QUANTI_g140.dat"
GENO = FUNC_PARSE_DAT(dat_fname)
# @time DelimitedFiles.writedlm(string(split(dat_fname, ".dat")[1], "_GENO.csv"), GENO, ',') # delay writing out and merge with loci info first
#################################################################################################################################


#################################################################################################################################
########################################################################
### merge loci info and genotype data for a consolidated GENO output ###
########################################################################
println("Merging the loci info and genotype data.")
LOCI_SPEC = CSV.read(genome_loci_spec_fname)
CHROM = repeat(LOCI_SPEC.CHROM, inner=nAlleles)
POS = repeat(LOCI_SPEC.POS, inner=nAlleles)
SNPs_alleles = ["A", "T", "C", "G", "DEL", "N"]
ALLELES = repeat(SNPs_alleles[1:nAlleles], outer=nrow(LOCI_SPEC))
GENO_PLUS_LOCI_SPEC = hcat(CHROM, POS, ALLELES, GENO')
@time DelimitedFiles.writedlm(string(split(dat_fname, ".dat")[1], "_GENO.csv"), GENO_PLUS_LOCI_SPEC, ',')
#################################################################################################################################


#################################################################################################################################
##############################
### collect phenotype data ###
##############################
println("Parsing the phenotype data.")
function FUNC_PARSE_PHENO(phe_fname, min_max_phe_fname)
    # phe_fname = "QUANTI_g140_p1.phe"
    phe = CSV.read(phe_fname, header=["IND", "y_untr"], datarow=3, delim="\t")
    nind = nrow(phe)
    phe.IND = collect(1:nind)
    # min_max_phe_fname = "../QTL_min_max_GEBV.spec"
    min_max_phe = CSV.read(min_max_phe_fname, delim="\t")
    MIN_y = minimum([min_max_phe.MIN_GEBV[1], minimum(phe.y_untr)]) #extreme values based on both the predicted and observed phenotypes,
    MAX_y = maximum([min_max_phe.MAX_GEBV[1], maximum(phe.y_untr)]) #i.e. observed phenotypes may be lower than the expected minimum or higher than the expecter maximum
    y = ( phe.y_untr .- MIN_y ) ./ (MAX_y - MIN_y)
    phe.y = y
    return(phe)
end
### INPUTS:
###### (1) phenotype (*.phe) file (filename)allele_counter[1] > 0
###### (2) QTL specifications (*.spec) (filename)
### OUTPUT:
###### (1) phenotye data frame: IND (individual ID number); y_untr (untransformed phenotypic value); y (transformed phenotypic value ranging from 0 to 1 based on the minimum and maximum possible genotypic value based on the allele effects)
###### (WITH HEADER)

### execute and save
# QTL_phen_min_max_fname = "QTL_min_max_GEBV.spec"
phe_fname = string(split(dat_fname, ".dat")[1], ".phe")
PHENO = FUNC_PARSE_PHENO(phe_fname, QTL_phen_min_max_fname)
CSV.write(string(split(phe_fname, ".phe")[1], "_PHENO.csv"), PHENO, delim=',')
#################################################################################################################################

#################################################################################################################################
##########################################################################
### build the sync file and the calculate the mean phenotypes per pool ###
##########################################################################
println("Pooling.")
function FUNC_POOLING(GENO, PHENO, npools, pool_size)
    ### extract the dataset dimensionsallele_counter[1] > 0
    nind = size(GENO)[1]
    nloci = convert(Int64, size(GENO)[2] / 5)
    if (nind != nrow(PHENO))
        println("Oh no! The number of individuals do match between the genotype and phenotype data!")
        println("Exiting the function!")
        return(999)
    end
    ### set the number of pools to divide the simulated dataset and the pools sizes if you want it to be manually set
    # npools = 5
    # pool_size = convert(Int64, round(nind / npools)) #equally-sized pools
    # # pool_size = 100

    ### group the individuals into npools groups, each equally sized: pool_size
    PERCENTILES = append!([0.0], Statistics.quantile(PHENO.y, [1.0/npools] .* collect(1:npools)))
    PHENO_POOLS = Array{Float64}(undef, npools, 1)
    GENO_POOLS = Array{Float64}(undef, npools, nloci*5)
    IDX_POOLING = Array{Int64}(undef, pool_size, npools)
    for i in 1:npools
        test = (PHENO.y .>= PERCENTILES[i]) .& (PHENO.y .<= PERCENTILES[i+1])
        idx = StatsBase.sample(collect(1:nind)[test], pool_size, replace=true) #simulate random sampling during leaf sampling, DNA extraction, library preparation, and sequencing
        PHENO_POOLS[i] = Statistics.mean(PHENO.y[idx])
        GENO_POOLS[i,:] = Statistics.mean(GENO[idx,:], dims=1) ./ 2 #diploid that's why we need to divide by 2
        IDX_POOLING[:, i] = idx #the index of individuals per pool
    end

    ### simulate a synchronized pileup file (*.sync) using these pool datasets
    DEPTH = 100 #just a fixed arbitrary simulted sequencing or read depth
    chr = LOCI_SPEC.CHROM
    pos = LOCI_SPEC.POS
    ref = []; alleles = ["A", "T", "C", "G", "DEL"]
    sync = Array{String}(undef, nloci, npools)
    pb = ProgressMeter.Progress(nloci, dt=1, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
    for i in 1:nloci
        _start_ = (5*(i-1)) + 1
        _end_ = 5*i
        subdat = convert(Array{Int64}, round.( DEPTH .* GENO_POOLS[:, _start_:_end_]' )) #5 (A,T,C,G,DEL) x npools
        for j in 1:npools
            sync[i, j] = string(string(join(string.(subdat[1:end-1,j]), ":"), ":0"), ":", string.(subdat[end,j])) #insert N allele so that: A:T:C:G:N:DEL
        end
        # save the reference allele
        allele_sums = sum(subdat, dims=2)
        push!(ref, alleles[(maximum(allele_sums) .== allele_sums)[:]][end])
        ProgressMeter.update!(pb, i)
    end
    SYNC = hcat(chr, pos, ref, sync)
    return(PERCENTILES, PHENO_POOLS, GENO_POOLS, SYNC, IDX_POOLING)
end
### INPUTS:
###### (1) Genotype matrix (nind x nloci*5) --> output of FUNC_PARSE_DAT()
###### (2) Phenotye data frame --> output of FUNC_PARSE_PHENO()
###### (3) number of pools
###### (4) size of each pool (equally-sized pools)
### OUTPUTS:
###### (1) phenotype percentiles per pool starting at 0.0 (npools+1 x 1) (NO HEADER)
###### (2) mean phenotype per pool (npools x 1) (NO HEADER)
###### (3) allele frequencies per pool matrix (npools x nloci*5) (NO HEADER)
###### (4) allele frequencies per pool in synchronized pileup format (nloci x npools+3) (*.sync) (NO HEADER)
###### (5) indices of individuals sampled per pool (poolsize x npools) (NO HEADER)

### execute and save
nind = size(GENO)[1] #or nind = size(PHENO)[1]
pool_size = convert(Int, round(nind/npools))
PERCENTILES, PHENO_POOLS, GENO_POOLS, SYNC, IDX_POOLING = FUNC_POOLING(GENO, PHENO, npools, pool_size)
DelimitedFiles.writedlm(string(split(dat_fname, ".dat")[1], "_POOLS_PERCENTILES.csv"), PERCENTILES, ',')
DelimitedFiles.writedlm(string(split(dat_fname, ".dat")[1], "_POOLS_PHENO.csv"), hcat(repeat([pool_size], npools), PHENO_POOLS), ',')
DelimitedFiles.writedlm(string(split(dat_fname, ".dat")[1], "_POOLS_GENO.csv"), GENO_POOLS, ',')
DelimitedFiles.writedlm(string(split(dat_fname, ".dat")[1], "_POOLS_GENO.sync"), SYNC, '\t')

### extra: build the GWAlpha.py phenotype input file:
Pheno_name = string("Pheno_name='pheno'", ";")
sig = string("sig=", sqrt(var(PHENO.y)), ";")
MIN = string("MIN=", minimum(PHENO.y), ";")
MAX = string("MAX=", maximum(PHENO.y), ";")
perc = string("perc=[", join(cumsum(repeat([1.0/npools], npools-1)), ","), "];")
q = string("q=[", join(PERCENTILES[2:end], ","), "];")
write_me_out = [Pheno_name, sig, MIN, MAX, perc, q]
DelimitedFiles.writedlm(string(split(dat_fname, ".dat")[1], "_POOLS_PHENO.py"), write_me_out, '\t')
#################################################################################################################################


#################################################################################################################################
############################################################
### build pileup files from FSTAT (*.dat) genotype files ###
############################################################
### allele format conversion from number string to character string and vice versa
function allele_num_2_char_viceversa_func(num_or_char, num2char_char2num)
    if num2char_char2num == 0
    ### for character to number string conversion
        if num_or_char == "A"
            return("1")
        elseif num_or_char == "T"
            return("2")
        elseif num_or_char == "C"
            return("3")
        elseif num_or_char == "G"
            return("4")
        end
    elseif num2char_char2num == 1
    ### for number to UPPERCASE character string conversion
        if num_or_char == "1"
            return("A")
        elseif num_or_char == "2"
            return("T")
        elseif num_or_char == "3"
            return("C")
        elseif num_or_char == "4"
            return("G")
        end
    elseif num2char_char2num == 2
        ### for number to LOWERCASE character string conversion
            if num_or_char == "1"
                return("a")
            elseif num_or_char == "2"
                return("t")
            elseif num_or_char == "3"
                return("c")
            elseif num_or_char == "4"
                return("g")
            end
    end
end
# INPUTS:
# num_or_char = [A,T,C,G] or ["1","2","3","4"]
# num2char_char2num = [0 for char2num,
#                       1 num2char UPPERCASE,
#                       2 num2char LOWERCASE]
### convert fstat numeric format into pileup reads format
function num_2digit_2_pileup_read_func(vec_num_2digit, ref, read_count)
    # #test
    # i = 1
    # vec_num_2digit = dat[:,i+1]
    # ref = FUNC_allele_num_2_char_viceversa(col3_ref[i], 0) #char to num string
    # read_count = col4_readcount[i]
    #initializ the ouput array of 5 elements corrensponding to the 5 possible SNP alleles: A,T,C,G and deletion
    out = repeat([".", ","], outer=convert(Int, read_count/2)) #alternating forward and reverse strands
    #convert the 2-digit number allele form into a string
    vec_str_2digit = string.(vec_num_2digit)
    mat_str_arr = split.(vec_str_2digit, "")
    # iterate across individuals
    for j in 1:convert(Int, read_count/2)
        # println(j)
        str_arr = mat_str_arr[j]
        idx_allele1 = (2*(j-1)) + 1
        idx_allele2 = 2*j
        # forward strand (first digit)
        if str_arr[1] != ref
            if str_arr[1] == "5" #if DEL
                out[idx_allele1] = "1" *  allele_num_2_char_viceversa_func(ref, 1) #uppercase
            else
                out[idx_allele1] = allele_num_2_char_viceversa_func(str_arr[1], 1) #uppercase
            end
        end
        # reverse strand (second digit)
        if str_arr[2] != ref
            if str_arr[2] == "5" #if DEL
                out[idx_allele2] = "1" *  allele_num_2_char_viceversa_func(ref, 2) #lowercase
            else
                out[idx_allele2] = allele_num_2_char_viceversa_func(str_arr[2], 2) #lowercase
            end
        end
    end
    return(join(out))
end
# INPUTS:
# vec_num_2digit = vector of the number string across individuals for the ith locus
# ref = reference allele in number string format for the ith locus
# read_count = read counts i.e. number of indiduals * 2 (diploid)
function dat_to_pileup_func(dat_fname, LOCI_SPEC, read_count)
    nLoci = nrow(LOCI_SPEC)
    dat = DelimitedFiles.readdlm(dat_fname, ' '; skipstart=nLoci+1)[:,1:(nLoci+1)] #remove the extra column weirdly included when reading with space as delimiter
    nInd = size(dat)[1]
    col1_chr = LOCI_SPEC.CHROM
    col2_pos = LOCI_SPEC.POS
    col3_ref = sample(["A", "T", "C", "G"], nLoci) #simulate random reference alleles: A,T,C & G
    ### prepare the reads info (columns 4, 5 and 6: reads, base_quality)
    # row1 = open(readlines, dat_fname)[1] # read the first row
    # header = ["IND"]
    # append!(header, ["loc_"] .* string.(collect(1:nLoci)))
    # append!(header, ["DELETE_ME"]) #for the extra column weirdly included when reading with space as delimiter
    # col4_readcount = repeat([nInd*2], inner=nLoci) #simulate read count as just as the number of individuals *2 (perfect representation of each individuals haploid genome)
    col4_readcount = repeat([read_count], inner=nLoci) #simulate read count
    ### convert fstat genotypes into pileup reads format
    col5_reads = Array{String}(undef, nLoci)
    progress_bar = ProgressMeter.Progress(nLoci, dt=1, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
    for i in 1:nLoci
        # println(i)
        # vec_num_2digit = dat[:,i+1]
        vec_num_2digit = dat[StatsBase.sample(collect(1:nInd), read_count) ,i+1] # random sampling read_count
        ref = allele_num_2_char_viceversa_func(col3_ref[i], 0) #char to num string
        col5_reads[i] = num_2digit_2_pileup_read_func(vec_num_2digit, ref, read_count)
        ProgressMeter.update!(progress_bar, i)
    end
    ### simulate base quality column as constant PHRED=90=="c"
    col6_quality = repeat([join(repeat(["c"], inner=read_count))], inner=nLoci)
    ### merge columns into the pileup output file and save as a tab-delimited pileup file
    PILEUP = hcat(col1_chr, col2_pos, col3_ref, col4_readcount, col5_reads, col6_quality)
    ### output
    return(PILEUP)
end

### build pileup file for the whole population
println("Building pileup file from whole population genotype data.")
WHOLE_POP_PILEUP = dat_to_pileup_func(dat_fname, LOCI_SPEC, read_count)
DelimitedFiles.writedlm(string(split(dat_fname, ".dat")[1], "_POPULATION.pileup"), WHOLE_POP_PILEUP, '\t')

### build pileup file for each pool in the population
println("Building pileup files for each pool per population.")
nLoci = nrow(LOCI_SPEC)
header = DelimitedFiles.readdlm(dat_fname, ' ')[1:(nLoci+1), 1:(nLoci+1)]
dat = DelimitedFiles.readdlm(dat_fname, ' '; skipstart=nLoci+1)[:,1:(nLoci+1)] #remove the extra column weirdly included when reading with space as delimiter
for i in 1:npools
    println(string("Pool: ", i))
    sub = dat[IDX_POOLING[:, i], :]
    pool_dat_fname = string(split(dat_fname, ".dat")[1], "_POOL", i, ".dat")
    DelimitedFiles.writedlm(pool_dat_fname, vcat(header, sub), ' ')
    POOL_PILEUP = dat_to_pileup_func(pool_dat_fname, LOCI_SPEC, read_count)
    DelimitedFiles.writedlm(string(split(pool_dat_fname, ".dat")[1], ".pileup"), POOL_PILEUP, '\t')
end
#################################################################################################################################
