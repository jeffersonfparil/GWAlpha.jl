####################################################################
###                                                              ###
###              Convert quantinemo2_individual                  ###
### genotype output files into a single synchronized pileup file ###
###         to simulation Population-level Pool-seq              ###
###                                                              ###
####################################################################

################
#### inputs ####
################
#### fname_list_txt = text file listing the filenames of the genotype data output of quantinemo2_geno_sim_01a_OUTPUTPARSE.jl
#### nLoci = number of loci in the genotype files
#### nAlleles = number of alleles per loci in the genotype files
#################
#### outputs ####
#################
#### (1) Synchronized pileup file (nLoci x nPop+3) (NO HEADER)

### load libraries
using CSV
using DataFrames
using DelimitedFiles
using Statistics


### load the text file listing the genotype data and the other input paramteres
fname_list_txt = ARGS[1]
nLoci = parse(Int64, ARGS[2])
nAlleles = parse(Int64, ARGS[3])
pheno_fname = ARGS[4]
pheno_min_max_fname = ARGS[5]
### test
# fname_list_txt = "/data/Lolium/Quantitative_Genetics/genomic_prediction_simulations/QUANTINEMO2_SIM/fname_list.txt"
# nLoci = 2000
# nAlleles = 5
# pheno_fname = "/data/Lolium/Quantitative_Genetics/genomic_prediction_simulations/QUANTINEMO2_SIM/test_1kloci_2019-08-23_04-16-33/PHENOTYPES_WITH_POPID.txt"
# pheno_min_max_fname = "/data/Lolium/Quantitative_Genetics/genomic_prediction_simulations/QUANTINEMO2_SIM/QTL_min_max_GEBV.spec"
geno_list = DelimitedFiles.readdlm(fname_list_txt)
nPop = length(geno_list)
phe = CSV.read(pheno_fname, delim="\t", header=["pop_id", "y_untr"])
min_max_phe = CSV.read(pheno_min_max_fname, delim="\t")

### mean phenotypic values per population and write out
MIN_y = minimum([min_max_phe.MIN_GEBV[1], minimum(phe.y_untr)]) #extreme values based on both the predicted and observed phenotypes,
MAX_y = maximum([min_max_phe.MAX_GEBV[1], maximum(phe.y_untr)]) #i.e. observed phenotypes may be lower than the expected minimum or higher than the expecter maximum
phe.y = ( phe.y_untr .- MIN_y ) ./ (MAX_y - MIN_y)
PHENO_AGGREGATE = aggregate(phe, :pop_id, mean)
pop_id_string = []
for i in 1:size(PHENO_AGGREGATE, 1) ### modify population ID (pop_id) into srtings with "p" prefix to align to cross-validation scripts
    nPop = maximum(PHENO_AGGREGATE.pop_id)
    nrepZERO = convert(Int, floor(log(10, nPop))) - convert(Int, floor(log(10, PHENO_AGGREGATE.pop_id[i])))
    push!(pop_id_string, string("p", join(push!(repeat([0], inner=nrepZERO), PHENO_AGGREGATE.pop_id[i]))))
end
PHENO_AGGREGATE.pop_id = string.(pop_id_string)
fname_split = split(basename(geno_list[1]), "_")
# CSV.write(string(dirname(geno_list[1]), "/", fname_split[1], "_", fname_split[2], "_ALLPOP_PHENO.pool"), PHENO_AGGREGATE, delim=',')
CSV.write(string(dirname(geno_list[1]), "/", join(fname_split[1:(end-3)], "_"), "_ALLPOP_PHENO.pool"), PHENO_AGGREGATE, delim=',')

### iterate across populations
println("Convert individual population genotype data into a single synchronised pileup file.")
SYNC = Array{String}(undef, nLoci, nPop)
for i in 1:nPop
    f = geno_list[i]
    dat = DelimitedFiles.readdlm(f, ',')
    sums = Statistics.sum(dat[:,4:end], dims=2) #generate a vector of sums across alleles and loci
    sums = reshape(sums, (nAlleles, nLoci))' #reshape into a matrix of m=nLoci and n=5 alleles: A-T-C-G-DEL
    sums = hcat(sums, sums[:,end]) #add the N allele such that we have A-T-C-G-N-DEL
    sums[:,5] .= 0 #set 0 to the whole N column
    ### concatenate allele counts per locus
    sync_column = Array{String}(undef, nLoci)
    for i in 1:length(sync_column)
        sync_column[i] = join(string.(sums[i,:]), ":")
    end
    SYNC[:,i] = sync_column
end

### output sync file
CHR = reshape(DelimitedFiles.readdlm(geno_list[1], ',')[:,1], (nAlleles, nLoci))'[:,1] #extract chromosome info
POS = reshape(DelimitedFiles.readdlm(geno_list[1], ',')[:,2], (nAlleles, nLoci))'[:,1] #extract position info
REF = repeat(["A"], inner=nLoci) #set the reference allele as just all A for simplicity eh!
OUT = hcat(CHR, POS, REF, SYNC)
fname_split = split(basename(geno_list[1]), "_")
# DelimitedFiles.writedlm(string(dirname(geno_list[1]), "/", fname_split[1], "_", fname_split[2], "_ALLPOP_GENO.sync"), OUT, '\t')
DelimitedFiles.writedlm(string(dirname(geno_list[1]), "/", join(fname_split[1:(end-3)], "_"), "_ALLPOP_GENO.sync"), OUT, '\t')
