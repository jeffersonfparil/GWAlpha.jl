module neuralnet_5layers

### load packages
using Statistics, Distributions, UnicodePlots, Flux

# ### naively simulate Pool-seq data and phenotypes
# n = 60
# # m = Int64(1e5)
# m = 100
# X = Array{Float32, 2}(undef, n, m)
# X_ref = abs.(rand(n, Int64(m/2))) #biallelic
# X_alt = 1 .- X_ref
# X[:, 1:2:m] = X_ref
# X[:, 2:2:m] = X_alt
# nQTL = 10
# posQTL = sort!(sample(collect(1:m), nQTL))
# b = zeros(m)
# b[posQTL] = rand(Distributions.Chisq(2), nQTL)
# err = rand(Distributions.Normal(0, 1), n)
# y = (X*b) .+ err
# UnicodePlots.histogram(y)

### tests:
using DelimitedFiles
cd("/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia")
X_raw = DelimitedFiles.readdlm("VCF/ACC_MAPQ20_BASQ20_MAF0.001_DEPTH10_ALLELEFREQ.csv", ',')
y_raw = DelimitedFiles.readdlm("ACC_pheno.csv", ',')
X_noFilter = convert(Array{Float32}, X_raw[:, 4:end]')
MAF = 0.01
allele_freqs = mean(X_noFilter, dims=1)
allele_filter_idx = ( (allele_freqs .> MAF) .& (allele_freqs .< (1.00 - MAF)) )[1,1:end]
LOCI_ID = X_raw[allele_filter_idx, 1:3]
X = X_noFilter[:, allele_filter_idx]
y = convert(Array{Float32}, y_raw[:, 2])
UnicodePlots.histogram(y)

function neuralnet_5layers(X::Array{Float32,2}, y::Array{Float32,1}, nEpochs::Int64=1000)
    ### build the neural network model
    nDataPoints = size(X, 1)
    if nDataPoints != size(y, 1)
        println("The number of datapoints in the genotype data is inconsitent with that of the phenotype data!")
        return(1)
    end
    nInputs = size(X, 2)
    nOutputs = size(y, 2)
    nHiddenLayers = 3 ### conceptually we are trying to simulate DNA -> RNA -> protein -> DNA-RNA-protein-interactions -> phenotype
    nNodes_l1 = Int(round(nInputs/nHiddenLayers))
    nNodes_l2 = Int(round(nNodes_l1/nHiddenLayers))
    nNodes_l3 = Int(round(nNodes_l2/nHiddenLayers))
    # nNodes_l4 = Int(round(nNodes_l3/nHiddenLayers))
    # nNodes_l5 = Int(round(nNodes_l4/nHiddenLayers))
    model = Flux.Chain(
                        Flux.Dense(nInputs, nNodes_l1, σ),
                        Flux.Dense(nNodes_l1, nNodes_l2, σ),
                        Flux.Dense(nNodes_l2, nNodes_l3, σ),
                        Flux.Dense(nNodes_l3, nOutputs)
            )
    ### define the cost function as a simple mean square error
    cost_function(input, output) = Flux.mse(model(input), output)
    ### define the trainable paramaters list
    params_list = Flux.params(model)
    ### define the optimizer
    optimizer = Flux.Descent(0.1)
    ### define the training data iterator
    data_iterator = Flux.Iterators.repeated((X', y'), Int64(ceil(size(X,1)/2)))
    @time Flux.@epochs nEpochs Flux.train!(cost_function, params_list, data_iterator, optimizer)
    return(model)
end
model = neuralnet_5layers(X, y, 100)

### test
model(X')
UnicodePlots.scatterplot(model(X')[1,:], y)
cor(model(X')[1,:], y)
Flux.mse(model(X')[1,:], y)


end ### end of the neuralnet_5layers module
