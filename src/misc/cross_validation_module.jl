#####################################################
###												  ###
### Cross-validation of genomic prediction models ###
###												  ###
#####################################################

module cross_validation_module

#####################
### load packages ###
#####################
using DelimitedFiles
using DataFrames
using Statistics
using GLM
using ProgressMeter
include("PoolGPAS_module.jl")
# JULIA_SCRIPT_HOME = "/data/Lolium/Softwares/genomic_prediction/src"
# push!(LOAD_PATH, JULIA_SCRIPT_HOME)
# using PoolGPAS_module
# ##########################################
# ### test:
# filename_sync = "/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/VCF/ACC_MAPQ20_BASQ20_MAF0.001_DEPTH10.sync"
# filename_phen_csv = "/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/ACC_pheno.csv"
# MAF = 0.001
# DEPTH = 10
# COVARIATE = DelimitedFiles.readdlm("/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/VCF/ACC_MAPQ20_BASQ20_MAF0.001_DEPTH10_COVARIATE_FST.csv", ',')
# NITER = 10
# NFOLD = 5
# ##########################################

############################
### function definitions ###
############################
### assess genomic prediction accuracy
function asses_model_accuracy_func(y_true, y_pred)
	if (var(y_pred) > 1.0e-10)
		# mean deviation of the predcited from the true phenotypic values
		MEAN_DEVIANCE = Statistics.mean(abs.(y_true .- y_pred)) #in percentage unit since the units we're using is in percentage eh!
		VAR_DEVIANCE = var(abs.(y_true .- y_pred)) #in percentage unit since the units we're using is in percentage eh!
		# Pearson' product moment correlation
		CORR = cor(y_true, y_pred)
		# modeling the true phenotypic values as a function of the predicted values
		data = DataFrames.DataFrame(y_true=convert(Array{Float64}, y_true), y_pred=convert(Array{Float64}, y_pred))
		# model = GLM.fit(GLM.LinearModel, GLM.@formula(y_true ~ 0 + y_pred), data) # fixed intercept at the origin because we assume a slope of one for perfect genomic prediction accuracy
		model = GLM.fit(GLM.LinearModel, GLM.@formula(y_true ~ y_pred), data)
		if length(coef(model)[:,1]) == 1
			local INTERCEPT = 0.0
		else
			local INTERCEPT = coef(model)[1,1]
		end
		SLOPE = coef(model)[end,1]
		R2 = var(GLM.predict(model)) / var(y_true) # not sure if this is correct
		RMSD = sqrt( (sum(y_pred .- y_true)^2)/length(y_true) ) #root mean suare deviation or RMSE (E for error)
	else
		MEAN_DEVIANCE = nothing
		VAR_DEVIANCE = nothing
		CORR = nothing
		INTERCEPT = nothing
		SLOPE = nothing
		R2 = nothing
		RMSD = nothing
	end
	# output: we want percent deviance to be zerol; correlation to be 1; intercept to be zero; and slope to be 1
	out = DataFrames.DataFrame(MEAN_DEVIANCE=MEAN_DEVIANCE, VAR_DEVIANCE=VAR_DEVIANCE, CORRELATION=CORR, INTERCEPT=INTERCEPT, SLOPE=SLOPE, R2=R2, RMSD=RMSD)
	return(out)
end

##########################
###					   ###
### EXECUTIVE FUNCTION ###
###					   ###
##########################
"""
# ____________________________________________
# Cross-validation of genomic predition models

`cross_validation()`

Cross-validate genomic prediction models...

# Input
1.

# Output
1.

# Examples
```
filename_sync = "test/test.sync"
filename_phen_csv = "test/test.csv"
MAF = 0.01
DEPTH = 10
@time OUT1 = PoolGPAS(filename_sync, filename_phen_py, MAF);
@time OUT2 = PoolGPAS(filename_sync, filename_phen_csv, MAF, DEPTH, MODEL="FIXED_LS", COVARIATE=nothing);
```
"""
function cross_validation(filename_sync::String, filename_phen_csv::String, MAF::Float64, DEPTH::Int64; COVARIATE=nothing, NITER=100, NFOLD=5)
	# # ##########################################
	# # ### TESTS:
	# filename_sync = "/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/VCF/ACC_MAPQ20_BASQ20_MAF0.001_DEPTH10.sync"
	# filename_phen_csv = "/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/ACC_pheno.csv"
	# MAF = 0.001
	# DEPTH = 10
	# COVARIATE = DelimitedFiles.readdlm("/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/VCF/ACC_MAPQ20_BASQ20_MAF0.001_DEPTH10_COVARIATE_FST.csv", ',')
	# NITER = 10
	# NFOLD = 5
	# # ##########################################

	### Load input data
	SYNC = DelimitedFiles.readdlm(filename_sync, '\t')
		SYNC_LOCI = SYNC[:, 1:3]
		SYNC_COUNTS = SYNC[:, 4:end]
	GENO = try
		DelimitedFiles.readdlm(string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv"), ',')
	catch
		sync_parsing_module.sync_parse(filename_sync)
		DelimitedFiles.readdlm(string(join(split(filename_sync, ".")[1:(end-1)], '.'), "_ALLELEFREQ.csv"), ',')
	end
	LOCI_ID = GENO[:, 1:3]
	FREQ = GENO[:, 4:end]
	PHEN = DelimitedFiles.readdlm(filename_phen_csv, ',')

	ITERATION_OUT = []
	# TRAIT_OUT = []
	FOLD_OUT = []
	TRAINING_SIZE_OUT = []
	VALIDATION_SIZE_OUT = []
	MODEL_OUT = []
	N_PREDICTORS_OUT = []
	MEAN_DEVIANCE_OUT = []
	VAR_DEVIANCE_OUT = []
	CORR_OUT = []
	SLOPE_OUT = []
	R2_OUT = []
	RMSD_OUT = []
	MODELS_LIST = ["FIXED_LS", "FIXED_RR", "FIXED_GLMNET", "FIXED_LASSO"]

	pb = ProgressMeter.Progress(NITER * NFOLD * length(MODELS_LIST), 1, "Cross-validation in progress...", 50)
	counter = [1]
	for j in 1:NITER
		# j = 1
		idx_nonmissing = collect(1:size(PHEN, 1))[.!isequal.(PHEN[:, 2], missing)]
		sub_geno = FREQ[:, idx_nonmissing]
		sub_phen = PHEN[idx_nonmissing, :]
		sub_cova = try
			COVARIATE[idx_nonmissing, :] ### NOTE: Z matrices are recomputed for the whole training and validation sets in teal-world applicatios since we're using the full relationships between training and cross-validation
		catch
			nothing
		end
		pop_id_vec = rand(collect(1:size(sub_phen,1)), size(sub_phen,1))[1:(end - (size(sub_phen,1) % NFOLD))] ### forcing the cross-validation matrix grouping rectangular
		n_sub = length(pop_id_vec)
		grouping_mat = reshape(pop_id_vec, NFOLD, convert(Int, n_sub/NFOLD))
		group_size = convert(Int, round(n_sub/NFOLD))
		for l in 1:NFOLD
			# l = 1
			RANDOM_ID = string(Int(round(rand() * rand() * time()/60)))
			idx_training = reshape(grouping_mat[(collect(1:NFOLD) .!= l), :], group_size*(NFOLD-1), 1)[:,1]
			idx_validation = grouping_mat[l, :]
			n_training = length(idx_training)
			n_validation = length(idx_validation)
			SYNC_training_fname = string("training_", j, "iter_", l, "fold_", RANDOM_ID,"_GENO.sync")
				DelimitedFiles.writedlm(SYNC_training_fname, hcat(SYNC_LOCI, SYNC_COUNTS[:, (idx_nonmissing)][:, (idx_training)]), '\t')
				DelimitedFiles.writedlm(replace(SYNC_training_fname, ".sync"=>"_ALLELEFREQ.csv"), hcat(LOCI_ID, sub_geno[:, (idx_training)]), ',')
			X_validation = sub_geno[:, (idx_validation)]'
			PHEN_training_fname = string("training_", j, "iter_", l, "fold_", RANDOM_ID,"_PHEN.csv")
				DelimitedFiles.writedlm(PHEN_training_fname, sub_phen[idx_training, :], ',')
			y_validation = sub_phen[idx_validation, 2]
			Z_training = try
				sub_cova[idx_training, :] ### NOTE: Z matrices are recomputed for the whole training and validation sets in teal-world applicatios since we're using the full relationships between training and cross-validation
			catch
				nothing
			end
			Z_validation = try
				sub_cova[idx_validation, :] ### NOTE: Z matrices are recomputed for the whole training and validation sets in teal-world applicatios since we're using the full relationships between training and cross-validation
			catch
				nothing
			end
			# TRAINING AND VALIDATION
			for model in MODELS_LIST
				# model = "FIXED_LS"
				ProgressMeter.update!(pb, counter[1])
				counter[1] = counter[1] + 1
				eff, covareff = PoolGPAS_module.PoolGPAS(SYNC_training_fname, PHEN_training_fname, MAF, DEPTH, MODEL=model, COVARIATE=Z_training);
				b = vcat(eff.BETA[1], covareff, eff.BETA[2:end])
				# io = open(string("CROSS_VALIDATION_5_FOLD_BETAS_", basename(filename_phen_csv)), "a")
				# DelimitedFiles.writedlm(io, b', ',')
				# close(io)
				y_pred = try
					hcat(repeat([1.0], inner=n_validation), Z_validation, X_validation[:, eff.LOCUS_ID[2:end]]) * b
				catch
					hcat(repeat([1.0], inner=n_validation), X_validation) * b
				end
				ASSESSMENT = asses_model_accuracy_func(y_validation, y_pred)
				push!(ITERATION_OUT, j)
				# push!(TRAIT_OUT, trait)
				push!(FOLD_OUT, l)
				push!(TRAINING_SIZE_OUT, n_training)
				push!(VALIDATION_SIZE_OUT, n_validation)
				push!(MODEL_OUT, model)
				push!(N_PREDICTORS_OUT, sum(abs.(b) .> 0.00)-1) #less intercept
				push!(MEAN_DEVIANCE_OUT, ASSESSMENT.MEAN_DEVIANCE[1])
				push!(VAR_DEVIANCE_OUT, ASSESSMENT.VAR_DEVIANCE[1])
				push!(CORR_OUT, ASSESSMENT.CORRELATION[1])
				push!(SLOPE_OUT, ASSESSMENT.SLOPE[1])
				push!(R2_OUT, ASSESSMENT.R2[1])
				push!(RMSD_OUT, ASSESSMENT.RMSD[1])
			end
			### clean up
			rm(SYNC_training_fname)
			rm(replace(SYNC_training_fname, ".sync"=>"_ALLELEFREQ.csv"))
			rm(PHEN_training_fname)
			rm.(string.(replace(SYNC_training_fname, "GENO.sync"=>"PHEN-"), MODELS_LIST, "_Manhattan.png"))
			rm.(string.(replace(SYNC_training_fname, "GENO.sync"=>"PHEN-"), MODELS_LIST, "_Alphas.csv"))
		end
	end
	# using JLD2
	# JLD2.@load "CROSS_VALIDATION_5_FOLD_OUTPUT.jld2"
	### convert nothing into missing
	for i in [ITERATION_OUT, FOLD_OUT, TRAINING_SIZE_OUT, VALIDATION_SIZE_OUT, MODEL_OUT, N_PREDICTORS_OUT, MEAN_DEVIANCE_OUT, VAR_DEVIANCE_OUT, CORR_OUT, SLOPE_OUT, R2_OUT, RMSD_OUT]
	 i[i .== nothing] .= missing
	end
	OUT = DataFrames.DataFrame(ITERATION=ITERATION_OUT, FOLD=FOLD_OUT, TRAINING_SIZE=TRAINING_SIZE_OUT, VALIDATION_SIZE=VALIDATION_SIZE_OUT,
								MODEL=MODEL_OUT, N_PREDICTORS=N_PREDICTORS_OUT, MEAN_DEVIANCE=MEAN_DEVIANCE_OUT, VAR_DEVIANCE=VAR_DEVIANCE_OUT,
								CORR=CORR_OUT, SLOPE=SLOPE_OUT, R2=R2_OUT, RMSD=RMSD_OUT)
	return(OUT)
end

end #end of cross_validation_module

# ############################################################################################
# ### SAMPLE EXECUTION
# filename_sync = "/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/VCF/ACC_MAPQ20_BASQ20_MAF0.001_DEPTH10.sync"
# filename_phen_csv = "/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/ACC_pheno.csv"
# MAF = 0.001
# DEPTH = 10
# COVARIATE = DelimitedFiles.readdlm("/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/VCF/ACC_MAPQ20_BASQ20_MAF0.001_DEPTH10_COVARIATE_FST.csv", ',')
# NITER = 1
# NFOLD = 5
# COVARIATE = DelimitedFiles.readdlm("/data/Lolium/Quantitative_Genetics/GWAS_GP_2019_SE_Australia/VCF/ACC_MAPQ20_BASQ20_MAF0.001_DEPTH10_COVARIATE_FST.csv", ',')
# @time OUT = cross_validation(filename_sync, filename_phen_csv, MAF, DEPTH, COVARIATE=COVARIATE, NITER=NITER, NFOLD=NFOLD)
# ### write-out
# using Dates
# using CSV
# fname_out = string(replace(filename_sync, ".sync" => "_CROSS_VALIDATION_"), Dates.now(), "_", Int(round(rand()*1e10)), ".csv")
# CSV.write(fname_out, OUT)
