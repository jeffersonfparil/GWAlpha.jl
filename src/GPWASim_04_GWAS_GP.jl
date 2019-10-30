##################################################
###						                       ###
### Genomic Prediction Models Cross-Validation ###
###						                       ###
##################################################

########################################
###									 ###
### load packages and custom modules ###
###									 ###
########################################
using DataFrames
using CSV
using DelimitedFiles
using Statistics
using GLM
using LinearAlgebra
using MultipleTesting
using ProgressMeter
using Random
using Plots; Plots.pyplot()
using ColorBrewer
JULIA_SCRIPT_HOME = @__DIR__
# # test:
# JULIA_SCRIPT_HOME = "/data/Lolium/Softwares/genomic_prediction/src"
push!(LOAD_PATH, JULIA_SCRIPT_HOME)
using GWAS_module
using GP_module
using sync_parsing_module
using filter_sync_module
using PoolGPAS_module

############################
###						 ###
### function definitions ###
###						 ###
############################

### genomic prediction iterative models
function predict_ITERATIVE_MODELS_func(X_test, loci_info, int_and_covariate_effects, beta_trained, y_trained, X_trained; COVARIATE_X_trained=nothing, COVARIATE_X_test=nothing)
	# ###########################################
	# # ### tests:
	# X_test=X_test
	# loci_info=TRAINING_OUT[i].LOCI
	# int_and_covariate_effects=TRAINING_OUT[i].INTCOVAR_EFF
	# beta_trained=TRAINING_OUT[i].EFF
	# y_trained=y_train
	# X_trained=X_train
	# COVARIATE_X_trained=TRAINING_OUT[i].COVARIATE
	# COVARIATE_X_test=COVARIATE_X_test=COVARIATE_X_test
	# ###########################################

	### data specs
	n_test = size(X_test)[1]
	n_trained = size(X_trained)[1]
	m_test = size(X_test)[2]
	m_trained = size(X_trained)[2]
	q_test = try
				size(COVARIATE_X_test)[2]
			catch
				0
			end
	q_trained = try
				size(COVARIATE_X_trained)[2]
			catch
				0
			end
	n_loci_with_effects = length(loci_info)
	### prepare X's and betas
	if (COVARIATE_X_trained == nothing) & (COVARIATE_X_test == nothing)
		INT_ANDOR_COVAR_trained = ones(n_trained)
		INT_ANDOR_COVAR_test = ones(n_test)
	else
		INT_ANDOR_COVAR_trained = hcat(ones(n_trained), COVARIATE_X_trained[:, 1:minimum([q_test, q_trained])])
		INT_ANDOR_COVAR_test = hcat(ones(n_test), COVARIATE_X_test[:, 1:minimum([q_test, q_trained])])
	end
	X_TRAINED = hcat(INT_ANDOR_COVAR_trained, X_trained[:, loci_info])
	X_TEST = hcat(INT_ANDOR_COVAR_test, X_test[:, loci_info])
	if int_and_covariate_effects == nothing
		X_TRAINED = X_TRAINED[:, 2:end]
		X_TEST = X_TEST[:, 2:end]
		BETA_TRAINED = beta_trained
	else
		BETA_TRAINED = vcat(int_and_covariate_effects[1:minimum([q_test, q_trained])+1], beta_trained) ### covariate means can be just the intercept mean across loci or the means of the intercept + covariates across loci
	end
	### prepare predict sum of products model
	y_trained_pred = X_TRAINED * BETA_TRAINED
	X_regress = convert(Array{Float64}, hcat(ones(n_trained), y_trained_pred)) #converting to a proper float type for the inv() function
	beta_regressor = inv(X_regress' * X_regress) * (X_regress' * y_trained)
	### predict
	y_pred_immature = X_TEST * BETA_TRAINED
	Y_PRED = convert(Array{Float64,1}, hcat(ones(n_test), y_pred_immature) * beta_regressor)
	return(Y_PRED)
end

### genomic prediction non-iterative models
function predict_NON_ITERATIVE_MODELS_func(X_test, loci_info, intcovar_eff, beta_trained; COVARIATE_X_test=nothing)
	# # ######################
	# # #test:
	# X_test=X_test
	# loci_info=TRAINING_OUT[i].LOCI
	# intcovar_eff=TRAINING_OUT[i].INTCOVAR_EFF
	# beta_trained=TRAINING_OUT[i].EFF
	# COVARIATE_X_test=COVARIATE_X_test
	# # ######################
	### data specs
	n = size(X_test)[1]
	m = length(loci_info)
	q_test = try
				size(COVARIATE_X_test)[2]
			catch
				0
			end
	q_trained = try
					length(intcovar_eff)-1
				catch
					0
				end
	### prepare X's and betas
	if (COVARIATE_X_test == nothing)
		INT_ANDOR_COVAR_test = ones(n)
	else
		INT_ANDOR_COVAR_test = hcat(ones(n), COVARIATE_X_test[:, 1:minimum([q_test, q_trained])])
	end
	X_TEST = hcat(INT_ANDOR_COVAR_test, X_test[:, loci_info])
	BETA_TRAINED = vcat(intcovar_eff[1:minimum([q_test, q_trained])+1], beta_trained) ### covariate means can be just the intercept mean across loci or the means of the intercept + covariates across loci
	### predict
	Y_PRED = X_TEST * BETA_TRAINED
	return(Y_PRED)
end

### assess genomic prediction accuracy
function asses_model_accuracy_func(y_true, y_pred; PLOT=false, PLOT_ID="")
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
		if PLOT==true
			Plots.plot([0,100], [0,100], seriestype=:line, color=:gray, xlab="Predicted", ylab="Observed") #1:1 line: perfection
			Plots.plot!(y_pred, y_true, seriestype=:scatter)
			newx = convert(Array{Float64}, collect(0:1:100))
			newy = INTERCEPT .+ (SLOPE .* newx)
			Plots.plot!([0,newx], [0,newy], seriestype=:line)
			Plots.savefig(string("Genomic_precition_accuracy_scatterplot_", PLOT_ID, ".png"))
		end
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

### extract QTL information ### NOTE: ADD QTL FREQUENCIES!
function QTL_SPEC_func(QTL_spec_fname)
	### getting the cummulative allele effect contribution
	qtl_spec = sort(CSV.read(QTL_spec_fname), (:CHROM, :POS))
	nqtl = length(unique(string.(qtl_spec.CHROM, qtl_spec.POS)))
	nalleles = convert(Int, nrow(qtl_spec) / nqtl)
	allele_frac_eff = qtl_spec.EFFECT ./ sum(qtl_spec.EFFECT)
	qtl_frac_eff = [sum(allele_frac_eff[(((i-1)*nalleles)+1):(i*nalleles)]) for i in 1:nqtl]
	qtl_max_eff = [maximum(allele_frac_eff[(((i-1)*nalleles)+1):(i*nalleles)]) for i in 1:nqtl]
	### QTL_SPEC dataframe
	QTL_SPEC = hcat(qtl_spec[qtl_spec.ALLELE .== "A", 1:2], qtl_frac_eff)
	DataFrames.rename!(QTL_SPEC, :x1 => :EFF_FRAC)
	QTL_SPEC.CHROM = [string(x) for x in QTL_SPEC.CHROM]
	return(QTL_SPEC)
end

### QTL detection power
function QTL_detection_power_func(TRAINING_OUT_i, id_merge, QTL_SPEC, LD_kb, pval_adj_method, alpha; PLOT=false)
	# BONFERRONI_THRESHOLD = -log10(0.05/length(TRAINING_OUT_i.LOCI))
	### p-value adjustments
	if (pval_adj_method == "Bonferroni")
		PVAL_adj = MultipleTesting.adjust(MultipleTesting.PValues(convert(Array{Float64}, TRAINING_OUT_i.PVAL)), Bonferroni())
	elseif (pval_adj_method == "BenjaminiHochberg")
		PVAL_adj = MultipleTesting.adjust(MultipleTesting.PValues(convert(Array{Float64}, TRAINING_OUT_i.PVAL)), BenjaminiHochberg())
	elseif (pval_adj_method == "BenjaminiYekutieli")
		PVAL_adj = MultipleTesting.adjust(MultipleTesting.PValues(convert(Array{Float64}, TRAINING_OUT_i.PVAL)), BenjaminiYekutieli())
	else
		println(string(pval_adj_method, " is not a valid p-value adjustment method. Using the Bonferroni method."))
		PVAL_adj = MultipleTesting.adjust(MultipleTesting.PValues(convert(Array{Float64}, TRAINING_OUT_i.PVAL)), Bonferroni())
	end
	LOD = -log.(10, PVAL_adj .+ 1e-20) ### adding a 1e-20 to avoid infinities
	# PUTATIVE_QTL = id_merge[TRAINING_OUT_i.LOCI[TRAINING_OUT_i.LOD .>= BONFERRONI_THRESHOLD], :]
	PUTATIVE_QTL = id_merge[TRAINING_OUT_i.LOCI[LOD .>= -log(10, alpha)], :]
	PUTATIVE_QTL.CHROM = Array{String}(PUTATIVE_QTL.CHROM)
	QTL_DETECTED_CHROM = []
	QTL_DETECTED_POS = []
	QTL_DETECTED_FREF = [] #QTL fraction effect
	for q_true in 1:size(QTL_SPEC)[1]
		for q_putative in 1:size(PUTATIVE_QTL)[1]
			if (QTL_SPEC.CHROM[q_true] == PUTATIVE_QTL.CHROM[q_putative]) & (abs(QTL_SPEC.POS[q_true] - PUTATIVE_QTL.POS[q_putative]) <= (LD_kb * 1000))
				push!(QTL_DETECTED_CHROM, PUTATIVE_QTL.CHROM[q_putative])
				push!(QTL_DETECTED_POS, PUTATIVE_QTL.POS[q_putative])
				push!(QTL_DETECTED_FREF, QTL_SPEC.EFF_FRAC[q_true])	# actually contain multple causal loci in LD
			end
		end
	end
	qtl_detected_id = unique(string.(QTL_DETECTED_CHROM, "_", QTL_DETECTED_POS))
	### Plot
	if PLOT == true
		nLoci_x_nAlleles = nrow(id_merge)
		idx_continuous = collect(1:nLoci_x_nAlleles)
		LOD_x = idx_continuous[TRAINING_OUT_i.LOCI]
		LOD_unadj = TRAINING_OUT_i.LOD
		LOD_adj = LOD
		QTL_x = []
		for i in 1:nrow(QTL_SPEC)
			append!(QTL_x, idx_continuous[(id_merge.CHROM .== QTL_SPEC.CHROM[i]) .& (id_merge.POS .== QTL_SPEC.POS[i])])
		end
		colours = ColorBrewer.palette("Set1", 9)
		Plots.plot([1, nLoci_x_nAlleles], [0, maximum(LOD_unadj)],
			seriestype=:scatter,
			marker=(:circle, 5, 0.5, :White, Plots.stroke(0, :white)),
			xlabel="SNP ID",
			ylabel="LOD",
			size=(1700, 500),
			lab=nothing)
		Plots.vline!(QTL_x,
			line=(7, :solid, 0.05, colours[5]),
			lab="QTL locations")
		Plots.plot!(LOD_x, LOD_unadj,
			seriestype=:scatter,
			marker=(:circle, 7, 0.7, colours[3], Plots.stroke(0, :white)),
			lab="LOD unadjusted")
		Plots.plot!(LOD_x, LOD_adj,
			seriestype=:scatter,
			marker=(:circle, 7, 0.7, colours[1], Plots.stroke(0, :white)),
			lab="LOD adjusted")
			Plots.savefig("Plot.png")
	end
	return(PUTATIVE_QTL, QTL_DETECTED_CHROM, QTL_DETECTED_POS, QTL_DETECTED_FREF, qtl_detected_id)
end

### within/across population + individual genotype data training and validataion
function CROSSVAL_INDI_func(POP_FNAME_INDI_DF, MODELS_INDI_DF, QTL_SPEC, pval_adj_method, LD_kb, alpha)
	# ### list of supported base models (base as in we can suffix these with _PC and _K, referring to PC and K covariates, respectively)
	# LIST_BASE_ITERATIVE_MODELS = ["GWAS"]
	# LIST_BASE_NON_ITERATIVE_MODELS = ["LS", "RR", "LASSO", "GLMNET"]
	### ouput arrays
	POP_TRAIN = []
	POP_TEST = []
	MODEL_ITERATION = []
	MODEL_COVARIATE = []
	MODEL_MODEL = []
	PREDICTORS = []
	NON_ZERO_PREDICTORS = []
	MEAN_DEVIANCE = []
	VAR_DEVIANCE = []
	CORRELATION = []
	INTERCEPT = []
	SLOPE = []
	R2 = []
	RMSD = []
	QTL_DETECTED_ID = []
	QTL_DETECTED_PERC = []
	QTL_FRAC_EFF = []
	FALSE_POSITIVE_RATE = []
	QTL_FREQS_TRAINING = []
	QTL_FREQS_VALIDATION = []
	### iterate across populations
	counter = [0]
	progress_bar = ProgressMeter.Progress((nrow(MODELS_INDI_DF) * length(POP_FNAME_INDI_DF.POP)^2), dt=1, barglyphs=BarGlyphs("[=> ]"), color=:yellow)
	for pop_train in POP_FNAME_INDI_DF.POP
		# pop_train="p001"
		geno_train_fname = POP_FNAME_INDI_DF.GENO[POP_FNAME_INDI_DF.POP .== pop_train][1]
		GENO_train = DelimitedFiles.readdlm(geno_train_fname, ',')
		id_merge = DataFrames.DataFrame(CHROM=string.(GENO_train[:,1]), POS=convert(Array{Int64}, GENO_train[:,2]), ALLELE=convert(Array{String}, GENO_train[:,3])) ###WE'RE ASSUMING SAME SETS OF LOCI FOR BOTH TRAINING AND VALIDATION POPULATIONS! 20190821
		X_train = convert(Array{Float64}, GENO_train[:, 4:end]')
		pheno_train_fname = POP_FNAME_INDI_DF.PHENO[POP_FNAME_INDI_DF.POP .== pop_train][1]
		pheno_train = CSV.read(pheno_train_fname, delim=",")
		y_train = convert(Array{Float64}, pheno_train.y)
		### QTL allele frequencies
		QTL_idx =[]
		for i in 1:nrow(QTL_SPEC)
			append!(QTL_idx, collect(1:nrow(id_merge))[(QTL_SPEC.CHROM[i] .== id_merge.CHROM) .& (QTL_SPEC.POS[i] .== id_merge.POS)])
		end
		QTL_FREQS_train = mean(X_train[:, QTL_idx], dims=1) ./ 2.0
		### Covariate estimations
		PC_train =  DelimitedFiles.readdlm(POP_FNAME_INDI_DF.COVARIATE_PC[POP_FNAME_INDI_DF.POP .== pop_train][1], ',')
		K_train =  DelimitedFiles.readdlm(POP_FNAME_INDI_DF.COVARIATE_K[POP_FNAME_INDI_DF.POP .== pop_train][1], ',')
		##################################
		### TRAINING AND QTL DETECTION ###
		##################################
		TRAINING_OUT = []
		for i in 1:nrow(MODELS_INDI_DF)
			# i = 6
			println(i)
			model = MODELS_INDI_DF.NAME[i] #e.g. model="ITERATIVE_NONE_FIXED_LS"
			# println(model)
			iteration_group = split(model, '_')[1]
			covariate_group = split(model, '_')[2]
			if covariate_group == "PC"
				COVARIATE_train = copy(PC_train)
			elseif covariate_group == "K"
				COVARIATE_train = copy(K_train)
			else
				COVARIATE_train = nothing
			end
			fmixedmodel_group = string(split(model, '_')[3], "_", split(model, '_')[4])
			# println(fmixedmodel_group)
			# nC = 10 ### NOTE: number of covariate columns (PC and K) to force for computational efficiency
			if iteration_group == "ITERATIVE"
				# if (fmixedmodel_group != "MIXED_EMMAX") & (COVARIATE_train != nothing)
				# 	_temp_ = COVARIATE_train[:,1:nC]
			    #     COVARIATE_train = copy(_temp_)
				# end
				LOCI, INTCOVAR_EFF, EFF, PVAL, LOD = GWAS_module.GWAS(X_train, y_train, 1.0/(2*length(y_train)); COVARIATE=COVARIATE_train, MODEL_TYPE=fmixedmodel_group)
			elseif iteration_group == "NON-ITERATIVE"
				# if (split(fmixedmodel_group, '_')[1] != "MIXED") & (COVARIATE_train != nothing)
				# 	_temp_ = COVARIATE_train[:,1:nC]
			    #     COVARIATE_train = copy(_temp_)
				# end
				LOCI, INTCOVAR_EFF, EFF, PVAL, LOD = GP_module.GP(X_train, y_train, 1.0/(2*length(y_train)); COVARIATE=COVARIATE_train, MODEL=fmixedmodel_group)
			else
				println(string("AwWwwW! SowWwWyyyy! We have not implemented the model: ", model, " yet."))
				println("¯\\_(๑❛ᴗ❛๑)_/¯ ʚ(´◡`)ɞ")
			end
			FIT = 	(MODEL_ITERATION=iteration_group, MODEL_COVARIATE=covariate_group, MODEL_MODEL=fmixedmodel_group,
					 LOCI=LOCI, INTCOVAR_EFF=INTCOVAR_EFF, EFF=EFF, PVAL=PVAL, LOD=LOD, COVARIATE=COVARIATE_train)
			push!(TRAINING_OUT, FIT)
		end
		for pop_test in POP_FNAME_INDI_DF.POP
			# pop_test="p002"
			geno_test_fname = POP_FNAME_INDI_DF.GENO[POP_FNAME_INDI_DF.POP .== pop_test][1]
			pheno_test_fname = POP_FNAME_INDI_DF.PHENO[POP_FNAME_INDI_DF.POP .== pop_test][1]
			# load 0,1,2 genotype matrices of pop1 and pop2
			# GENO_train = DelimitedFiles.readdlm(geno_train_fname, ',')
			GENO_test = DelimitedFiles.readdlm(geno_test_fname, ',')
			### merge the two genotype data so we start off with common SNPs and alleles
			# id_merge, X_train, X_test = merge_two_geno_func(GENO_train, GENO_test, "csv") ### EXCLUDED: WE'RE ASSUMING SAME SETS OF LOCI FOR BOTH TRAINING AND VALIDATION POPULATIONS! 20190821
			X_test = convert(Array{Float64}, GENO_test[:, 4:end]')
			# pheno_train = CSV.read(pheno_train_fname, delim=",")
			pheno_test = CSV.read(pheno_test_fname, delim=",")
			### preparte the phenotype data
			# y_train = convert(Array{Float64}, pheno_train.y)
			y_test = convert(Array{Float64}, pheno_test.y)
			### QTL allele frequencies
			QTL_idx =[]
			for i in 1:nrow(QTL_SPEC)
				append!(QTL_idx, collect(1:nrow(id_merge))[(QTL_SPEC.CHROM[i] .== id_merge.CHROM) .& (QTL_SPEC.POS[i] .== id_merge.POS)])
			end
			QTL_FREQS_test = mean(X_test[:, QTL_idx], dims=1) ./ 2.0
			### Covariate estimations
			PC_test = DelimitedFiles.readdlm(POP_FNAME_INDI_DF.COVARIATE_PC[POP_FNAME_INDI_DF.POP .== pop_test][1], ',')
			K_test = DelimitedFiles.readdlm(POP_FNAME_INDI_DF.COVARIATE_K[POP_FNAME_INDI_DF.POP .== pop_test][1], ',')
			######################################
			### PREDICTION AND QTL VALIUDATION ###
			######################################
			for i in 1:size(TRAINING_OUT, 1)
				# i = 1
				# println(i)
				model_iteration = TRAINING_OUT[i].MODEL_ITERATION
				model_covariate = TRAINING_OUT[i].MODEL_COVARIATE
				model_model = TRAINING_OUT[i].MODEL_MODEL
				##############################
				### EXTRACT THE COVARIATES ###
				##############################
				if TRAINING_OUT[i].COVARIATE == nothing
					# COVARIATE_X_trained = nothing
					COVARIATE_X_test = nothing
				elseif model_covariate =="PC"
					COVARIATE_X_test = PC_test
				elseif model_covariate =="K"
					COVARIATE_X_test = K_test
				else
					println(string("AwWwwW! SowWwWyyyy! We have not implemented using the covariate: ", model_covariate, " yet."))
					println("¯\\_(๑❛ᴗ❛๑)_/¯ ʚ(´◡`)ɞ")
				end
				if model_iteration == "ITERATIVE"
					########################
					### ITERATIVE MODELS ###
					########################
					y_pred = predict_ITERATIVE_MODELS_func(X_test, TRAINING_OUT[i].LOCI, TRAINING_OUT[i].INTCOVAR_EFF, TRAINING_OUT[i].EFF, y_train, X_train; COVARIATE_X_trained=TRAINING_OUT[i].COVARIATE, COVARIATE_X_test=COVARIATE_X_test)
				elseif model_iteration == "NON-ITERATIVE"
					############################
					### NON-ITERATIVE MODELS ###
					############################
					y_pred = predict_NON_ITERATIVE_MODELS_func(X_test, TRAINING_OUT[i].LOCI, TRAINING_OUT[i].INTCOVAR_EFF, TRAINING_OUT[i].EFF; COVARIATE_X_test=COVARIATE_X_test)
				else
					println(string("AwWwwW! SowWwWyyyy! We have not implemented the ", model_iteration, " model yet."))
					println("¯\\_(๑❛ᴗ❛๑)_/¯ ʚ(´◡`)ɞ")
				end
				### assess acuracy of the genomic prediction model
				accuracy = asses_model_accuracy_func(y_test, y_pred)
				push!(POP_TRAIN, pop_train)
				push!(POP_TEST, pop_test)
				push!(MODEL_ITERATION, model_iteration)
				push!(MODEL_COVARIATE, model_covariate)
				push!(MODEL_MODEL, model_model)
				push!(PREDICTORS, length(TRAINING_OUT[i].EFF))
				push!(NON_ZERO_PREDICTORS, sum(TRAINING_OUT[i].EFF .!= 0.0))
				push!(MEAN_DEVIANCE, accuracy.MEAN_DEVIANCE[1])
				push!(VAR_DEVIANCE, accuracy.VAR_DEVIANCE[1])
				push!(CORRELATION, accuracy.CORRELATION[1])
				push!(INTERCEPT, accuracy.INTERCEPT[1])
				push!(SLOPE, accuracy.SLOPE[1])
				push!(R2, accuracy.R2[1])
				push!(RMSD, accuracy.RMSD[1])
				### Power to detect QTL ### NOTE: NEED TO TRANSFORM MIXED MODEL LOD SCORES TO BONFERRONI-THESHOLDA-ABLE!!! 20190912
				PUTATIVE_QTL, QTL_DETECTED_CHROM, QTL_DETECTED_POS, QTL_DETECTED_FREF, qtl_detected_id = QTL_detection_power_func(TRAINING_OUT[i], id_merge, QTL_SPEC, LD_kb, pval_adj_method, alpha)
				push!(QTL_DETECTED_ID, join(qtl_detected_id, ';'))	#list of QTL chrom_position separated by ';'
				push!(QTL_DETECTED_PERC, length(qtl_detected_id) / size(QTL_SPEC)[1])
				if length(qtl_detected_id) > 0
					push!(QTL_FRAC_EFF, sum(QTL_DETECTED_FREF) / length(qtl_detected_id))
				else
					push!(QTL_FRAC_EFF, 0.0)
				end
				push!(FALSE_POSITIVE_RATE, (size(PUTATIVE_QTL, 1) - length(QTL_DETECTED_CHROM)) / size(PUTATIVE_QTL, 1))
				push!(QTL_FREQS_TRAINING, join(QTL_FREQS_train, ';'))
				push!(QTL_FREQS_VALIDATION, join(QTL_FREQS_test, ';'))
				### update progress bar
				counter[1] = counter[1] + 1
				ProgressMeter.update!(progress_bar, counter[1])
			end
		end
	end
	println("Summarising output.")
	OUT = DataFrames.DataFrame( POP_TRAIN=POP_TRAIN, POP_TEST=POP_TEST, MODEL_ITERATION=MODEL_ITERATION, MODEL_COVARIATE=MODEL_COVARIATE, MODEL_MODEL=MODEL_MODEL,
								PREDICTORS=PREDICTORS, NON_ZERO_PREDICTORS=NON_ZERO_PREDICTORS, MEAN_DEVIANCE=MEAN_DEVIANCE, VAR_DEVIANCE=VAR_DEVIANCE, CORRELATION=CORRELATION, INTERCEPT=INTERCEPT, SLOPE=SLOPE, R2=R2, RMSD=RMSD,
								QTL_DETECTED_PERC=QTL_DETECTED_PERC, QTL_FRAC_EFF=QTL_FRAC_EFF, QTL_DETECTED_ID=QTL_DETECTED_ID, FALSE_POSITIVE_RATE=FALSE_POSITIVE_RATE, QTL_FREQS_TRAINING=QTL_FREQS_TRAINING, QTL_FREQS_VALIDATION=QTL_FREQS_VALIDATION)
	return(OUT)
end

### within/across population + pool genotype data training and validataion
function CROSSVAL_POOL_func(POP_FNAME_POOL_DF, MODELS_POOL_DF, QTL_SPEC, pval_adj_method, LD_kb, alpha)
	### ouput arrays
	POP_TRAIN = []
	POP_TEST = []
	MODEL_ITERATION = []
	MODEL_COVARIATE = []
	MODEL_MODEL = []
	PREDICTORS = []
	NON_ZERO_PREDICTORS = []
	MEAN_DEVIANCE = []
	VAR_DEVIANCE = []
	CORRELATION = []
	INTERCEPT = []
	SLOPE = []
	R2 = []
	RMSD = []
	QTL_DETECTED_ID = []
	QTL_DETECTED_PERC = []
	QTL_FRAC_EFF = []
	FALSE_POSITIVE_RATE = []
	QTL_FREQS_TRAINING = []
	QTL_FREQS_VALIDATION = []
	### iterate across populations
	counter = [0]
	progress_bar = ProgressMeter.Progress((nrow(MODELS_POOL_DF) * length(POP_FNAME_POOL_DF.POP)^2), dt=1, barglyphs=BarGlyphs("[=> ]"), color=:yellow)
	for pop_train in POP_FNAME_POOL_DF.POP
		# pop_train="p02"
		fname_geno_sync_train = POP_FNAME_POOL_DF.GENO_SYNC[POP_FNAME_POOL_DF.POP .== pop_train][1]
		fname_geno_csv_train = POP_FNAME_POOL_DF.GENO_CSV[POP_FNAME_POOL_DF.POP .== pop_train][1]
		fname_pheno_py_train = POP_FNAME_POOL_DF.PHENO_PY[POP_FNAME_POOL_DF.POP .== pop_train][1]
		fname_pheno_csv_train = POP_FNAME_POOL_DF.PHENO_CSV[POP_FNAME_POOL_DF.POP .== pop_train][1]
		X_train = try
			DelimitedFiles.readdlm(fname_geno_csv_train, ',')[:, 4:end]'
		catch
			rm(fname_geno_csv_train, force=true)
			sync_parsing_module.sync_parse(fname_geno_sync_train) #writes out this output allele frequency csv file: string(split(geno_train_fname, ".")[end-1], "_ALLELEFREQ.csv")
			DelimitedFiles.readdlm(fname_geno_csv_train, ',')[:, 4:end]'
		end
		X_sync_train  = DelimitedFiles.readdlm(fname_geno_sync_train, '\t')
		y_df_train = CSV.read(fname_pheno_csv_train, datarow=1)
		y_train = convert(Array{Float64}, y_df_train[:,2])
		nLoci = size(X_train, 2)
		nPools = nrow(y_df_train)
		poolSizes = y_df_train[:,1]
		MAF = 1.0/(2 * sum(poolSizes .* nPools))
		LOCI = collect(1:nLoci)[ filter_sync_module.filter_sync(filename_sync=fname_geno_sync_train, MAF=MAF, DEPTH=1) ]
		id_merge = DataFrames.DataFrame(CHROM=repeat(string.(X_sync_train[:,1]), inner=6), POS=repeat(convert(Array{Int64}, X_sync_train[:,2]), inner=6), ALLELE=repeat(["A", "T", "C", "G", "N", "DEL"], outer=size(X_sync_train,1))) ###WE'RE ASSUMING SAME SETS OF LOCI FOR BOTH TRAINING AND VALIDATION POPULATIONS! 20190821
		### QTL allele frequencies
		QTL_idx =[]
		for i in 1:nrow(QTL_SPEC)
			append!(QTL_idx, collect(1:nrow(id_merge))[(QTL_SPEC.CHROM[i] .== id_merge.CHROM) .& (QTL_SPEC.POS[i] .== id_merge.POS)])
		end
		QTL_FREQS_train = mean(X_train[:, QTL_idx], dims=1)
		### COVARIATES: Fst vs summary_statistics
		NPSTAT_train = DelimitedFiles.readdlm(POP_FNAME_POOL_DF.COVARIATE_NPSTAT[POP_FNAME_POOL_DF.POP .== pop_train][1], ',')
		FST_train = DelimitedFiles.readdlm(POP_FNAME_POOL_DF.COVARIATE_FST[POP_FNAME_POOL_DF.POP .== pop_train][1], ',')
		### QTL detection and predictor estimation
		println("Model training. Population: ", pop_train)
		TRAINING_OUT = []
		for i in 1:nrow(MODELS_POOL_DF)
			# i = 1
			println(i)
			model = MODELS_POOL_DF.NAME[i] #e.g. model="ITERATIVE_NONE_GWAlpha" and model="NON-ITERATIVE_NPSTAT_MIXED_GLMNET"
			# println(model)
			iteration_group = split(model, '_')[1]
			covariate_group = split(model, '_')[2]
			if covariate_group == "NPSTAT"
				COVARIATE_train = copy(NPSTAT_train)
			elseif covariate_group == "FST"
				COVARIATE_train = copy(FST_train)
			else
				COVARIATE_train = nothing
			end
			fmixedmodel_group = string(split(model, '_')[3], "_", split(model, '_')[4])
			if fmixedmodel_group == "FIXED_GWAlpha"
				OUT, COVAR_EFF = PoolGPAS_module.PoolGPAS(fname_geno_sync_train, fname_pheno_py_train, MAF, 1; MODEL=fmixedmodel_group, COVARIATE=COVARIATE_train)
				LOCI = OUT.LOCUS_ID
				INTCOVAR_EFF = nothing
				EFF = OUT.ALPHA
				PVAL = OUT.PVALUES
				LOD = OUT.LOD
			else
				OUT, COVAR_EFF = PoolGPAS_module.PoolGPAS(fname_geno_sync_train, fname_pheno_csv_train, MAF, 1; MODEL=fmixedmodel_group, COVARIATE=COVARIATE_train)
				LOCI = OUT.LOCUS_ID[2:end]
				COVAR_EFF == nothing ? INTCOVAR_EFF = [OUT.BETA[1]] : INTCOVAR_EFF = vcat(OUT.BETA[1], COVAR_EFF)
				EFF = OUT.BETA[2:end]
				PVAL = OUT.PVALUES[2:end]
				LOD = OUT.LOD[2:end]
			end
			FIT = 	(MODEL_ITERATION=iteration_group, MODEL_COVARIATE=covariate_group, MODEL_MODEL=fmixedmodel_group,
					 LOCI=LOCI, INTCOVAR_EFF=INTCOVAR_EFF, EFF=EFF, PVAL=PVAL, LOD=LOD, COVARIATE=COVARIATE_train)
			push!(TRAINING_OUT, FIT)
		end
		for pop_test in POP_FNAME_POOL_DF.POP
			println("Model validation. Population: ", pop_test)
			# pop_test="p03"
			fname_geno_sync_test = POP_FNAME_POOL_DF.GENO_SYNC[POP_FNAME_POOL_DF.POP .== pop_test][1]
			fname_geno_csv_test = POP_FNAME_POOL_DF.GENO_CSV[POP_FNAME_POOL_DF.POP .== pop_test][1]
			fname_pheno_py_test = POP_FNAME_POOL_DF.PHENO_PY[POP_FNAME_POOL_DF.POP .== pop_test][1]
			fname_pheno_csv_test = POP_FNAME_POOL_DF.PHENO_CSV[POP_FNAME_POOL_DF.POP .== pop_test][1]
			X_test = try
				DelimitedFiles.readdlm(fname_geno_csv_test, ',')[:, 4:end]'
			catch
				rm(fname_geno_csv_test, force=true)
				sync_parsing_module.sync_parse(fname_geno_sync_test) #writes out this output allele frequency csv file: string(split(geno_test_fname, ".")[end-1], "_ALLELEFREQ.csv")
				DelimitedFiles.readdlm(fname_geno_csv_test, ',')[:, 4:end]'
			end
			X_sync_test  = DelimitedFiles.readdlm(fname_geno_sync_test, '\t')
			y_df_test = CSV.read(fname_pheno_csv_test, datarow=1)
			y_test = convert(Array{Float64}, y_df_test[:,2])
			nLoci = size(X_test, 2)
			nPools = nrow(y_df_test)
			poolSizes = y_df_test[:,1]
			MAF = 1.0/(2 * sum(poolSizes .* nPools))
			LOCI = collect(1:nLoci)[ filter_sync_module.filter_sync(filename_sync=fname_geno_sync_test, MAF=MAF, DEPTH=1) ]
			id_merge = DataFrames.DataFrame(CHROM=repeat(string.(X_sync_test[:,1]), inner=6), POS=repeat(convert(Array{Int64}, X_sync_test[:,2]), inner=6), ALLELE=repeat(["A", "T", "C", "G", "N", "DEL"], outer=size(X_sync_test,1))) ###WE'RE ASSUMING SAME SETS OF LOCI FOR BOTH TRAINING AND VALIDATION POPULATIONS! 20190821
			### QTL allele frequencies
			QTL_idx =[]
			for i in 1:nrow(QTL_SPEC)
				append!(QTL_idx, collect(1:nrow(id_merge))[(QTL_SPEC.CHROM[i] .== id_merge.CHROM) .& (QTL_SPEC.POS[i] .== id_merge.POS)])
			end
			QTL_FREQS_test = mean(X_test[:, QTL_idx], dims=1)
			### COVARIATES: Fst vs summary_statistics
			NPSTAT_test = DelimitedFiles.readdlm(POP_FNAME_POOL_DF.COVARIATE_NPSTAT[POP_FNAME_POOL_DF.POP .== pop_test][1], ',')
			FST_test = DelimitedFiles.readdlm(POP_FNAME_POOL_DF.COVARIATE_FST[POP_FNAME_POOL_DF.POP .== pop_test][1], ',')
			### Prediction
			for i in 1:nrow(MODELS_POOL_DF)
				model_iteration = TRAINING_OUT[i].MODEL_ITERATION
				model_covariate = TRAINING_OUT[i].MODEL_COVARIATE
				model_model = TRAINING_OUT[i].MODEL_MODEL
				##############################
				### EXTRACT THE COVARIATES ###
				##############################
				if TRAINING_OUT[i].COVARIATE == nothing
					COVARIATE_X_trained = nothing
					COVARIATE_X_test = nothing
				elseif model_covariate =="NPSTAT"
					COVARIATE_X_test = NPSTAT_test
				elseif model_covariate =="FST"
					COVARIATE_X_test = FST_test
				else
					println(string("AwWwwW! SowWwWyyyy! We have not implemented using the covariate: ", model_covariate, " yet."))
					println("¯\\_(๑❛ᴗ❛๑)_/¯ ʚ(´◡`)ɞ")
				end
				if model_iteration == "ITERATIVE"
					########################
					### ITERATIVE MODELS ###
					########################
					y_pred = predict_ITERATIVE_MODELS_func(X_test, TRAINING_OUT[i].LOCI, TRAINING_OUT[i].INTCOVAR_EFF, TRAINING_OUT[i].EFF, y_train, X_train; COVARIATE_X_trained=TRAINING_OUT[i].COVARIATE, COVARIATE_X_test=COVARIATE_X_test)
				elseif model_iteration == "NON-ITERATIVE"
					############################
					### NON-ITERATIVE MODELS ###
					############################
					y_pred = predict_NON_ITERATIVE_MODELS_func(X_test, TRAINING_OUT[i].LOCI, TRAINING_OUT[i].INTCOVAR_EFF, TRAINING_OUT[i].EFF; COVARIATE_X_test=COVARIATE_X_test)
				else
					println(string("AwWwwW! SowWwWyyyy! We have not implemented the ", model_iteration, " model yet."))
					println("¯\\_(๑❛ᴗ❛๑)_/¯ ʚ(´◡`)ɞ")
				end
				### assess acuracy of the genomic prediction model
				accuracy = asses_model_accuracy_func(y_test, y_pred)
				push!(POP_TRAIN, pop_train)
				push!(POP_TEST, pop_test)
				push!(MODEL_ITERATION, model_iteration)
				push!(MODEL_COVARIATE, model_covariate)
				push!(MODEL_MODEL, model_model)
				push!(PREDICTORS, length(TRAINING_OUT[i].EFF))
				push!(NON_ZERO_PREDICTORS, sum(TRAINING_OUT[i].EFF .!= 0.0))
				push!(MEAN_DEVIANCE, accuracy.MEAN_DEVIANCE[1])
				push!(VAR_DEVIANCE, accuracy.VAR_DEVIANCE[1])
				push!(CORRELATION, accuracy.CORRELATION[1])
				push!(INTERCEPT, accuracy.INTERCEPT[1])
				push!(SLOPE, accuracy.SLOPE[1])
				push!(R2, accuracy.R2[1])
				push!(RMSD, accuracy.RMSD[1])
				### Power to detect QTL ### NOTE: NEED TO TRANSFORM MIXED MODEL LOD SCORES TO BONFERRONI-THESHOLDA-ABLE!!! 20190912
				PUTATIVE_QTL, QTL_DETECTED_CHROM, QTL_DETECTED_POS, QTL_DETECTED_FREF, qtl_detected_id = QTL_detection_power_func(TRAINING_OUT[i], id_merge, QTL_SPEC, LD_kb, pval_adj_method, alpha)
				push!(QTL_DETECTED_ID, join(qtl_detected_id, ';'))	#list of QTL chrom_position separated by ';'
				push!(QTL_DETECTED_PERC, length(qtl_detected_id) / size(QTL_SPEC)[1])
				if length(qtl_detected_id) > 0
					push!(QTL_FRAC_EFF, sum(QTL_DETECTED_FREF) / length(qtl_detected_id))
				else
					push!(QTL_FRAC_EFF, 0.0)
				end
				push!(FALSE_POSITIVE_RATE, (size(PUTATIVE_QTL, 1) - length(QTL_DETECTED_CHROM)) / size(PUTATIVE_QTL, 1))
				push!(QTL_FREQS_TRAINING, join(QTL_FREQS_train, ';'))
				push!(QTL_FREQS_VALIDATION, join(QTL_FREQS_test, ';'))
				### update progress bar
				counter[1] = counter[1] + 1
				ProgressMeter.update!(progress_bar, counter[1])
			end
		end
	end
	println("Summarising output.")
	OUT = DataFrames.DataFrame( POP_TRAIN=POP_TRAIN, POP_TEST=POP_TEST, MODEL_ITERATION=MODEL_ITERATION, MODEL_COVARIATE=MODEL_COVARIATE, MODEL_MODEL=MODEL_MODEL,
								PREDICTORS=PREDICTORS, NON_ZERO_PREDICTORS=NON_ZERO_PREDICTORS, MEAN_DEVIANCE=MEAN_DEVIANCE, VAR_DEVIANCE=VAR_DEVIANCE, CORRELATION=CORRELATION, INTERCEPT=INTERCEPT, SLOPE=SLOPE, R2=R2, RMSD=RMSD,
								QTL_DETECTED_PERC=QTL_DETECTED_PERC, QTL_FRAC_EFF=QTL_FRAC_EFF, QTL_DETECTED_ID=QTL_DETECTED_ID, FALSE_POSITIVE_RATE=FALSE_POSITIVE_RATE, QTL_FREQS_TRAINING=QTL_FREQS_TRAINING, QTL_FREQS_VALIDATION=QTL_FREQS_VALIDATION)
	return(OUT)
end

#####################
###				  ###
### main function ###
###				  ###
#####################

#############
### input ###
#############
println("Loading input files")
DIR=ARGS[1]
cd(DIR)
POP_FNAME_DF_filename = ARGS[2]
POP_FNAME_DF = CSV.read(POP_FNAME_DF_filename, datarow=2) 		### csv file of lists of individual or pool data per population or across populations; HEADER: POP(population ID), GENO(full path + filename of the genotype *.csv file), PHENO(full path + filename of phenotype *.csv file)
### HEADERS:
### POP,GENO,PHENO,COVARIATE_PC,COVARIATE_K for: POP_FNAME_WITHIN_INDI_DF.csv and POP_FNAME_ACROSS_INDI_DF.csv
### POP,GENO_SYNC,GENO_CSV,PHENO_PY,PHENO_CSV,COVARIATE_NPSTAT,COVARIATE_FST for: POP_FNAME_WITHIN_POOL_DF.csv or POP_FNAME_ACROSS_POOL_DF.csv
TEST = ARGS[3]											### cross-validations to perform: {"WITHIN_INDI", "WITHIN_POOL", "ACROSS_INDI", "ACROSS_POOL"}
QTL_spec_fname = ARGS[4] 								### QTL identities (CHROM,POS,ALLELE,EFFECT)
LD_kb = parse(Float64, ARGS[5])							### size of linkage block in kilobases
pval_adj_method = ARGS[6]								### p-value adjustment method: "Bonferroni", "BenjaminiHochberg" or "BenjaminiYekutieli"
alpha = parse(Float64, ARGS[7])							### significance level (decimal): prbability of incorrectly  rejecting the null hypothesis
# ### TESTS START #####################
# POP_FNAME_DF_filename = "POP_FNAME_WITHIN_INDI_DF.csv"
# 	# POP_FNAME_DF_filename = "POP_FNAME_ACROSS_INDI_DF.csv"
# 	# POP_FNAME_DF_filename = "POP_FNAME_WITHIN_POOL_DF.csv"
# 	# POP_FNAME_DF_filename = "POP_FNAME_ACROSS_POOL_DF.csv"
# POP_FNAME_DF = CSV.read(POP_FNAME_DF_filename, datarow=2)
# TEST = "INDI"
# 	# TEST = "POOL"
# QTL_spec_fname = "../QTL_SPEC.csv"
# LD_kb = 1
# pval_adj_method = "Bonferroni"
# alpha = 0.01
# ### TESTS END #######################

### GWAS & GP models: INDIVIDUAL GENOTPYE DATA
COL1_INDI = vcat(repeat(["ITERATIVE"], inner=7), repeat(["NON-ITERATIVE"], inner=18))
COL2_INDI = vcat(["NONE"], repeat(["PC", "K"], inner=3), repeat(["NONE"], inner=4), repeat(["PC", "K"], inner=7))
COL3_INDI = vcat(["FIXED"], repeat(["FIXED", "MIXED", "MIXED"], outer=2), repeat(["FIXED"], inner=4), repeat(vcat(repeat(["FIXED"], inner=4), repeat(["MIXED"], inner=3)), outer=2))
COL4_INDI = vcat(["LS"], repeat(["LS", "FAST", "EMMAX"], outer=2), ["LS", "RR", "GLMNET", "LASSO"], repeat(["LS", "RR", "GLMNET", "LASSO", "RR", "GLMNET", "LASSO"], outer=2))
MODEL_NAMES = string.(COL1_INDI, "_", COL2_INDI, "_", COL3_INDI, "_", COL4_INDI)
MODELS_INDI_DF = DataFrames.DataFrame(ALGORITHM=COL1_INDI, COVARIATE=COL2_INDI, FIXED_MIXED=COL3_INDI, MODEL=COL4_INDI, NAME=MODEL_NAMES)
### excluding highly computationally inefficient iterative mixed models (need better implementation of EMMAx)
MODELS_INDI_DF = MODELS_INDI_DF[(MODELS_INDI_DF.ALGORITHM .== "NON-ITERATIVE") .| ( (MODELS_INDI_DF.ALGORITHM .== "ITERATIVE") .& (MODELS_INDI_DF.COVARIATE .== "NONE") ), :]
### prohibiting PC covariates (1 PC) in MIXED MODELS (because they require a square Z matrix)
MODELS_INDI_DF = MODELS_INDI_DF[ (MODELS_INDI_DF.COVARIATE .!= "PC") .| (MODELS_INDI_DF.FIXED_MIXED .!= "MIXED"), : ]

### Pool-GWAS models: POOL-SEQ DATA
COL1_POOL = vcat(["ITERATIVE"], repeat(["NON-ITERATIVE"], inner=14))
COL2_POOL = vcat(repeat(["NONE"], inner=3), repeat(["NPSTAT", "FST"], inner=4), repeat(["NPSTAT", "FST"], inner=2))
COL3_POOL = vcat(repeat(["FIXED"], inner=11), repeat(["MIXED"], inner=4))
COL4_POOL = vcat(["GWAlpha", "LS", "RR"], repeat(["LS", "RR", "GLMNET", "LASSO"], outer=2), repeat(["GLMNET", "LASSO"], outer=2))
MODEL_NAMES = string.(COL1_POOL, "_", COL2_POOL, "_", COL3_POOL, "_", COL4_POOL)
MODELS_POOL_DF = DataFrames.DataFrame(ALGORITHM=COL1_POOL, COVARIATE=COL2_POOL, FIXED_MIXED=COL3_POOL, MODEL=COL4_POOL, NAME=MODEL_NAMES)

### extract the QTL information
QTL_SPEC = QTL_SPEC_func(QTL_spec_fname)

### execute cross-validations
if TEST == "INDI"
	##################################
	### WITHIN POPULATION TRAINING ### INDIVIDUAL GENOTYPE & PHENOTYPE DATA
	##################################
	println("Training and testing using individual genotype data")
	global OUT = CROSSVAL_INDI_func(POP_FNAME_DF, MODELS_INDI_DF, QTL_SPEC, pval_adj_method, LD_kb, alpha)
elseif TEST == "POOL"
	##################################
	### WITHIN POPULATION TRAINING ### POOL GENOTYPE & PHENOTYPE DATA
	##################################
	println("Training and testing using Pool-seq data")
	global OUT = CROSSVAL_POOL_func(POP_FNAME_DF, MODELS_POOL_DF, QTL_SPEC, pval_adj_method, LD_kb, alpha)
else
	println(string("Unrecognised cross-validation test: ", TEST, "."))
	println("Please select from: WIHIN_INDI, WITHIN_POOL, ACROSS_INDI, and ACROSS_POOL")
	exit()
end

### write the output
if dirname(POP_FNAME_DF_filename) == ""
	OUT_FNAME = string("CROSS_VALIDATION_OUTPUT_", POP_FNAME_DF_filename)
else
	OUT_FNAME = string(dirname(POP_FNAME_DF_filename), "/CROSS_VALIDATION_OUTPUT_", basename(POP_FNAME_DF_filename))
end

try
	CSV.write(OUT_FNAME, OUT)
catch
	### convert nothing into "NA"
	for i in 1:ncol(OUT)
		println(i)
		n_missing = sum(OUT[:, i] .== nothing)
		println(n_missing)
		OUT[(OUT[:, i] .== nothing), i] = repeat(["NA"], inner=n_missing)
	end
	CSV.write(OUT_FNAME, OUT)
end
println("====================================================")
println("Everything went well! The output file is:")
println(OUT_FNAME)
println("====================================================")

### SAMPLE EXECUTION:
# OUTDIR=/data/Lolium/Quantitative_Genetics/genomic_prediction_simulations/QUANTINEMO2_SIM/test_1kloci_2019-09-24_09-00-55
# cd ${OUTDIR}
# GEN_PRED_SRC_DIR=/data/Lolium/Softwares/genomic_prediction/src
# time \
# julia ${GEN_PRED_SRC_DIR}/QUANTINEMO_03_GWAS_GP.jl \
# 	${OUTDIR} \
# 	POP_FNAME_WITHIN_POOL_DF.csv \
# 	POOL \
# 	../QTL_SPEC.csv \
# 	1 \
# 	Bonferroni \
# 	0.01
