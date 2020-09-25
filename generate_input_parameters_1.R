# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to generate a matrix of epidemiological, cost, and utility inputs

generate_input_parameters <- function(n_samples, sensitivity="None") {
	n_parameters <- 44
	input_parameters <- matrix(NA, nrow = n_samples, ncol = n_parameters)

	parameter_names <- c("Log hazard Stroke", "Log hazard MI", "Log hazard Bleed", "Log hazard ICH",
	   "Log hazard Death", "Log hazard SE", "Log hazard TIA", "Log HR Apixaban Stroke", "Log HR Apixaban MI",
	   "Log HR Apixaban Bleed", "Log HR Apixaban ICH", "Log HR Apixaban Death", "Log HR Apixaban SE",
	   "Log HR Apixaban TIA", "Log HR Dabigatran Stroke", "Log HR Dabigatran MI", "Log HR Dabigatran Bleed",
	   "Log HR Dabigatran ICH", "Log HR Dabigatran Death", "Log HR Dabigatran SE", "Log HR Dabigatran TIA",
			"Stroke acute cost", "MI acute cost", "Bleed acute cost", "ICH acute cost", "SE acute cost",
			"TIA acute cost", "Coumarin cost", "Apixaban cost", "Dabigatran cost",
			"Stroke management cost", "MI management cost", "ICH management cost",
			"AF Well utility", "Stroke utility", "MI utility", "Bleed utility", "ICH utility",
			"Stroke disutility", "MI disutility", "Bleed disutility", "ICH disutility", "SE disutility",
			"TIA disutility")
	colnames(input_parameters) <- parameter_names

	# Load the baseline/coumarin log hazard and treatment log hazard ratios estimated by a 
	# network meta-analysis in OpenBUGS
	bugs_log_hr <- x<-as.matrix(read.csv(file = "bugs_loghr.csv"))
	bugs_log_hazard <- as.matrix(read.csv(file = "bugs_baseline.csv"))
	# Replace . with space to ease referencing
	# . is an escape character so needs double backslash in gsub
	colnames(bugs_log_hr) <- gsub("\\.", " ", colnames(bugs_log_hr))
	colnames(bugs_log_hazard) <- gsub("\\.", " ", colnames(bugs_log_hazard))
	
	# Load the baseline (coumarin) log hazards
	input_parameters[, paste("Log hazard", event_names)] <- bugs_log_hazard[1:n_samples, ]
	# Load the log hazard ratios
	# First create the vector of treatment and event names
	log_hr_names <- paste("Log HR", rep(treatment_names[-1], 
	                     each = length(event_names)),
	                     rep(event_names, 3))
	# Then extract the appropriate log hazard ratios
	input_parameters[, log_hr_names] <- bugs_log_hr[1:n_samples, log_hr_names]
	
	
	# SE acute costs
	input_parameters[, "SE acute cost"] <- runif(n_samples, 1186.5, 3559.5)
	# TIA acute costs
	input_parameters[, "TIA acute cost"] <- runif(n_samples, 532, 1596)
	# Bleed acute costs
	input_parameters[, "Bleed acute cost"] <- runif(n_samples, 875.75, 2627.25)
	# MI acute costs
	input_parameters[, "MI acute cost"] <- runif(n_samples, 2415.24, 7245.72)
	# Ischemic stroke acute costs
	S_acute_cost_mean = 11626
	S_acute_cost_SE = 16868 / sqrt(162)
	input_parameters[, "Stroke acute cost"] <- rnorm(n_samples, mean = S_acute_cost_mean, 
	                                                 sd = S_acute_cost_SE)
	# ICH acute costs
	I_acute_cost_mean=11453
	I_acute_cost_SE = 13815 / sqrt(17)
	input_parameters[, "ICH acute cost"] <- rnorm(n_samples, mean = I_acute_cost_mean, 
	                                              sd = I_acute_cost_SE)

  # Treatment costs are for 3 monthly cycles
	# Uniform distribution on Coumarin costs (these are mostly due to management costs)
	input_parameters[, "Coumarin cost"] <- runif(n_samples, 52.57, 157.70)
	# The DOAC costs are fixed
	input_parameters[, "Apixaban cost"] <- rep(200.42, n_samples)
	input_parameters[, "Dabigatran cost"] <- rep(200.42, n_samples)

	# Health state costs (divided by four to go from annual to quarterly costs)
	# Stroke
	S_cost_mean = 3613
	S_cost_SE = 4235 / sqrt(136)
	input_parameters[, "Stroke management cost"] <- rnorm(n_samples, mean = S_cost_mean, 
	                                                      sd = S_cost_SE) 
	# ICH (Assume it is similar to stroke; divided by four to go to quarterly costs)
	I_cost_mean = 3613
	I_cost_SE=4235 / sqrt(136)
	input_parameters[, "ICH management cost"] <- rnorm(n_samples, mean = I_cost_mean, sd = I_cost_SE)
	# MI adds only an instant cost, so this post-state has 0 management cost
	input_parameters[, "MI management cost"] <- rep(0, n_samples)

	# Health state utilities (from sources identified in Bayer Table 49)
	# These are combined proportionally
	# All utilities capped above at 1
	# All utilities divided by 4 to make them 3-monthly
	AF_baseline_utility <- rnorm(n_samples, 0.779, 0.0045)
	input_parameters[, "AF Well utility"] <- min(AF_baseline_utility, 1) / 4
	input_parameters[, "Stroke utility"] <- min(rnorm(n_samples, 0.69, 0.025), 1) / 4
	input_parameters[, "MI utility"] <- min(rnorm(n_samples, 0.718, 0.0163), 1) / 4
	input_parameters[, "ICH utility"] <- min(rbeta(n_samples, 3.941, 1.385), 1) / 4
	# Assume bleed utility the same as Stroke
	input_parameters[, "Bleed utility"] <- input_parameters[, "Stroke utility"]

	# Disutilities applied for only 3 months
	input_parameters[, "MI disutility"] <- (rnorm(n_samples, 0.683, 0.0156) - AF_baseline_utility) * 0.25
	input_parameters[,"Stroke disutility"] <- (runif(n_samples, min = 1.5 * (-0.59), max = 0.5 * (-0.59))) * 0.25
	input_parameters[,"TIA disutility"] <- (runif(n_samples, min= 1.5 * (-0.131), max = 0.5 * (-0.131))) * 0.25
	input_parameters[,"Bleed disutility"] <- - (rnorm(n_samples, 0.03, 0.001531)) * 0.25
	input_parameters[,"SE disutility"] <- (runif(n_samples, min = 1.5 * (-0.131), max = 0.5 * (-0.131))) * 0.25
	input_parameters[,"ICH disutility"] <- (rnorm(n_samples, 0.6, 0.064) - AF_baseline_utility) * 0.25

	return(input_parameters)
} # End function