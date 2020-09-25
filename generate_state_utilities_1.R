# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert input parameters to age dependent state utilities

generate_state_utilities <- function(input_parameters, age) {
  
	# Construct array of state utilities with default value 0
	state_utilities <- array(0, dim=c(n_samples, n_treatments, n_states), 
	                         dimnames = list(NULL, treatment_names, state_names))

	# Construct hazards of transient events
	baseline_hazard_SE <- exp(input_parameters[, "Log hazard SE"])
	baseline_hazard_TIA <- exp(input_parameters[, "Log hazard TIA"])
	
	
	# State utilities are a mixture of state event disutilities, state utilities, and disutilities of transient events
	
	for(i_treatment in 1:n_treatments) {
		# Add state utility (if any - death has none)
		for(i_state in 1:n_states) {
			# Check if there is a utility for this state
			if(sum(grepl(state_names[i_state], colnames(input_parameters)) & 
			       grepl("utility", colnames(input_parameters)) & 
			       !grepl("disutility", colnames(input_parameters))) > 0){
				# Then add the appropriate acute cost to the total state cost
				state_utilities[, i_treatment, i_state] <- state_utilities[, i_treatment, i_state] + 
				  input_parameters[, grepl(state_names[i_state], colnames(input_parameters)) & 
				                     grepl("utility", colnames(input_parameters)) & 
				                     !grepl("disutility", colnames(input_parameters))]
			}
		} # End loop over states

		# Now add state event disutility (if any)
		for(i_state in 1:n_states) {
			# Check if there is a disutility for this state
			if(sum(grepl(state_names[i_state], colnames(input_parameters)) & 
			       grepl("disutility", colnames(input_parameters))) > 0){
				# Then add the appropriate acute cost to the total state cost
				state_utilities[, i_treatment, i_state] <- state_utilities[, i_treatment,i_state] + 
				  input_parameters[, grepl(state_names[i_state], colnames(input_parameters)) & 
				                     grepl("disutility", colnames(input_parameters))]
			}
		} # End loop over states

		# Finally add the transient event disutilities (SE and TIA) which differ by treatment but not state
		# Each state has a proportion who experience transient events, use this to scale the disutility for these events
		# Only add hazard ratio if not Coumarin
		hazard_SE <- baseline_hazard_SE
		hazard_TIA <- baseline_hazard_TIA
		if(i_treatment != 1) hazard_SE <- baseline_hazard_SE * 
		  exp(input_parameters[, paste0("Log HR ",treatment_names[i_treatment], " SE")])
		# Calculate probabilities of SE and TIA
		prob_SE <- (1 - exp(-0.25 * hazard_SE))
		prob_TIA <- (1 - exp(-0.25 * hazard_TIA))
		# Now add disutility for SE and TIA (same applies to all states except death)
		state_utilities[, i_treatment, -n_states] <- state_utilities[, i_treatment, -n_states] + 
		  prob_SE * input_parameters[, "SE disutility"] + prob_TIA * input_parameters[, "TIA disutility"]
	} # End loop over treatments

	return(state_utilities)
} # End function