# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert input parameters to state costs

# TODO Use more efficient combination of grepl() to remove loops and if statements adding acute and management costs
# TODO Maybe fix probabilities of SE and TIA so they account for competing risks?

# TODO
# Change acute costs to be the same for all states
# They'll be a weighted average of events
# Hazards and hazard ratios are separately used to determine transition probabilities
# If using competing risks formula, this will solve competing risks issue above
# Do the same for utilities

generate_state_costs<-function(input_parameters) {
  
	# Construct array of state costs with default value of zero
	state_costs <- array(0, dim = c(n_samples, n_treatments, n_states), 
	                     dimnames = list(NULL, treatment_names, state_names))

	# Construct hazards of transient events
	baseline_hazard_SE <- exp(input_parameters[, "Log hazard SE"])
	baseline_hazard_TIA <- exp(input_parameters[, "Log hazard TIA"])
	
	# State costs are a mixture of treatment costs, acute event costs, event managment costs, 
	# and costs of transient events
	
	for(i_treatment in 1:n_treatments) {
		# First add the treatment costs (applies to all states except death)
		state_costs[, i_treatment, 1:(n_states - 1)] <-
		  input_parameters[, grep(paste(treatment_names[i_treatment], "cost"), 
		                          colnames(input_parameters))]

		# Add acute event costs (if any)
		for(i_state in 1:n_states) {
			# Check if there are acute costs for this state
			if(sum(grepl(state_names[i_state], colnames(input_parameters)) & 
			       grepl("acute cost", colnames(input_parameters))) > 0){
				# Then add the appropriate acute cost to the total state cost
				state_costs[, i_treatment, i_state] <- state_costs[, i_treatment, i_state] +
				  input_parameters[, grepl(state_names[i_state], colnames(input_parameters)) & 
				                     grepl("acute cost",colnames(input_parameters))]
			}
		} # End loop over states

		# Now add management costs (if any)
		for(i_state in 1:n_states){
			# Check if there are management costs for this state
			if(sum(grepl(state_names[i_state], colnames(input_parameters)) & 
			       grepl("management cost", colnames(input_parameters))) > 0){
				# Then add the appropriate management cost to the total state cost
				state_costs[, i_treatment, i_state] <- state_costs[, i_treatment, i_state] +
				  input_parameters[, grepl(state_names[i_state], colnames(input_parameters)) & 
				                     grepl("management cost", colnames(input_parameters))]
			}
		}

		# Finally add the transient event costs (SE and TIA) which differ by treatment but not state
		# Each state has a proportion who experience transient events, use this to scale the acute costs for these events
		# Only add hazard ratio if not Coumarin
		hazard_SE <- baseline_hazard_SE
		hazard_TIA <- baseline_hazard_TIA
		if(i_treatment != 1) hazard_SE <- baseline_hazard_SE * 
		  exp(input_parameters[, paste0("Log HR ", treatment_names[i_treatment], " SE")])
		# Calculate probabilities of SE and TIA
		prob_SE <- (1 - exp(-0.25 * hazard_SE))
		prob_TIA <- (1 - exp(-0.25 * hazard_TIA))
		# Now add acute costs for SE and TIA (same applies to all states except dead)
		state_costs[, i_treatment, -n_states] <- state_costs[, i_treatment, -n_states] + 
		  prob_SE * input_parameters[, "SE acute cost"] +
		    prob_TIA * input_parameters[, "TIA acute cost"]
	} # End loop over treatments

	return(state_costs)
} # End function
