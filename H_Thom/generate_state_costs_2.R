# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert input parameters to state costs


generate_state_costs<-function(input_parameters) {
  
	# Construct array of state costs with default value of zero
	state_costs <- array(0, dim = c(n_samples, n_treatments, n_states), 
	                     dimnames = list(NULL, treatment_names, state_names))

	# State costs are a mixture of treatment costs, acute event costs, event managment costs, 
	# and costs of transient events
	
	for(i_treatment in 1:n_treatments) {
		# First add the treatment costs (applies to all states except death)
		state_costs[, i_treatment, 1:(n_states - 1)] <-
		  input_parameters[, grep(paste(treatment_names[i_treatment], "cost"), 
		                          colnames(input_parameters))]

		# Then add acute event costs (same for all states except death)
		# Acute event costs are averaged over all events
		hazards_temp <- generate_hazards(input_parameters, i_treatment)
		# Use competing risks formula to convert to probabilities
		probs_temp <- (1 - exp(-0.25 * rowSums(hazards_temp))) * hazards_temp / rowSums(hazards_temp)
		# Weighted average of acute costs using probabilities of event as weight
		# Cost of AF well and death are zero
		acute_costs <- rowSums(input_parameters[, paste(event_names[-5], "acute cost")] * probs_temp[, -c(1,6)])
		state_costs[, i_treatment, 1:(n_states - 1)] <- 
		  state_costs[, i_treatment, 1:(n_states - 1)] + acute_costs
		
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

	} # End loop over treatments

	return(state_costs)
} # End function
