# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert input parameters state utilities

generate_state_utilities <- function(input_parameters) {
  
	# Construct array of state utilities with default value 0
	state_utilities <- array(0, dim=c(n_samples, n_treatments, n_states), 
	                         dimnames = list(NULL, treatment_names, state_names))

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
	  
	  # Then add acute event disutilities (same for all states except death)
	  # Acute event disutilities are averaged over all events
	  hazards_temp <- generate_hazards(input_parameters, i_treatment)
	  # Use competing risks formula to convert to probabilities
	  probs_temp <- (1 - exp(-0.25 * rowSums(hazards_temp))) * hazards_temp / rowSums(hazards_temp)
	  # Weighted average of acute disutilities using probabilities of event as weight
	  # Disutility for AF well and death are zero
	  acute_disutilities <- rowSums(input_parameters[, paste(event_names[-5], "disutility")] * probs_temp[, -c(1,6)])
	  state_utilities[, i_treatment, 1:(n_states - 1)] <- state_utilities[, i_treatment, 1:(n_states - 1)] +  acute_disutilities
	} # End loop over treatments

	return(state_utilities)
} # End function