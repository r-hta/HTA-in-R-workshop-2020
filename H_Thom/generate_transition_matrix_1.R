# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert input parameters to age-dependent transition matrices
# Output is a 3-dimensional array of n_samples x n_treatments x n_states x n_states
# Age dependent as is adjusted for (increasing) mortality

# TODO check competing risks formula is correct
# TODO make probability of death decrease with age (Scale with life tables)

# TODO Force hierarchy of state transitions
# AF well < bleed < MI < ICH/stroke < dead


generate_transition_matrix <- function(input_parameters, age = age, sensitivity = "None") {
  
  # Need a matrix specifying less severe states
  less_severe_states <- matrix(0, nrow = n_states, ncol = n_states,
                               dimnames = list(state_names, state_names))
  less_severe_states["Stroke", ] <- c(1, 0, 0, 1, 1, 0)
  less_severe_states["ICH", ] <- c(1, 0, 0, 1, 1, 0)
  less_severe_states["MI", ] <- c(1, 0, 0, 0, 1, 0)
  
  # Construct named array of transition probabilities
	transition_matrix <- array(dim = c(n_samples, n_treatments, n_states, n_states), 
	                           dimnames = list(NULL, treatment_names, state_names, state_names))
	for(i_treatment in 1:n_treatments) {
	  # Create temporary hazard and probability objects before putting them in the transition_matrix
	  # Baseline hazards for Coumarin
	  hazards_temp <- exp(input_parameters[, grepl("Log hazard", colnames(input_parameters))])
	  # If treatment is not Coumarin need to add appropriate log hazard ratio
	  if(i_treatment != 1) {
	    hazards_temp <- hazards_temp * 
	      exp(input_parameters[, grepl(paste("Log HR", treatment_names[i_treatment]), 
	                                   colnames(input_parameters))])
	  }
	  # Give the hazards appropriate names
	  colnames(hazards_temp) <- paste("Hazard", event_names)
	  # Don't need SE or TIA as they do not affect state occupancy
	  hazards_temp <- hazards_temp[, 1:5]
	  # Calculate competing risks probabiliites (cycle length is 0.25 years)
	  probs_temp <- (1 - exp(-0.25 * rowSums(hazards_temp))) * hazards_temp / rowSums(hazards_temp)
	  # Add a probability of staying in AF well (initialised to zero)
	  probs_temp <- cbind(0, probs_temp)
	  # And name correctly (note event is 'death' but state is 'dead')
	  colnames(probs_temp) <- paste("Probability", c("AF Well", event_names[1:4], "Dead"))
	  # Reorder to match state_names
	  probs_temp[, paste("Probability", state_names)]
	  # Now store in correct location (all probabilities are the same except for death)
	  for(i_state in 1:n_states) {
	    # Need to put competing risks formula here
	    # TODO: And at this stage impose the severity hierarchy
	    # AF well < bleed < MI < ICH/stroke < dead
	    transition_matrix[, i_treatment, i_state, ] <- probs_temp[, paste("Probability",state_names)]
	    
	    # TODO Force hierarchy of state transitions
	    # This is a crude approach. Method above (incorporating competing risks) is preferable.
	    #transition_matrix[, i_treatment, i_state, less_severe_states[i_state, ] == 1] <- 0
	    
	    
	    # And ensure those without event stay in the same state
	    # TODO :This shouldn't be summed over probs_temp but over transition_matrix itself
	    transition_matrix[, i_treatment, i_state, i_state] <- 
	      transition_matrix[, i_treatment, i_state, i_state] +
	      1 - rowSums(probs_temp[, paste("Probability",state_names)])
	  } # End loop over states
	  # Set transitions away from dead to 0 and from dead to dead to 1 (absorbing state)
    transition_matrix[, i_treatment, "Dead", ] <- 0
    transition_matrix[, i_treatment, "Dead", "Dead"] <- 1
	} # End loop over treatments
	return(transition_matrix)
}


