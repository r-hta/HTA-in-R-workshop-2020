# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert input parameters to age-dependent transition matrices
# Output is a 3-dimensional array of n_samples x n_treatments x n_states x n_states


generate_transition_matrix <- function(input_parameters) {
  
  # Construct named array of transition probabilities
  transition_matrix <- array(dim = c(n_samples, n_treatments, n_states, n_states), 
                             dimnames = list(NULL, treatment_names, state_names, state_names))
  
  # Need a matrix specifying the severity hierarchy of events
  # i.e. which transitions are permitted
  less_severe_states <- matrix(0, nrow = n_states, ncol = n_states,
                               dimnames = list(state_names, state_names))
  # Stroke can have further stroke or an ICH, or die
  less_severe_states["Stroke", ] <- c(1, 0, 0, 1, 1, 0)
  # ICH can have further ICH or die (stroke is about the same cost but higher utility)
  less_severe_states["ICH", ] <- c(1, 1, 0, 1, 1, 0)
  # MI can have any event except bleed
  less_severe_states["MI", ] <- c(1, 0, 0, 0, 1, 0)
  
  
	for(i_treatment in 1:n_treatments) {
	  # Generate treatment specific hazards of each of the events
	  hazards_temp <- generate_hazards(input_parameters, i_treatment)
	  # Don't need SE or TIA as they do not affect state occupancy
	  hazards_temp <- hazards_temp[, 1:6]
	  # Reorder to match state_names
	  hazards_temp <- hazards_temp[, paste("Hazard", state_names)]
	  
	  # Now loop through states converting hazards to probabilities
	  for(i_state in 1:n_states) {
	    # Using a severity hierarchy to avoid an event improving utility or reducing costs
	    # AF well < bleed < MI < ICH/stroke < dead
	    # State specific hazards
	    hazards_temp_state <- hazards_temp
	    # Set hazard of having less severe state to zero
	    hazards_temp_state[, less_severe_states[i_state, ] == 1] <- 0
	    # Calculate competing risks probabiliites (cycle length is 0.25 years)
	    probs_temp <- (1 - exp(-0.25 * rowSums(hazards_temp_state))) * hazards_temp_state / rowSums(hazards_temp_state)
	    transition_matrix[, i_treatment, i_state, ] <- probs_temp
	    
	    # And ensure those without event stay in the same state
	    transition_matrix[, i_treatment, i_state, i_state] <- 
	      transition_matrix[, i_treatment, i_state, i_state] +
	      1 - rowSums(transition_matrix[, i_treatment, i_state, ])
	      	  } # End loop over states
	  # Set transitions away from dead to 0 and from dead to dead to 1 (absorbing state)
    transition_matrix[, i_treatment, "Dead", ] <- 0
    transition_matrix[, i_treatment, "Dead", "Dead"] <- 1
	} # End loop over treatments
	return(transition_matrix)
}


