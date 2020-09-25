# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert model inputs to ouputs using Markov model

# TODO: See if apply statement (line 50) for total costs works

generate_model_outputs <- function(input_parameters, initial_age = 70, final_age = 100) {
  # Calculate number of cycles
  n_cycles <- (final_age - initial_age) * 4
  # Sample transition probability matrix
	transition_matrix <- generate_transition_matrix(input_parameters)
	# Sample the utilities associated with each state
	state_utilities <- generate_state_utilities(input_parameters)
	# Sample the costs associated with each state
	state_costs <- generate_state_costs(input_parameters)
	
	# Create an array of cohort vectors. 
	cohort_array <- array(0, dim = c(n_samples, n_treatments, n_states), 
	                      dimnames = list(NULL, treatment_names, state_names))
	# Initialise the cohort vector with all patients starting in AF Well
	cohort_array[, , "AF Well"] <- 1
	
	# Store the total costs and QALYs
	total_costs <- total_qalys <- matrix(0, nrow = n_samples, ncol = n_treatments)
	colnames(total_costs) <- colnames(total_qalys) <- treatment_names
	
	for(i_cycle in 1:n_cycles) {
	  # Calculate the discount factor (decreases by 1/1.035 every 4 cycles)
	  discount_factor <- (1 / 1.035) ^ floor((i_cycle / 4))
	  for(i_treatment in 1:n_treatments) {
	    # Calculate state occupancy costs and QALYs for this cycle and treatment
	    total_costs[, i_treatment] <- total_costs[, i_treatment] +
	        discount_factor * rowSums(cohort_array[, i_treatment, ] * state_costs[, i_treatment, ])
	    total_qalys[, i_treatment] <- total_qalys[, i_treatment] +
	        discount_factor * rowSums(cohort_array[, i_treatment, ] * state_utilities[, i_treatment, ])
	   
	    # Half cycle correction (subtract half of initial QALYs and costs)
	    if(i_cycle==1){
	      total_costs[, i_treatment] <- total_costs[, i_treatment] * 0.5
	      total_qalys[, i_treatment] <- total_qalys[, i_treatment] * 0.5
	    }
	    # Now update the cohort vector using the transition probability matrix.
	    # Preferable to do this using array arithmetic, rather than a slow for loop.
	    for(i_sample in 1:n_samples) {
	      cohort_array[i_sample, i_treatment, ] <-
	        cohort_array[i_sample, i_treatment, ] %*% transition_matrix[i_sample, i_treatment, , ]
	    }
	  } # End loop over treatments
	} # End loop over cycles
	
	# Half cycle correction (add half of final QALYs and costs)
	total_costs[, i_treatment] <- total_costs[, i_treatment] + 
	  0.5 * discount_factor * rowSums(cohort_array[, i_treatment, ] * state_costs[, i_treatment, ])
	total_qalys[, i_treatment] <- total_qalys[, i_treatment] + 
	  0.5 * discount_factor * rowSums(cohort_array[, i_treatment,] * state_utilities[, i_treatment, ])
	
	model_outputs <- list("total_costs" = total_costs, "total_qalys" = total_qalys)
	return(model_outputs)
} # End function