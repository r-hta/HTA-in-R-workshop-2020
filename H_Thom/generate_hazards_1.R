# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to generate treatment specific hazards of acute events
# Uses the coumarin hazards and hazard ratios from input_parameters
# Output is a matrix with n_samples rows and 8 columns (one for each event)

generate_hazards <- function(input_parameters, i_treatment) {
  # Create temporary hazard and probability objects before putting them in the transition_matrix
  # Baseline hazards for Coumarin
  hazards_temp <- exp(input_parameters[, grepl("Log hazard", colnames(input_parameters))])
  # If treatment is not Coumarin need to add appropriate log hazard ratio
  if(i_treatment != 1) {
    hazards_temp <- hazards_temp * 
      exp(input_parameters[, grepl(paste("Log HR", treatment_names[i_treatment]), 
                                   colnames(input_parameters))])
  }
  # Append a zero for "AF well"
  hazards_temp <- cbind(0, hazards_temp)
  # And name correctly (note event is 'death' but state is 'dead')
  colnames(hazards_temp) <- paste("Hazard", c("AF Well", event_names[1:4], "Dead", event_names[6:7])) 
  return(hazards_temp)
}