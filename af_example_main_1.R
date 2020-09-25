# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Script conduct cost-effectiveness analysis of atrial fibrillation using example Markov model

library("remotes")
library("fs")

my_lib <- path_home_r("R/win-library/my-cran-versions/")

# detach("package:BCEA", unload = TRUE)
# Load BCEA library to help analyse and visualise results
library("BCEA", lib.loc = my_lib)


# Load the necessary modules
source("generate_input_parameters_1.R")
source("generate_hazards_1.R")
source("generate_transition_matrix_2.R")
source("generate_state_utilities_2.R")
source("generate_state_costs_2.R")
source("generate_model_outputs_1.r")

# Define global simulation parameters
n_samples <- 1000

# Define global model structure parameters
n_states <- 6
state_names <- c("AF Well", "Stroke", "ICH", "MI", "Bleed", "Dead")

n_treatments <- 3
treatment_names <- c("Coumarin", "Apixaban", "Dabigatran")

event_names <- c("Stroke", "MI", "Bleed", "ICH", "Death", "SE", "TIA")

# Define global scenario parameters
initial_age <- 70
final_age <- 100

# Generate the input parameters
# This will be converted into transition matrix, state costs, and state utilities
input_parameters <- generate_input_parameters(n_samples = n_samples)

# Run the Markov model to get the model outputs
model_outputs <- generate_model_outputs(input_parameters, 
                                        initial_age = initial_age, 
                                        final_age = final_age)

####################################################################################################
## Quick manual check of resutls ###################################################################
####################################################################################################

# Expected net benefit at 20,000
with(model_outputs, colMeans(20000 * total_qalys - total_costs))


##################################################################################################################
## Now use BCEA to analyse the outputs   #########################################################################
##################################################################################################################
# Create a bcea object for the doac model
doac_bcea <- bcea(e = model_outputs$total_qalys,
                  c = model_outputs$total_costs, ref = 1,
                  interventions = treatment_names)

# Summarise the results
summary(doac_bcea, wtp = 20000)

# Plot cost-effectiveness plane using base graphics

## Coumarin vs apixaban 
#jpeg(file = paste0("ce_plane_coumarin_vs_apixaban_", n_samples, ".jpg"))

ceplane.plot(doac_bcea, comparison = 1, wtp = 20000, graph = "base")

#dev.off()

## Coumarin vs dabigatran
#jpeg(file = paste0("ce_plane_coumarin_vs_dabigatran_", n_samples, ".jpg"))

ceplane.plot(doac_bcea, comparison = 2, wtp = 20000, graph = "base")

#dev.off()

# For multiple treatment comparison
doac_multi_ce <- multi.ce(doac_bcea)

## Cost-effectiveness acceptability curve
#jpeg(file = paste0("ceac_", n_samples, ".jpg"))

mce.plot(doac_multi_ce, pos = c(1, 0.5), graph = c("base", "ggplot2"))

#dev.off()

## Cost-effectiveness frontier curve
#jpeg(file = paste0("ceaf_", n_samples, ".jpg"))

ceaf.plot(doac_multi_ce, graph = c("base", "ggplot2"))
legend("bottomright", lty = c(1), legend = c("One"))

#dev.off()


# contour(doac_bcea, comparison = 1, scale = 0.5, nlevels = 4, levels = NULL,
#         pos = c(1, 0), xlim = NULL, ylim = NULL, graph=c("base", "ggplot2"))
# 
# contour2(doac_bcea,       # uses the results of the economic evalaution 
#          wtp = 20000,  # selects the willingness-to-pay threshold
#          xl = NULL,    # assumes default values
#          yl = NULL     # assumes default values
# )


ib.plot(doac_bcea,        # uses the results of the economic evalaution 
        comparison = 1, # if more than 2 interventions, selects the 
        #  pairwise comparison 
        wtp = 20000,    # selects the relevant willingness 
        #  to pay (default: 25,000)
        graph = "base"  # uses base graphics
)











