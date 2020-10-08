################################################################################
## Project: Introduction to PLS
## Script name: .Rprofile
## Script purpose: This script loads the libaries required, checks the R version and welcomes the client to the model
## Date: August 2020
## Author: Joe Moss (joe.moss@york.ac.uk)
## Organisation: York Health Economics Consortium (YHEC)
################################################################################


## welcome message ====================================================================================================================================


  message("\n\nWelcome to the \"Intro to Patient Level Simulation Model\" built by York Health Economics Consortium.\n\n")
  message("This is a cost-effectiveness model design to assess \"Drug X\" for the treatment of Chronic kidney disease.\n\n")
  message("This is version 1.0 of the model.\n\n")
  message("The model was last updated on the 15th of August 2020.\n\n")
  message("The model was built using R 4.0.0 (2020-04-24).\n\n")
  message(paste("You are using ", R.Version()$version.string, ".\n\n", sep = ""))

  message("If you are using a version of R that that is not 4.0.0, We cannot guarantee the full functionality of the model.\n\n")
  




## Clean R ============================================================================================================================================

# Before running the model, it is a good idea to clear all previous inputs in the model
# This will prevent the model running with outdated inputs if they have been changed between runs of the model

rm(list = ls())


## Load functions =====================================================================================================================================

source("R Scripts/Functions.R")
