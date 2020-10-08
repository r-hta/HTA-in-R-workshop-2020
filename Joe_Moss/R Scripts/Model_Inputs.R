################################################################################
## Project: Introduction to PLS
## Script name: Model_Inputs.R
## Script purpose: This script stores all of the model inputs
## Date: August 2020
## Author: Joe Moss (joe.moss@york.ac.uk)
## Organisation: York Health Economics Consortium (YHEC)
################################################################################

## Message ======================================================================================================================================================

# If this model was being used for a HTA submission or for a client  I would recommend having the model inputs supplied
# in some kind of interface, either Excel file or an Rshiny frontend.

# For this demonstration, all model inputs are supplied to the model in this script.


## Basic background =============================================================================================================================================

# The model simulates 5 treatments (1 intervention "Drug X" and 4 comparators "Drug A", "Drug B", "Drug C" and "Drug D)
# Inputs that contain a value for each treatment are coded in the order above


## Model Setup ==================================================================================================================================================

Input_MS_Patient_Number <- 5000 # Set the number of patients to be simulated in the model
Input_MS_WPay <- 20000 # The willingness to pay threshold
Input_MS_Time_Horizon <- 60 # The number of years to simulate the model for
Input_MS_Treat_Length <- 5 # The number of years the average individual can tollerate treatment
Input_MS_Dialysis_Threshold <- 30 # GFR threshold at which dialysis will begin


## Discount rates ===============================================================================================================================================

Input_DC_Cost <- 0.035
Input_DC_Benefits <- 0.035


## Baseline Characteristics =====================================================================================================================================

Input_BL_Age <- c(55.4, 8.4) # First input is the mean age, second input is the standard deviation
Input_BL_Prop_Female <- 0.578 # % of the population that are female
Input_BL_GFR <- c(48, 9.2) # Mean and standard deviation for baseline glomerular filtration rate


## Treatment effect =============================================================================================================================================

# Load a coda derived from a NMA to simulate variance in treatment effects

Input_TE_Coda <- as.list(readRDS("Inputs/Treatment_Effect_Coda.rds"))


## Time varying inputs ==========================================================================================================================================

Input_TV_GFR_Decline <- c(2.358, 0.144) # Mean and SD for decline in GFR decline if a patient is not on treatment


## Mortality ====================================================================================================================================================

Input_Mort_Rate <- c(1.04, 0.01) # Mean increase in relative risk of mortality applied for every 1 GFR lost, used to adjust general population mortality rates
Input_Mort_Gen <- as.list(readRDS("Inputs/ONS_Life_Table.rds")) # Read in ONS life table but convert to a list in order to improve performace later


## TRAE rates ===================================================================================================================================================

Input_TRAE_Stroke <-  list(
  Drug_X = 0.055,
  Drug_A = 0.034,
  Drug_B = 0.036,
  Drug_C = 0.044,
  Drug_D = 0.048
) # Rate of a stroke in all 5 treatments use in the model

Input_TRAE_MI  <- list(
  Drug_X = 0.023,
  Drug_A = 0.016,
  Drug_B = 0.019,
  Drug_C = 0.033,
  Drug_D = 0.022
) # Rate of MI in all 5 treatment arms


## Number of injections a year ==================================================================================================================================

# Injections follow a gamma distribution so we need a shape and scale parameter for each drug

Input_Inj_Shape <- c(80, 90, 120, 122, 111)
Input_Inj_Scale <- c(18, 12, 17, 16, 12)


## Resource use =================================================================================================================================================

Input_RU_Dialysis <- 29 # If a patient ever reaches this value as GFR declines they will begin dialysis


## Resource costs ===============================================================================================================================================

Input_RC_Treatment <- list(
  Drug_X = 230,
  Drug_A = 200,
  Drug_B = 260,
  Drug_C = 310,
  Drug_D = 180
) # Cost per injection of treatment
Input_RC_Stroke <- 4800 # Initial cost of stroke
Input_RC_Stroke_LT <- 4800 # Cost per year after event
Input_RC_MI <- 2500 # Initial cost of MI
Input_RC_MI_LT <- 700 # Cost per year after
Input_RC_Dialysis <- 25000 # cost per year


## HRQoL ========================================================================================================================================================

# A previous regression analysis has shown that HRQoL can be predicted by 2 variables log(GFR) and Sex
# Equation below
# HRQoL ~ 0.178646993 * log(GFR) + -0.038606670 * Sex(M) + 0.005240453

Input_HRQOL_Int <- 0.005240453
Input_HRQOL_LGFR <- 0.178646993
Input_HRQOL_SexM <- -0.038606670
Input_HRQOL_Chol <- readRDS("Inputs/HRQoL_Regression_Cholesky.rds") # Load the cholskey matrix saved previously to include variance in the HRQoL estimate


## Model Set up code ============================================================================================================================================

# These steps generate a unique number for the model inputs, as a result it allows us to save and load previous simulation results

# First clear old code if present

if(exists("Input_Model_Code")){
  rm("Input_Model_Code")
}

Input_Model_Code <- grep("Input_", ls(.GlobalEnv)) # Find all values in the environment that start with "Input_"
Input_Model_Code <- ls(.GlobalEnv)[Input_Model_Code][-c(6, 12, 26)] # restrict the list to all non-table based inputs
Input_Model_Code <- unlist(map(Input_Model_Code, ~get(.x))) # Extract all the values
Input_Model_Code <- sum(as.numeric(Input_Model_Code)) # Take sum of all inputs
Input_Model_Code <- str_replace_all(Input_Model_Code, "[^[:alnum:]]", "") # Remove decimal point to generate unique code
