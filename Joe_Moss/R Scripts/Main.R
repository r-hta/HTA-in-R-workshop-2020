################################################################################
## Project: Introduction to PLS
## Script name: Main.R
## Script purpose: This script is used to run the simulation
## Date: August 2020
## Author: Joe Moss (joe.moss@york.ac.uk)
## Organisation: York Health Economics Consortium (YHEC)
################################################################################


##  Load model inputs =======================================================================================================================

source("R Scripts/Libraries.R") # You may need to remove the comments in the Libraries.R script in order to install the required packages
source("R Scripts/Model_Inputs.R")


## Patient Characteristics ==================================================================================================================

# Generate the characteristics that are unique to each patient for the required number of patients

Out_Patient_Char <- func_Patient_Char_List(Input_MS_Patient_Number, age = Input_BL_Age, sex = Input_BL_Prop_Female, gfr = Input_BL_GFR,
                       coda = Input_TE_Coda, gfr.decline = Input_TV_GFR_Decline, mortality.rate = Input_Mort_Rate,
                       injection.info = c(Input_Inj_Shape, Input_Inj_Scale), hrqol.info = list(Input_HRQOL_Int, Input_HRQOL_LGFR,
                                                                                               Input_HRQOL_SexM, Input_HRQOL_Chol))


## Run simulation ===========================================================================================================================

# Generate or load simulation results

Out_Sim_Results <- func_Simulation()

#View(as.data.frame(Out_Sim_Results[["Drug_X"]][[1]]))

# Condense summary results

Out_Sim_Summary <- func_Condense_Results()

# Produce a results table (use "save.table = TRUE" to save a copy of the table to the outputs directory)

func_CE_Result()

# Plot a CEP for Drug_X versus Drug_C

func_Plot_CEP(treat.name = "Drug_X", treat.name2 = "Drug_C",
              y.limits = c(-500000, 500000), y.breaks = seq(-500000, 500000, 50000),
              x.limits = c(-20, 20), x.breaks = seq(-20, 20, 2),
              txt.size = 16)

# Plot a CEAC

func_Plot_CEAC(treat.names = c("Drug_X", "Drug_A", "Drug_B", "Drug_C", "Drug_D"),
               wtp = seq(0,100000,1000), y.limits = c(0,0.5), y.breaks = seq(0,0.5,0.05))


## Check model stability ===================================================================================================================

func_Model_Stability("Cost", y.limits = c(0, 200000), y.breaks = seq(0, 200000, 20000), interactive = TRUE)
func_Model_Stability("HRQoL", y.limits = c(5, 15), y.breaks = seq(5, 15, 1), interactive = TRUE)
