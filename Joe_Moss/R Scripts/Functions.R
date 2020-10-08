################################################################################
## Project: Introduction to PLS
## Script name: Functions.R
## Script purpose: This script is where all of the model functions are found
## Date: August 2020
## Author: Joe Moss (joe.moss@york.ac.uk)
## Organisation: York Health Economics Consortium (YHEC)
################################################################################


## Baseline character function ====================================================================================================================================

# Generates random patient characteristics

func_Patient_Char_Gen <- function(age, # BL age (mean and SD)
                                  sex, # gender
                                  gfr, # BL gfr (mean and SD)
                                  coda, # Treatment effect coda
                                  gfr.decline, # rate of decline in gfr when not on treatment
                                  mortality.rate, # individual mortality rate
                                  injection.info, # shape and scale info for injection numbers
                                  hrqol.info){ # coefficients and cholesky
  
  ## Step 1 - Create a list to fill
  
  Temp_List <- list(
    BL_Age = 0,
    Sex = 0,
    GFR = 0,
    GFR_Decline = 0,
    Coda = 0,
    Treatment_Effect = list(
      Drug_X = 0,
      Drug_A = 0,
      Drug_B = 0,
      Drug_C =0,
      Drug_D =0
    ),
    Mortality_Rate = 0,
    Injection_N = list(
      Drug_X = 0,
      Drug_A = 0,
      Drug_B = 0,
      Drug_C =0,
      Drug_D =0
    ),
    HRQoL_Info = list(
      Intercept = 0,
      Log_GFR = 0,
      Sex = 0,
      Cholesky = 0
    )
  )
  
  ## Step 2 - generate results
  
  Temp_List[["BL_Age"]] <- rnorm(1, age[1], age[2])
  Temp_List[["Sex"]] <- sample(c(0, 1), 1, prob = c(sex, 1 - sex))
  Temp_List[["GFR"]] <- rnorm(1, gfr[1], gfr[2])
  Temp_List[["GFR_Decline"]] <- rnorm(1, gfr.decline[1], gfr.decline[2])
  
  Temp_List[["Coda"]] <- sample(c(1:length(coda$Drug_A)), 1)
  Temp_List[["Treatment_Effect"]][["Drug_X"]] <- coda$Drug_X[Temp_List[["Coda"]]]
  Temp_List[["Treatment_Effect"]][["Drug_A"]] <- coda$Drug_A[Temp_List[["Coda"]]]
  Temp_List[["Treatment_Effect"]][["Drug_B"]] <- coda$Drug_B[Temp_List[["Coda"]]]
  Temp_List[["Treatment_Effect"]][["Drug_C"]] <- coda$Drug_C[Temp_List[["Coda"]]]
  Temp_List[["Treatment_Effect"]][["Drug_D"]] <- coda$Drug_D[Temp_List[["Coda"]]]
  
  Temp_List[["Mortality_Rate"]] <- rnorm(1, mortality.rate[1], mortality.rate[2])
  
  Temp_List[["Injection_N"]][["Drug_X"]] <- rgamma(n = 1, shape = injection.info[1], rate = injection.info[6])
  Temp_List[["Injection_N"]][["Drug_A"]] <- rgamma(n = 1, shape = injection.info[2], rate = injection.info[7])
  Temp_List[["Injection_N"]][["Drug_B"]] <- rgamma(n = 1, shape = injection.info[3], rate = injection.info[8])
  Temp_List[["Injection_N"]][["Drug_C"]] <- rgamma(n = 1, shape = injection.info[4], rate = injection.info[9])
  Temp_List[["Injection_N"]][["Drug_D"]] <- rgamma(n = 1, shape = injection.info[5], rate = injection.info[10])
  
  
  Temp_List[["HRQoL_Info"]][["Cholesky"]] <- as.vector(t(hrqol.info[[4]]))
  
  
  ## Step 3 - Generate some random coefficient for the HRQoL regression
  
  z <- rnorm(n = 3, mean = 0, sd = 1)
  
  Temp_List[["HRQoL_Info"]][["Intercept"]] <- hrqol.info[[1]] + sum(Temp_List[["HRQoL_Info"]][["Cholesky"]][1:3]*z)
  Temp_List[["HRQoL_Info"]][["Log_GFR"]] <- hrqol.info[[2]] + sum(Temp_List[["HRQoL_Info"]][["Cholesky"]][4:6]*z)
  Temp_List[["HRQoL_Info"]][["Sex"]] <- hrqol.info[[3]] + sum(Temp_List[["HRQoL_Info"]][["Cholesky"]][7:9]*z)
  
  ## Step 4 - Return list
  
  return(Temp_List)
  
}


## List ============================================================================================================================================================

# Generates random patient characteristics for the required number of patients to be simulated

func_Patient_Char_List <- function(patient.n, use.seed = FALSE, seed.value = 1, ...){
  
  if(use.seed){
    set.seed(seed.value)
  }
  
  Temp_list <- map(.x = 1:patient.n, .f = function(x) func_Patient_Char_Gen(...)
                   )
  
  return(Temp_list)
  
}


## GFR over time ===================================================================================================================================================

# Simple loop to track a patient's GFR depending on treatment status

func_Treat_Effect <- function(bl.gfr, gfr.decline, treat.effect, treat.name, cycle.num, treat.check,...){
  
  # Calc Baseline GFR and treatment effect once instead of every time its needed
  
  GFR_Increase <- bl.gfr + treat.effect[[treat.name]]
  
  # Loop and logic
  
  V_GFR_Change = vector(mode = "numeric", length = max(cycle.num)+1)
  
  for(i in 1:(max(cycle.num)+1)){
    
    if(cycle.num[i] == 0){V_GFR_Change[i] <- bl.gfr}
    
    else if(cycle.num[i] > 0 & treat.check[i] == 1){V_GFR_Change[i] <- GFR_Increase}
    
    else if(cycle.num[i] > 0 & treat.check[i] == 0){V_GFR_Change[i] <- max(V_GFR_Change[i-1] - gfr.decline, 0)}
    
  }
  
  # return output
  
  return(V_GFR_Change)
  
}


## Moratlity ========================================================================================================================================================

# Random number check against mortality rate

func_Mortality <- function(age, sex, mort.gen, gfr, mort.rr){
  
  # Find general morality rate for the age and sex
  
  Gen_Mort_Rate <- mort.gen[[3-sex]][floor(map_dbl(age, ~min(.x, 105)))+1]
  
  # Adjust mortality rate based on per 5 GFR change
  
  Adj_Mort_RR <- map_dbl(gfr, ~mort.rr ^ ((100-.x)/5))
  
  # Adjust general morality rates using the rr
  
  Adj_Gen_Mort_Rate <- Gen_Mort_Rate * Adj_Mort_RR
  
  # Check whether individual lives of dies
  
  Mort_Out <- ifelse(runif(length(age)) < Adj_Gen_Mort_Rate, 1, 0)
  
  return(Mort_Out)
  
}


## Discounting ======================================================================================================================================================

# Formula to discount costs and benefits

func_Discount <- function(item, discon.rate, cycle.n){
  
  #Temport store of item
  
  Temp_Cost <- item
  
  # Apply discount
  
  Temp_Cost <- Temp_Cost / (1+discon.rate)^cycle.n
  
  # return item
  
  return(Temp_Cost)
  
}


## HRQoL ===========================================================================================================================================================

# Formula to calculate utility values based on gfr

func_HRQoL <- function(hrqol.info, gfr, sex){
  
  # Check for 0 gfr as it cannot be logged
  
  if(gfr == 0){
    Temp_HRQoL <- -0.5
  }
  
  # HRQoL ~ 0.178646993 * log(GFR) + -0.038606670 * Sex(M) + 0.005240453
  
  else{
    Temp_HRQoL <- as.numeric(hrqol.info[2]) * log(gfr) + as.numeric(hrqol.info[3]) * sex + as.numeric(hrqol.info[1])
  }
  
  return(Temp_HRQoL)
  
}


## Total costs =====================================================================================================================================================

# SImple loop to find the total costs in each cycle

func_Total_Cost <- function(cost.treat, cost.stroke, cost.mi, cost.dialysis){
  
  V_Total_Cost <- vector("numeric", length = Input_MS_Time_Horizon+1)
  
  for(i in 1:(Input_MS_Time_Horizon+1)){
    
    V_Total_Cost[i] <- sum(cost.treat[i], cost.stroke[i], cost.mi[i], cost.dialysis[i])
  }
  
  return(V_Total_Cost)
  
}


## Main logic function =============================================================================================================================================

# Runs each patient through the simulation by Treatment name

func_Paient_Results_Single <- function(treat.name, time.horizon, bl.age, treat.len, inj.num, bl.gfr, gfr.decline, treat.effect,
                                       mort.rate, sex, mort.gen, trae.stroke, trae.mi, dialysis.start, cost.trt, cost.stroke,
                                       cost.mi, cost.dialysis, discon.rate, hrqol.info){
  
  # Create a series of empty vectors to fill
  
  Temp_List <- list(
    V_Cycle_N = vector(mode = "numeric", length = time.horizon+1),
    V_Age = vector(mode = "numeric", length = time.horizon+1),
    V_Treat_Check = vector(mode = "numeric", length = time.horizon+1),
    V_Inj_N = vector(mode = "numeric", length = time.horizon+1),
    V_GFR = vector(mode = "numeric", length = time.horizon+1),
    
    V_Mortality = vector(mode = "numeric", length = time.horizon+1),
    V_Mortality_Lag = vector(mode = "numeric", length = time.horizon+1),
    
    V_AE_Stroke = vector(mode = "numeric", length = time.horizon+1),
    V_AE_Stroke_Lag = vector(mode = "numeric", length = time.horizon+1),
    V_AE_MI = vector(mode = "numeric", length = time.horizon+1),
    V_AE_MI_Lag = vector(mode = "numeric", length = time.horizon+1),
    
    V_Dialysis = vector(mode = "numeric", length = time.horizon+1),
    
    V_Cost_Treat = vector(mode = "numeric", length = time.horizon+1),
    V_Cost_Stroke = vector(mode = "numeric", length = time.horizon+1),
    V_Cost_MI = vector(mode = "numeric", length = time.horizon+1),
    V_Cost_Dialysis = vector(mode = "numeric", length = time.horizon+1),
    
    V_Cost_Disc_Treat = vector(mode = "numeric", length = time.horizon+1),
    V_Cost_Disc_Stroke = vector(mode = "numeric", length = time.horizon+1),
    V_Cost_Disc_MI = vector(mode = "numeric", length = time.horizon+1),
    V_Cost_Disc_Dialysis = vector(mode = "numeric", length = time.horizon+1),
    
    V_Cost_Disc_Total = vector(mode = "numeric", length = time.horizon+1),
    
    V_HRQoL = vector(mode = "numeric", length = time.horizon+1),
    V_HRQoL_Disc = vector(mode = "numeric", length = time.horizon+1)
    
  )
  
  # Fill vectors based on inputs and results of other functions
  
  Temp_List[["V_Cycle_N"]] <- 0:time.horizon
  
  Temp_List[["V_Age"]] <- bl.age + Temp_List[["V_Cycle_N"]]
  
  Temp_List[["V_Treat_Check"]] <- ifelse(Temp_List[["V_Cycle_N"]] <= treat.len, 
                                         ifelse(Temp_List[["V_Cycle_N"]] == 0, 0, 1),
                                         0)
  
  Temp_List[["V_Inj_N"]] <- ifelse(Temp_List[["V_Treat_Check"]] == 1, inj.num[[treat.name]], 0)
  
  Temp_List[["V_GFR"]] <- func_Treat_Effect(cycle.num = 0:time.horizon, treat.check = Temp_List[["V_Treat_Check"]], bl.gfr, gfr.decline, treat.effect, treat.name)
  
  Temp_List[["V_Mortality"]] <- func_Mortality(age = Temp_List[["V_Age"]], sex, mort.gen, gfr = Temp_List[["V_GFR"]], mort.rate)
  
  Temp_List[["V_Mortality_Lag"]] <- ifelse(cumsum(cumsum(Temp_List[["V_Mortality"]])) > 1, NA, 1)
  
  Temp_List[["V_AE_Stroke"]] <- map_dbl(Temp_List[["V_Treat_Check"]], ~ifelse(.x == 1,
                                                                              ifelse(runif(1) < trae.stroke[[treat.name]], 1, 0),
                                                                              0))
  
  Temp_List[["V_AE_Stroke_Lag"]] <- ifelse(cumsum(cumsum(Temp_List[["V_AE_Stroke"]])) == 0, 0,
                                           ifelse(cumsum(cumsum(Temp_List[["V_AE_Stroke"]])) == 1, 1,
                                                  ifelse(cumsum(cumsum(Temp_List[["V_AE_Stroke"]])) > 1, 2, NA
                                                  ))
  )
  
  Temp_List[["V_AE_MI"]] <- map_dbl(Temp_List[["V_Treat_Check"]], ~ifelse(.x == 1,
                                                                          ifelse(runif(1) < trae.mi[[treat.name]], 1, 0),
                                                                          0))
  
  Temp_List[["V_AE_MI_Lag"]] <- ifelse(cumsum(cumsum(Temp_List[["V_AE_MI"]])) == 1, 1,
                                       ifelse(cumsum(cumsum(Temp_List[["V_AE_MI"]])) > 1, 2, 0
                                       )
  )
  
  Temp_List[["V_Dialysis"]] <- ifelse(Temp_List[["V_GFR"]] <= dialysis.start, 1, 0)
  
  
  # Fill cost vectors
  
  Temp_List[["V_Cost_Treat"]] <- Temp_List[["V_Inj_N"]] * cost.trt[[treat.name]]
  
  Temp_List[["V_Cost_Stroke"]] <- ifelse(Temp_List[["V_AE_Stroke_Lag"]] == 1, cost.stroke[1],
                                         ifelse(Temp_List[["V_AE_Stroke_Lag"]] == 2, cost.stroke[2],0))
  
  Temp_List[["V_Cost_MI"]] <- ifelse(Temp_List[["V_AE_MI_Lag"]] == 1, cost.mi[1],
                                     ifelse(Temp_List[["V_AE_MI_Lag"]] == 2, cost.mi[2],0))
  
  Temp_List[["V_Cost_Dialysis"]] <- ifelse(Temp_List[["V_Dialysis"]] == 1, cost.dialysis, 0)
  
  # Discounting costs
  
  Temp_List[["V_Cost_Disc_Treat"]] <- cumsum(map2_dbl(.x = Temp_List[["V_Cost_Treat"]],.y = Temp_List[["V_Cycle_N"]], ~func_Discount(.x, discon.rate[1], .y)))
  
  Temp_List[["V_Cost_Disc_Stroke"]] <- cumsum(map2_dbl(.x = Temp_List[["V_Cost_Stroke"]],.y = Temp_List[["V_Cycle_N"]], ~func_Discount(.x, discon.rate[1], .y)))
  
  Temp_List[["V_Cost_Disc_MI"]] <- cumsum(map2_dbl(.x = Temp_List[["V_Cost_MI"]],.y = Temp_List[["V_Cycle_N"]], ~func_Discount(.x, discon.rate[1], .y)))
  
  Temp_List[["V_Cost_Disc_Dialysis"]] <- cumsum(map2_dbl(.x = Temp_List[["V_Cost_Dialysis"]],.y = Temp_List[["V_Cycle_N"]], ~func_Discount(.x, discon.rate[1], .y)))
  
  Temp_List[["V_Cost_Disc_Total"]] <- func_Total_Cost(cost.treat = Temp_List[["V_Cost_Disc_Treat"]], cost.stroke = Temp_List[["V_Cost_Disc_Stroke"]],
                                                      cost.mi = Temp_List[["V_Cost_Disc_MI"]], cost.dialysis = Temp_List[["V_Cost_Disc_Dialysis"]])
  
  # HRQoL
  
  Temp_List[["V_HRQoL"]] <- map_dbl(.x = Temp_List[["V_GFR"]], ~func_HRQoL(hrqol.info, .x, sex))
  
  Temp_List[["V_HRQoL_Disc"]] <- cumsum(map2_dbl(.x = Temp_List[["V_HRQoL"]],.y = Temp_List[["V_Cycle_N"]], ~func_Discount(.x, discon.rate[2], .y)))
  
  
  # Remove any results post death
  
  Temp_List <- map(Temp_List, ~Temp_List[["V_Mortality_Lag"]]*.x)
  
  
  
  # return output
  
  return(Temp_List)
  
}


## Main simulation function =============================================================================================================================================

# This function runs the entire simulation
# It generates results for all patients in all 5 treatment arms
# It also has an option to laod previously generated results

func_Simulation <- function(use.seed = FALSE,
                            seed.value = 1,
                            ...) {
  # Set seed if required
  
  if (use.seed) {
    set.seed(seed.value)
  }
  
  # remove old sim results to prevent large list being passed to parallel loops
  if (exists("Out_Sim_Results", envir = .GlobalEnv)) {
    rm("Out_Sim_Results", envir = .GlobalEnv)
  }
  
  # Check if previous sim has been saved
  
  # Unique input string
  
  File_String <- Input_Model_Code
  
  Dir_Test <-
    file.exists(paste("Outputs/Sim_", File_String, "/Output.rds", sep = ""))
  
  # If previous results do not exist run the simulation
  
  if (Dir_Test == FALSE) {
    # Ask user is they want to save results for next time
    
    Save_Result <-
      askYesNo("Would you like to save the outputs of the simulation for a quick load next time?")
    
    # Set up parallel
    
    cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
    registerDoSNOW(cl)
    
    # Fill list ===================================================
    
    message("Starting simulation")
    pb = txtProgressBar(
      min = 0,
      max = 5,
      initial = 0,
      style = 3
    )
    
    message("\n\nGenerating results for Drug_X")
    # Drug X
    
    Temp_List_Drug_X <- foreach(
      i = 1:Input_MS_Patient_Number,
      .packages = c("purrr"),
      .export = c(ls(.GlobalEnv)),
      .inorder = FALSE
    ) %dopar% {
      func_Paient_Results_Single(
        treat.name = "Drug_X",
        time.horizon = Input_MS_Time_Horizon,
        bl.age =  Out_Patient_Char[[i]][["BL_Age"]],
        treat.len = Input_MS_Treat_Length,
        inj.num = Out_Patient_Char[[i]][["Injection_N"]],
        bl.gfr = Out_Patient_Char[[i]][["GFR"]],
        gfr.decline = Out_Patient_Char[[i]][["GFR_Decline"]],
        treat.effect = Out_Patient_Char[[i]][["Treatment_Effect"]],
        mort.rate = Out_Patient_Char[[i]][["Mortality_Rate"]],
        sex = Out_Patient_Char[[i]][["Sex"]],
        mort.gen = Input_Mort_Gen,
        trae.stroke = Input_TRAE_Stroke,
        trae.mi = Input_TRAE_MI,
        dialysis.start = Input_MS_Dialysis_Threshold,
        cost.trt = Input_RC_Treatment,
        cost.stroke = c(Input_RC_Stroke, Input_RC_Stroke_LT),
        cost.mi = c(Input_RC_MI, Input_RC_MI_LT),
        cost.dialysis = Input_RC_Dialysis,
        discon.rate = c(Input_DC_Cost, Input_DC_Benefits),
        hrqol.info = Out_Patient_Char[[i]][["HRQoL_Info"]][-4]
      )
    }
    
    setTxtProgressBar(pb, 1)
    
    # Drug A
    
    message("\n\nGenerating results for Drug_A")
    
    Temp_List_Drug_A <- foreach(
      i = 1:Input_MS_Patient_Number,
      .packages = c("purrr"),
      .export = c(ls(.GlobalEnv)),
      .inorder = FALSE
    ) %dopar% {
      func_Paient_Results_Single(
        treat.name = "Drug_A",
        time.horizon = Input_MS_Time_Horizon,
        bl.age =  Out_Patient_Char[[i]][["BL_Age"]],
        treat.len = Input_MS_Treat_Length,
        inj.num = Out_Patient_Char[[i]][["Injection_N"]],
        bl.gfr = Out_Patient_Char[[i]][["GFR"]],
        gfr.decline = Out_Patient_Char[[i]][["GFR_Decline"]],
        treat.effect = Out_Patient_Char[[i]][["Treatment_Effect"]],
        mort.rate = Out_Patient_Char[[i]][["Mortality_Rate"]],
        sex = Out_Patient_Char[[i]][["Sex"]],
        mort.gen = Input_Mort_Gen,
        trae.stroke = Input_TRAE_Stroke,
        trae.mi = Input_TRAE_MI,
        dialysis.start = Input_MS_Dialysis_Threshold,
        cost.trt = Input_RC_Treatment,
        cost.stroke = c(Input_RC_Stroke, Input_RC_Stroke_LT),
        cost.mi = c(Input_RC_MI, Input_RC_MI_LT),
        cost.dialysis = Input_RC_Dialysis,
        discon.rate = c(Input_DC_Cost, Input_DC_Benefits),
        hrqol.info = Out_Patient_Char[[i]][["HRQoL_Info"]][-4]
      )
    }
    
    setTxtProgressBar(pb, 2)
    
    # Drug B
    
    message("\n\nGenerating results for Drug_B")
    
    Temp_List_Drug_B <- foreach(
      i = 1:Input_MS_Patient_Number,
      .packages = c("purrr"),
      .export = c(ls(.GlobalEnv)),
      .inorder = FALSE
    ) %dopar% {
      func_Paient_Results_Single(
        treat.name = "Drug_B",
        time.horizon = Input_MS_Time_Horizon,
        bl.age =  Out_Patient_Char[[i]][["BL_Age"]],
        treat.len = Input_MS_Treat_Length,
        inj.num = Out_Patient_Char[[i]][["Injection_N"]],
        bl.gfr = Out_Patient_Char[[i]][["GFR"]],
        gfr.decline = Out_Patient_Char[[i]][["GFR_Decline"]],
        treat.effect = Out_Patient_Char[[i]][["Treatment_Effect"]],
        mort.rate = Out_Patient_Char[[i]][["Mortality_Rate"]],
        sex = Out_Patient_Char[[i]][["Sex"]],
        mort.gen = Input_Mort_Gen,
        trae.stroke = Input_TRAE_Stroke,
        trae.mi = Input_TRAE_MI,
        dialysis.start = Input_MS_Dialysis_Threshold,
        cost.trt = Input_RC_Treatment,
        cost.stroke = c(Input_RC_Stroke, Input_RC_Stroke_LT),
        cost.mi = c(Input_RC_MI, Input_RC_MI_LT),
        cost.dialysis = Input_RC_Dialysis,
        discon.rate = c(Input_DC_Cost, Input_DC_Benefits),
        hrqol.info = Out_Patient_Char[[i]][["HRQoL_Info"]][-4]
      )
    }
    
    setTxtProgressBar(pb, 3)
    
    # Drug C
    
    message("\n\nGenerating results for Drug_C")
    
    Temp_List_Drug_C <- foreach(
      i = 1:Input_MS_Patient_Number,
      .packages = c("purrr"),
      .export = c(ls(.GlobalEnv)),
      .inorder = FALSE
    ) %dopar% {
      func_Paient_Results_Single(
        treat.name = "Drug_C",
        time.horizon = Input_MS_Time_Horizon,
        bl.age =  Out_Patient_Char[[i]][["BL_Age"]],
        treat.len = Input_MS_Treat_Length,
        inj.num = Out_Patient_Char[[i]][["Injection_N"]],
        bl.gfr = Out_Patient_Char[[i]][["GFR"]],
        gfr.decline = Out_Patient_Char[[i]][["GFR_Decline"]],
        treat.effect = Out_Patient_Char[[i]][["Treatment_Effect"]],
        mort.rate = Out_Patient_Char[[i]][["Mortality_Rate"]],
        sex = Out_Patient_Char[[i]][["Sex"]],
        mort.gen = Input_Mort_Gen,
        trae.stroke = Input_TRAE_Stroke,
        trae.mi = Input_TRAE_MI,
        dialysis.start = Input_MS_Dialysis_Threshold,
        cost.trt = Input_RC_Treatment,
        cost.stroke = c(Input_RC_Stroke, Input_RC_Stroke_LT),
        cost.mi = c(Input_RC_MI, Input_RC_MI_LT),
        cost.dialysis = Input_RC_Dialysis,
        discon.rate = c(Input_DC_Cost, Input_DC_Benefits),
        hrqol.info = Out_Patient_Char[[i]][["HRQoL_Info"]][-4]
      )
    }
    
    setTxtProgressBar(pb, 4)
    
    # Drug D
    
    message("\n\nGenerating results for Drug_D")
    
    Temp_List_Drug_D <- foreach(
      i = 1:Input_MS_Patient_Number,
      .packages = c("purrr"),
      .export = c(ls(.GlobalEnv)),
      .inorder = FALSE
    ) %dopar% {
      func_Paient_Results_Single(
        treat.name = "Drug_D",
        time.horizon = Input_MS_Time_Horizon,
        bl.age =  Out_Patient_Char[[i]][["BL_Age"]],
        treat.len = Input_MS_Treat_Length,
        inj.num = Out_Patient_Char[[i]][["Injection_N"]],
        bl.gfr = Out_Patient_Char[[i]][["GFR"]],
        gfr.decline = Out_Patient_Char[[i]][["GFR_Decline"]],
        treat.effect = Out_Patient_Char[[i]][["Treatment_Effect"]],
        mort.rate = Out_Patient_Char[[i]][["Mortality_Rate"]],
        sex = Out_Patient_Char[[i]][["Sex"]],
        mort.gen = Input_Mort_Gen,
        trae.stroke = Input_TRAE_Stroke,
        trae.mi = Input_TRAE_MI,
        dialysis.start = Input_MS_Dialysis_Threshold,
        cost.trt = Input_RC_Treatment,
        cost.stroke = c(Input_RC_Stroke, Input_RC_Stroke_LT),
        cost.mi = c(Input_RC_MI, Input_RC_MI_LT),
        cost.dialysis = Input_RC_Dialysis,
        discon.rate = c(Input_DC_Cost, Input_DC_Benefits),
        hrqol.info = Out_Patient_Char[[i]][["HRQoL_Info"]][-4]
      )
    }
    
    setTxtProgressBar(pb, 5)
    
    message("\n\nSimulation complete")
    
    # Results list
    
    Output <- list(
      Drug_X = Temp_List_Drug_X,
      Drug_A = Temp_List_Drug_A,
      Drug_B = Temp_List_Drug_B,
      Drug_C = Temp_List_Drug_C,
      Drug_D = Temp_List_Drug_D
    )

    # stop parallel
    
    stopCluster(cl)
    
    # Save results
    
    if (Save_Result) {
      if (!dir.exists(paste("Outputs/Sim_", File_String, sep = ""))) {
        dir.create(paste("Outputs/Sim_", File_String, sep = ""))
      }
      
      saveRDS(Output,
              paste("Outputs/Sim_", File_String, "/Output.rds", sep = ""))
    }
  }
  
  if (Dir_Test == TRUE) {
    Use_Prev_Result_Test <-
      askYesNo(
        "A simulation has previously been run using these input parameters, would you like to load the previous results?"
      )
    
    if (Use_Prev_Result_Test) {
      message("Loading previous simulation")
      
      Output <-
        readRDS(paste("Outputs/Sim_", File_String, "/Output.rds", sep = ""))
      
      message("Results loaded")
      
    }
    
    else{
      # Ask user is they want to save results for next time
      
      Save_Result <-
        askYesNo("Would you like to save the outputs of the simulation for a quick load next time?")
      
      # Remove existing file if new results are to be saved
      
      if(Save_Result){
        file.remove(paste("Outputs/Sim_", File_String, "/Output.rds", sep = ""))
      }
      
      # Set up parallel
      
      cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
      registerDoSNOW(cl)
      
      # Fill list ===================================================
      
      message("Starting simulation")
      pb = txtProgressBar(
        min = 0,
        max = 5,
        initial = 0,
        style = 3
      )
      
      message("\n\nGenerating results for Drug_X")
      # Drug X
      
      Temp_List_Drug_X <- foreach(
        i = 1:Input_MS_Patient_Number,
        .packages = c("purrr"),
        .export = c(ls(.GlobalEnv)),
        .inorder = FALSE
      ) %dopar% {
        func_Paient_Results_Single(
          treat.name = "Drug_X",
          time.horizon = Input_MS_Time_Horizon,
          bl.age =  Out_Patient_Char[[i]][["BL_Age"]],
          treat.len = Input_MS_Treat_Length,
          inj.num = Out_Patient_Char[[i]][["Injection_N"]],
          bl.gfr = Out_Patient_Char[[i]][["GFR"]],
          gfr.decline = Out_Patient_Char[[i]][["GFR_Decline"]],
          treat.effect = Out_Patient_Char[[i]][["Treatment_Effect"]],
          mort.rate = Out_Patient_Char[[i]][["Mortality_Rate"]],
          sex = Out_Patient_Char[[i]][["Sex"]],
          mort.gen = Input_Mort_Gen,
          trae.stroke = Input_TRAE_Stroke,
          trae.mi = Input_TRAE_MI,
          dialysis.start = Input_MS_Dialysis_Threshold,
          cost.trt = Input_RC_Treatment,
          cost.stroke = c(Input_RC_Stroke, Input_RC_Stroke_LT),
          cost.mi = c(Input_RC_MI, Input_RC_MI_LT),
          cost.dialysis = Input_RC_Dialysis,
          discon.rate = c(Input_DC_Cost, Input_DC_Benefits),
          hrqol.info = Out_Patient_Char[[i]][["HRQoL_Info"]][-4]
        )
      }
      
      setTxtProgressBar(pb, 1)
      
      # Drug A
      
      message("\n\nGenerating results for Drug_A")
      
      Temp_List_Drug_A <- foreach(
        i = 1:Input_MS_Patient_Number,
        .packages = c("purrr"),
        .export = c(ls(.GlobalEnv)),
        .inorder = FALSE
      ) %dopar% {
        func_Paient_Results_Single(
          treat.name = "Drug_A",
          time.horizon = Input_MS_Time_Horizon,
          bl.age =  Out_Patient_Char[[i]][["BL_Age"]],
          treat.len = Input_MS_Treat_Length,
          inj.num = Out_Patient_Char[[i]][["Injection_N"]],
          bl.gfr = Out_Patient_Char[[i]][["GFR"]],
          gfr.decline = Out_Patient_Char[[i]][["GFR_Decline"]],
          treat.effect = Out_Patient_Char[[i]][["Treatment_Effect"]],
          mort.rate = Out_Patient_Char[[i]][["Mortality_Rate"]],
          sex = Out_Patient_Char[[i]][["Sex"]],
          mort.gen = Input_Mort_Gen,
          trae.stroke = Input_TRAE_Stroke,
          trae.mi = Input_TRAE_MI,
          dialysis.start = Input_MS_Dialysis_Threshold,
          cost.trt = Input_RC_Treatment,
          cost.stroke = c(Input_RC_Stroke, Input_RC_Stroke_LT),
          cost.mi = c(Input_RC_MI, Input_RC_MI_LT),
          cost.dialysis = Input_RC_Dialysis,
          discon.rate = c(Input_DC_Cost, Input_DC_Benefits),
          hrqol.info = Out_Patient_Char[[i]][["HRQoL_Info"]][-4]
        )
      }
      
      setTxtProgressBar(pb, 2)
      
      # Drug B
      
      message("\n\nGenerating results for Drug_B")
      
      Temp_List_Drug_B <- foreach(
        i = 1:Input_MS_Patient_Number,
        .packages = c("purrr"),
        .export = c(ls(.GlobalEnv)),
        .inorder = FALSE
      ) %dopar% {
        func_Paient_Results_Single(
          treat.name = "Drug_B",
          time.horizon = Input_MS_Time_Horizon,
          bl.age =  Out_Patient_Char[[i]][["BL_Age"]],
          treat.len = Input_MS_Treat_Length,
          inj.num = Out_Patient_Char[[i]][["Injection_N"]],
          bl.gfr = Out_Patient_Char[[i]][["GFR"]],
          gfr.decline = Out_Patient_Char[[i]][["GFR_Decline"]],
          treat.effect = Out_Patient_Char[[i]][["Treatment_Effect"]],
          mort.rate = Out_Patient_Char[[i]][["Mortality_Rate"]],
          sex = Out_Patient_Char[[i]][["Sex"]],
          mort.gen = Input_Mort_Gen,
          trae.stroke = Input_TRAE_Stroke,
          trae.mi = Input_TRAE_MI,
          dialysis.start = Input_MS_Dialysis_Threshold,
          cost.trt = Input_RC_Treatment,
          cost.stroke = c(Input_RC_Stroke, Input_RC_Stroke_LT),
          cost.mi = c(Input_RC_MI, Input_RC_MI_LT),
          cost.dialysis = Input_RC_Dialysis,
          discon.rate = c(Input_DC_Cost, Input_DC_Benefits),
          hrqol.info = Out_Patient_Char[[i]][["HRQoL_Info"]][-4]
        )
      }
      
      setTxtProgressBar(pb, 3)
      
      # Drug C
      
      message("\n\nGenerating results for Drug_C")
      
      Temp_List_Drug_C <- foreach(
        i = 1:Input_MS_Patient_Number,
        .packages = c("purrr"),
        .export = c(ls(.GlobalEnv)),
        .inorder = FALSE
      ) %dopar% {
        func_Paient_Results_Single(
          treat.name = "Drug_C",
          time.horizon = Input_MS_Time_Horizon,
          bl.age =  Out_Patient_Char[[i]][["BL_Age"]],
          treat.len = Input_MS_Treat_Length,
          inj.num = Out_Patient_Char[[i]][["Injection_N"]],
          bl.gfr = Out_Patient_Char[[i]][["GFR"]],
          gfr.decline = Out_Patient_Char[[i]][["GFR_Decline"]],
          treat.effect = Out_Patient_Char[[i]][["Treatment_Effect"]],
          mort.rate = Out_Patient_Char[[i]][["Mortality_Rate"]],
          sex = Out_Patient_Char[[i]][["Sex"]],
          mort.gen = Input_Mort_Gen,
          trae.stroke = Input_TRAE_Stroke,
          trae.mi = Input_TRAE_MI,
          dialysis.start = Input_MS_Dialysis_Threshold,
          cost.trt = Input_RC_Treatment,
          cost.stroke = c(Input_RC_Stroke, Input_RC_Stroke_LT),
          cost.mi = c(Input_RC_MI, Input_RC_MI_LT),
          cost.dialysis = Input_RC_Dialysis,
          discon.rate = c(Input_DC_Cost, Input_DC_Benefits),
          hrqol.info = Out_Patient_Char[[i]][["HRQoL_Info"]][-4]
        )
      }
      
      setTxtProgressBar(pb, 4)
      
      # Drug D
      
      message("\n\nGenerating results for Drug_D")
      
      Temp_List_Drug_D <- foreach(
        i = 1:Input_MS_Patient_Number,
        .packages = c("purrr"),
        .export = c(ls(.GlobalEnv)),
        .inorder = FALSE
      ) %dopar% {
        func_Paient_Results_Single(
          treat.name = "Drug_D",
          time.horizon = Input_MS_Time_Horizon,
          bl.age =  Out_Patient_Char[[i]][["BL_Age"]],
          treat.len = Input_MS_Treat_Length,
          inj.num = Out_Patient_Char[[i]][["Injection_N"]],
          bl.gfr = Out_Patient_Char[[i]][["GFR"]],
          gfr.decline = Out_Patient_Char[[i]][["GFR_Decline"]],
          treat.effect = Out_Patient_Char[[i]][["Treatment_Effect"]],
          mort.rate = Out_Patient_Char[[i]][["Mortality_Rate"]],
          sex = Out_Patient_Char[[i]][["Sex"]],
          mort.gen = Input_Mort_Gen,
          trae.stroke = Input_TRAE_Stroke,
          trae.mi = Input_TRAE_MI,
          dialysis.start = Input_MS_Dialysis_Threshold,
          cost.trt = Input_RC_Treatment,
          cost.stroke = c(Input_RC_Stroke, Input_RC_Stroke_LT),
          cost.mi = c(Input_RC_MI, Input_RC_MI_LT),
          cost.dialysis = Input_RC_Dialysis,
          discon.rate = c(Input_DC_Cost, Input_DC_Benefits),
          hrqol.info = Out_Patient_Char[[i]][["HRQoL_Info"]][-4]
        )
        
      }
      
      setTxtProgressBar(pb, 5)
      
      message("\n\nSimulation complete")
      
      # Results list
      
      Output <- list(
        Drug_X = Temp_List_Drug_X,
        Drug_A = Temp_List_Drug_A,
        Drug_B = Temp_List_Drug_B,
        Drug_C = Temp_List_Drug_C,
        Drug_D = Temp_List_Drug_D
      )
      
      # stop parallel
      
      stopCluster(cl)
      
      # Save results
      
      if (Save_Result) {
        if (!dir.exists(paste("Outputs/Sim_", File_String, sep = ""))) {
          dir.create(paste("Outputs/Sim_", File_String, sep = ""))
        }
        
        saveRDS(Output,
                paste("Outputs/Sim_", File_String, "/Output.rds", sep = ""))
        
      }
    }
  }
  
  return(Output)
  
}


## Condense results ==================================================================================================================================================

# Function to extract key results to display

func_Condense_Results <- function(...){
  
  message("Condensing results")
  
  Output <- list(
    Max = list( # Find the max cost and HRQoL for each patient
      Cost = map(.x = Out_Sim_Results, ~map_dbl(.x = .x, ~max(.x[["V_Cost_Disc_Total"]], na.rm = TRUE))),
      HRQoL = map(.x = Out_Sim_Results, ~map_dbl(.x = .x, ~max(.x[["V_HRQoL_Disc"]], na.rm = TRUE)))
    ),
    Mean = list(
      Cost = 0,
      HRQoL = 0
    ),
    Rolling_Mean = list(
      Cost = 0,
      HRQoL = 0
    )
  )
  
  Output[["Mean"]] <- map(.x = Output[["Max"]], ~map_dbl(.x = .x, ~mean(.x, na.rm = TRUE))) # Find the average for each treatment arm
  
  Output[["Rolling_Mean"]][["Cost"]] <- map(.x = Output[["Max"]][["Cost"]], ~cummean(.x))
  Output[["Rolling_Mean"]][["HRQoL"]] <- map(.x = Output[["Max"]][["HRQoL"]], ~cummean(.x))
  
  message("Results condensed")
  
  return(Output)
  
}


## Fully incremental cost effectiveness analysis ====================================================================================================================

# Function calculates all of the results and produces a table
# User can save the output by using save.table = TRUE

func_CE_Result <- function(save.table = FALSE){
  
  # Condusct fully incremental cost-effectivness analysis
  
  ICERS <- calculate_icers(c(Out_Sim_Summary[["Mean"]][["Cost"]][["Drug_X"]], Out_Sim_Summary[["Mean"]][["Cost"]][["Drug_A"]], Out_Sim_Summary[["Mean"]][["Cost"]][["Drug_B"]], Out_Sim_Summary[["Mean"]][["Cost"]][["Drug_C"]], Out_Sim_Summary[["Mean"]][["Cost"]][["Drug_D"]]
  ),
  c(Out_Sim_Summary[["Mean"]][["HRQoL"]][["Drug_X"]], Out_Sim_Summary[["Mean"]][["HRQoL"]][["Drug_A"]], Out_Sim_Summary[["Mean"]][["HRQoL"]][["Drug_B"]], Out_Sim_Summary[["Mean"]][["HRQoL"]][["Drug_C"]], Out_Sim_Summary[["Mean"]][["HRQoL"]][["Drug_D"]]
  ),
  c("Drug_X", "Drug_A", "Drug_B", "Drug_C", "Drug_D")
  )
  
  # Recode status variable to make it clearer
  
  ICERS <- ICERS %>% 
    mutate(Status = recode(Status, ND = "Non-dominated"),
           Status = recode(Status, D = "Dominated"),
           Status = recode(Status, ED = "Extended dominated"))
  
  # Output a table
  
  Output <- ICERS %>% 
    mutate(Cost = comma(Cost),
           Effect = round(Effect, 2),
           Inc_Cost = ifelse(is.na(Inc_Cost), NA, comma(Inc_Cost)),
           Inc_Effect = ifelse(is.na(Inc_Effect), NA, round(Inc_Effect, 2)),
           ICER = ifelse(is.na(ICER), NA, comma(ICER))) %>%
    kable()
  
  # Save table as excel file if required
  
  if(save.table){
    
    # Create directory if it doesn't exist
    
    if (!dir.exists(paste("Outputs/Sim_", Input_Model_Code, sep = ""))) {
      dir.create(paste("Outputs/Sim_", Input_Model_Code, sep = ""))
    }
    
    ICERS %>% 
      mutate(Cost = comma(Cost),
             Effect = round(Effect, 2),
             Inc_Cost = ifelse(is.na(Inc_Cost), NA, comma(Inc_Cost)),
             Inc_Effect = ifelse(is.na(Inc_Effect), NA, round(Inc_Effect, 2)),
             ICER = ifelse(is.na(ICER), NA, comma(ICER))) %>%
      write_csv(paste("Outputs/Sim_", Input_Model_Code, "/Result_Table.csv", sep = ""))
  }
  
  # Return table
  
  return(Output)
  
}


## Create a CEP ==============================================================================================================================

func_Plot_CEP <- function(treat.name,
                          treat.name2,
                          heat.map = FALSE,
                          y.limits = NULL,
                          y.breaks = NULL,
                          x.limits = NULL,
                          x.breaks = NULL,
                          txt.size = 12,
                          currency = "£") {
  
  ce_df <- data.frame("Cost" = Out_Sim_Summary[["Max"]][["Cost"]][[treat.name]] - Out_Sim_Summary[["Max"]][["Cost"]][[treat.name2]],
                      "Effectiveness" = Out_Sim_Summary[["Max"]][["HRQoL"]][[treat.name]]- Out_Sim_Summary[["Max"]][["HRQoL"]][[treat.name2]])
  
  Output <- ce_df %>% 
    ggplot(aes(x = Effectiveness, y = Cost))+
    geom_point(col = "#6FB4F4")+
    geom_point(aes(x = mean(Effectiveness), y = mean(Cost)), col = "red", size = 3)+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    geom_abline(intercept = 0, slope = Input_MS_WPay, lty = 2, size = 1)+
    scale_x_continuous("QALY Difference", expand=c(0,0), limits = x.limits, breaks = x.breaks)+
    scale_y_continuous("Cost Difference (\u00A3)", expand=c(0,0), label = comma, limits = y.limits, breaks = y.breaks)+
    ggtitle(paste("Cost-effectiveness plane for", treat.name, "versus", treat.name2, sep = " "))+
    theme_classic(base_size = txt.size)
  
  # Add a heat map where there are lots of points
  
  if(heat.map){
    Output <- Output+
      stat_density2d(aes(alpha=..level..),geom='polygon',fill="#3e678c", colour='#084e8c', show.legend = FALSE, bins=10)
  }
  
  return(Output)
  
}


## Calculate outcome function =====================================================================================================================

# Taken from dampack package so it can be used in the modified plot ceac function

calculate_outcome <- function(outcome = c("nhb", "nmb", "eff", "cost", "nhb_loss",
                                          "nmb_loss", "nhb_loss_voi", "nmb_loss_voi"),
                              cost, effect, wtp) {
  outcome <- match.arg(outcome)
  n_sim <- nrow(cost)
  if (outcome == "eff") {
    y <- effect
  } else if (outcome == "cost") {
    y <- cost
  } else {
    if (is.null(wtp)) {
      # the call. = FALSE makes the error message more clear
      stop("wtp must be provided for NHB and NMB",  call. = FALSE)
    }
    if (is.null(cost)) {
      stop("must provide cost for NHB and NMB.",  call. = FALSE)
    }
    if (outcome == "nhb") {
      y <- effect - cost / wtp
    }
    if (outcome == "nmb") {
      y <- effect * wtp - cost
    }
    if (outcome == "nhb_loss" | outcome == "nmb_loss") {
      if (outcome == "nhb_loss") {
        net_outcome <- "nhb"
      }
      if (outcome == "nmb_loss") {
        net_outcome <- "nmb"
      }
      netben <- calculate_outcome(net_outcome, cost, effect, wtp)
      max_str_rowwise <- max.col(netben)
      y <-  netben[cbind(1:n_sim, max_str_rowwise)] - netben
    }
    if (outcome == "nhb_loss_voi" | outcome == "nmb_loss_voi") {
      if (outcome == "nhb_loss_voi") {
        net_outcome <- "nhb"
      }
      if (outcome == "nmb_loss_voi") {
        net_outcome <- "nmb"
      }
      netben <- calculate_outcome(net_outcome, cost, effect, wtp)
      max_str <- which.max(colMeans(netben))
      y <- netben - netben[cbind(1:n_sim), max_str]
    }
  }
  return(y)
}


## Plot CEAC =========================================================================================================================================

# Adapted from dampack

func_Plot_CEAC <- function(treat.names, wtp, interactive = FALSE,
                           y.limits = c(0,1),
                           y.breaks = seq(0,1,0.1),
                           x.limits = c(0,100),
                           x.breaks = seq(0,100,10),
                           txt.size = 12,
                           currency = "\u00A3"){
  
  # define needed variables
  strategies <- treat.names
  n_strategies <- length(treat.names)
  effectiveness <- as.data.frame(Out_Sim_Summary[["Max"]][["HRQoL"]])
  cost <- as.data.frame(Out_Sim_Summary[["Max"]][["Cost"]])
  n_sim <- Input_MS_Patient_Number
  
  # number of willingness to pay thresholds
  n_wtps <- length(wtp)
  
  # matrix to store probability optimal for each strategy
  cea <- matrix(0, nrow = n_wtps, ncol = n_strategies)
  colnames(cea) <- strategies
  
  # vector to store strategy at the cost-effectiveness acceptability frontier
  frontv <- rep(0, n_wtps)
  
  for (l in 1:n_wtps) {
    # calculate net monetary benefit at wtp[l]
    lth_wtp <- wtp[l]
    nmb <-  calculate_outcome("nmb", cost, effectiveness, lth_wtp)
    
    # find the distribution of optimal strategies
    max.nmb <- max.col(nmb)
    opt <- table(max.nmb)
    cea[l, as.numeric(names(opt))] <- opt / n_sim
    
    # calculate point on CEAF
    # the strategy with the highest expected nmb
    frontv[l] <- which.max(colMeans(nmb))
  }
  
  # make cea df
  cea_df <- data.frame(wtp, cea, strategies[frontv],
                       stringsAsFactors = FALSE)
  colnames(cea_df) <- c("WTP", strategies, "fstrat")
  
  # make ceac df
  ceac <- melt(cea_df, id.vars = c("WTP", "fstrat"),
               variable.name = "Strategy", value.name = "Proportion")
  
  # replace factors with strings (melt creates factors)
  ceac$Strategy <- as.character(ceac$Strategy)
  
  # boolean for on frontier or not
  ceac$On_Frontier <- (ceac$fstrat == ceac$Strategy)
  
  # drop fstrat column
  ceac$fstrat <- NULL
  
  # order by WTP
  ceac <- ceac[order(ceac$WTP), ]
  
  # remove rownames
  rownames(ceac) <- NULL
  
  # plot =====================================================================================================================================
  
  Output <- ceac %>% 
    ggplot(aes(x = WTP/1000, y = Proportion, colour = Strategy))+
    geom_line() +
    scale_x_continuous(paste("Willingness to Pay (Thousand ", currency, " / QALY)", sep = ""), expand=c(0,0), limits = x.limits, breaks = x.breaks)+
    scale_y_continuous("Probability of being Cost-Effective", expand=c(0,0), limits = y.limits, breaks = y.breaks)+
    theme_classic(base_size = txt.size)
  
  if(interactive){
    Output <- ggplotly(Output)
  }
  
  return(Output)
  
}


## Test model stability =====================================================================================================================

# This is a key function as it allows you to visually see if the model has reached stability

func_Model_Stability <- function(type = c("Cost", "HRQoL"), interactive = FALSE,
                                 y.limits = NULL,
                                 y.breaks = NULL,
                                 x.breaks = seq(0,Input_MS_Patient_Number,Input_MS_Patient_Number/10),
                                 txt.size = 12){
  
  # Extract the required data into a dataframe
  
  Temp_Data <- data.frame(as.data.frame(Out_Sim_Summary$Rolling_Mean[[type]]), N = 1:Input_MS_Patient_Number)
  
  Temp_Data <- melt(Temp_Data, id.vars = "N",
                    variable.name = "Strategy", value.name = "Average")
  
  # plot
  
  Output <- Temp_Data %>% 
    ggplot(aes(x = N, y = Average, col = Strategy))+
    geom_line() +
    scale_x_continuous(paste("Number of simulations", sep = ""), expand=c(0,0), limits = c(0, Input_MS_Patient_Number), breaks = x.breaks)+
    scale_y_continuous(paste("Rolling Average (", type, ")", sep = ""), expand=c(0,0), limits = y.limits, breaks = y.breaks, label = comma)+
    theme_classic(base_size = txt.size)
  
  if(interactive){
    Output <- ggplotly(Output)
  }
  
  return(Output)
  
  
}
