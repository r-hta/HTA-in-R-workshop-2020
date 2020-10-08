

    library(dplyr)
    library(tidyr)
    library(expss)                                                                            # vlookup_df function
    library(MASS)                                                                             # mvrnorm function for simulation
    library(ggplot2)                                                                          # plot the PSA results
    
    
    
    ###############################################################################
    ####################### Part 0: Parameter set up ##############################
    ###############################################################################

    ### 0.1 user input parameters
      age <- 60
      SIMD <- 60.8                                                                            # Scottish Index of Multiple Deprivation (continuous, range: ?)
      Diabetes <- 0                                                                           # diabetes: 0 - no diabetes; 1 - with diabetes
      FH <- 0                                                                                 # family history: 0 - no history ; 1 - with history
      CPD <- 20                                                                               # No. of cigarette per day (continuous)
      SBP <- 160                                                                              # Systolic blood pressure
      TC <- 7                                                                                 # Total cholestorol 
      HDL <- 1
      sex <-1                                                                                 # male. sex <-2 female
      
      tm <- 100                                                                               # As consistent with the Excel model - time_horizon set as 115 cycles. Alternatively, if setting tm to be bounded by age 100: tm <- 100 - age
      disc <-0.035                                                                            # disc <-0 for undiscounted

      statin_ldleffect <- 0.26                                                                # treatment effects on ldl and hdl
      statin_hdleffect <- 1.04
      statin_cost <- 13
      statin_disu <-0.001
      
      
    
    ### 0.2 load model parameters from survival analysis, cost analysis and utility. 
      source("CVD_parameters.R",local = TRUE)                                                 # This script contains two options (EXCEL or script) of loading data. If run the full script, data are loaded from EXCEL first and then the objects are re-written by the script. But they are the same anyway.
      
      if(sex == 1){                                                                           # set a switch to choose to use whether male or female parameters.
        first_event_coef <- first_event_coef_m
        post_event_coef <- post_event_coef_m
        secondary_event_inputs <- secondary_event_inputs_m
        coef_c1 <- coef_c1_m
        coef_c2 <- coef_c2_m
      }else{
        first_event_coef <- first_event_coef_f
        post_event_coef <- post_event_coef_f
        secondary_event_inputs <- secondary_event_inputs_f
        coef_c1 <- coef_c1_f
        coef_c2 <- coef_c2_f
      }
    
      
      ###############################################################################
      ####################### Part I: Life expectancy ###############################
      ###############################################################################


    model <- function(age,SIMD,Diabetes,FH,CPD,SBP,TC,HDL,sex,first_event_coef,post_event_coef,secondary_event_inputs, coef_c1,coef_c2, disc, tm) {                               
    
      ### 1.0 additional adjustment for the coefficients
      
      calibration_f1 <- 0.96
      multi<-ifelse(sex==1, 0.99, 1.05)
      first_event_coef[9,]<-first_event_coef[9,]*multi                                      # Constant
      
      if (sex==1) {
        first_event_coef[7:8,4] <- 0                                                        # coefficient adjustment for TC and HDL for nonCVDdeath first event for male
      } else {
        first_event_coef[7,2]<-0                                                            # coefficient adjustment for TC CBVD and nonCVDdeath first events for female
        first_event_coef[7,4]<-0
      }
      


    ### 1.1 First equation: Hazard (transition probability), Cumulative incidence (trace)
    
    ## Linear predictors
      riskf_v <- c(age,SIMD,Diabetes,FH,CPD,SBP,TC,HDL)
      riskf_v2 <- c(age,SIMD,FH)
      
      linpred <- rep(0,4)
      for (num in 1:4) {
        linpred[num]<- (sum(riskf_v*first_event_coef[1:8,num])+first_event_coef[9,num])*calibration_f1
      }
    
    
    ## Cycle and age vectors 
      cycle <- seq(from = 1, to = tm,by = 1)                                                  # From cycle 1 to end, not including cycle 0
      age_cyc <- age + cycle
    
    
    ## Hazards of each first event
      hazard <- matrix(nrow=tm, ncol=4,0)
      for (num in 1:4) {
        hazard[,num]<- exp(linpred[num])*exp(first_event_coef[10,num]*cycle)
      }
    
    
    ## Cumulative incidence
      CumInc <- matrix(data = 0,nrow = tm,ncol = 5)                                           # Set up. 1st column is CVD free
      CumInc[1,] <- c((1-sum(hazard[1,])),hazard[1,])
      
      for (t in 2:tm){
        CumInc[t,2] <-  CumInc[t-1,2]+CumInc[t-1,1]*hazard[t,1]
        CumInc[t,3] <-  CumInc[t-1,3]+CumInc[t-1,1]*hazard[t,2]
        CumInc[t,4] <-  CumInc[t-1,4]+CumInc[t-1,1]*hazard[t,3]
        CumInc[t,5] <-  CumInc[t-1,5]+CumInc[t-1,1]*hazard[t,4]
        CumInc[t,1] <- ifelse(CumInc[t-1,1]<0.001, 0, (1- sum(hazard[t,]))*CumInc[t-1,1])     # The fudge is to avoid negative numbers for alive CVD free
      }
      
      
    ## Probability of getting an event at each cycle
      CumInc_temp <- rbind(rep(0,4), CumInc[-nrow(CumInc),2:5])                               # Create a matrix to enable the cumulative events deduction.
      CumInc_all<- rowSums(CumInc[,2:5])
      prob_event_yr <- CumInc[,2:5]- CumInc_temp
      prob_event_yr[CumInc_all>0.998]<-0                                                      # Fudge. if the total hazard >0.998, then change all the probability to 0, assuming all people have died.
    
      
    ### 1.2 - 2nd equation: Remaining life years for the non-fatal first events 
      age_auc<- age:(age+tm-1)                                                                # 115 years time horizon including time 0, so finish at 114+age.
      agem<-t(replicate(tm, age_auc)) 
      cyclem<-matrix(rep(cycle, tm), ncol=tm)
      discym <- agem-age
      
      
      
    ## Post event survival
      post_LE <- matrix(0, nrow=tm, ncol = 2)
      for (num in 1:2)  {                                                                     # For both CHD and CBVD
        Post_surv <- exp(-exp(post_event_coef[1,num]*(agem)+post_event_coef[2,num]*SIMD+post_event_coef[3,num]*FH+post_event_coef[4,num])*
                                  (1/post_event_coef[5,num])*(exp(post_event_coef[5,num]*cyclem)-1))
        Post_surv_sub1 <- rbind(rep(1,tm), Post_surv)                                         # Add one row of 1 at the top for calculation convenience
        Post_surv_sub2 <- rbind(Post_surv, rep(0, tm))                                        # Add one row of 0 at the bottom for calculation convenience
        Post_surv_AUC <- (Post_surv_sub1+Post_surv_sub2)/2                                    # Trapezoid
        Post_surv_AUC_und<-Post_surv_AUC[-(tm),]                                              # Delete the last row for calculation convenience.
        if (num==1) {Post_surv_AUC_CHD_und<-Post_surv_AUC_und }                               # Keep the result for IHD and CBVD for the survival condition for costs.
        Post_surv_AUC <- Post_surv_AUC_und / (1+disc)^(cyclem+discym)
        post_LE[,num] <- colSums(Post_surv_AUC)                                               # Remaining life expectancy
      }
        Post_surv_AUC_CBVD_und<-Post_surv_AUC_und
        
        
    ## All survival since start of the model.
        discy <- 1/(1+disc)^cycle
        preeventLE <- NULL
        for (i in 1: tm) {
          preeventLE[i] <- sum(discy[1:i])
        }
        LY_add <- matrix(data = 0,nrow = tm,ncol = 4)                                         # Calculate additional LY (undiscounted)
        LY_add[,1] <- post_LE[1:tm,1] + preeventLE
        LY_add[,2] <- post_LE[1:tm,2] + preeventLE
        LY_add[,3] <- preeventLE
        LY_add[,4] <- preeventLE

    
    ## 1.3 - calculate the weighted LE per period.
        Total_LY = prob_event_yr*LY_add
      
        
    ## 1.4 Results ---------->>  Total LE
        sum(Total_LY)+age                                                                     # This should be cross checked with 'Male - Parameters' - cell B3.
    
        
        
        
        
  
    ###############################################################################
    ################################# Part II: Cost ###############################
    ###############################################################################  
  
    ### 2.1 Cost setup                                                                        # 't1, t2 all 1-100, whereas the cycles are 1-115. At the later cycles, t1, t2 are assumed to be 0 in the excel model.'
    m0<-matrix(0, nrow=15, ncol=5)
    coef_ct<-data.matrix(coef_c2, rownames.force = NA)
    coef_ct<-rbind(coef_ct,m0)
    coef_ct<-coef_ct[1:tm,]                                                                   # take the rows up to the defined time horizon
    
    
    ## 2.2 Calculate the cost for pre event
    c_preCHD_und <- sum(riskf_v2*coef_c1[3:5, 1]) + coef_c1[6,1]+coef_c1[1,1]*coef_ct[,2]+coef_c1[2,1]*coef_ct[,5]
    c_preCBVD_und <- sum(riskf_v2*coef_c1[3:5, 2]) + coef_c1[6,2]+coef_c1[1,2]*coef_ct[,2]+coef_c1[2,2]*coef_ct[,5]
    
    ## Incorporating the discounting factor
    c_preCHD_disc<-c_preCHD_und*discy
    c_preCBVD_disc <- c_preCBVD_und*discy
    
    ## Cumulative sums
    c_preCHD_disc_cul <-cumsum(c_preCHD_disc)                                                 # If disc=0, then this will give the result of undiscounted.
    c_preCBVD_disc_cul <- cumsum(c_preCBVD_disc)
    

    ### 2.3 Calculate the cost after CHD and DBVD non fatal events
    
    ## Create a parameter matrix for t1, t2, and age (these parameters change vertically or horizontally in the matrix)
    mc_ti1<-matrix(rep(coef_ct[,2], tm), ncol=tm)
    mc_postCHD_ti2<-matrix(rep(coef_ct[,3], tm), ncol=tm)
    mc_postCBVD_ti2<-matrix(rep(coef_ct[,4], tm), ncol=tm)
    m_age<- age:(age+tm-1)
    m_agem<-t(replicate(tm, m_age))
    
    ## Calculate the cost
    c_postCHD<-sum(riskf_v2[2:3]*coef_c1[4:5, 3]) + coef_c1[3,3]*m_agem+ coef_c1[6,3]+coef_c1[1,3]*mc_ti1+coef_c1[2,3]*mc_postCHD_ti2
    c_postCBVD<-sum(riskf_v2[2:3]*coef_c1[4:5, 4]) + coef_c1[3,4]*m_agem+ coef_c1[6,4]+coef_c1[1,4]*mc_ti1+coef_c1[2,4]*mc_postCBVD_ti2
    
    ## Kaplan-Meier sample average (KMSA): Cost conditional on survival at each year after the event
    c_postCHD_sur <- c_postCHD*Post_surv_AUC_CHD_und
    c_postCBVD_sur<- c_postCBVD*Post_surv_AUC_CBVD_und
    
    ## Incorprating the discount factor
    c_postCHD_disc <- c_postCHD_sur / (1+disc)^(cyclem+discym)
    c_postCBVD_disc <- c_postCBVD_sur / (1+disc)^(cyclem+discym)
    
    
    ## Calculate lifetime cost
    c_postCHD_life<- colSums(c_postCHD_disc)
    c_postCBVD_life<-colSums(c_postCBVD_disc)
    
    
    ### 2.4 Calculate the cost prior to CVD death and non-CVD death  
    c_fatalcvd_ann<-sum(riskf_v2*coef_c1[3:5,5])+coef_c1[6,5]+coef_c1[1,5]*coef_ct[,2]+coef_c1[2,5]*coef_ct[,5]
    c_fatalnoncvd_ann<-sum(riskf_v2*coef_c1[3:5,6])+coef_c1[6,6]+coef_c1[1,6]*coef_ct[,2]+coef_c1[2,6]*coef_ct[,5]
    c_fatalcvd <- cumsum(c_fatalcvd_ann)
    c_fatalnoncvd <- cumsum(c_fatalnoncvd_ann)
    
    c_fatalcvd_disc <-c_fatalcvd*discy
    c_fatalnoncvd_disc <- c_fatalnoncvd*discy
  
    
    ### 2.5 Calculate the total costs.
    c_chd<-c_preCHD_disc_cul + c_postCHD_life
    c_cbvd<-c_preCBVD_disc_cul + c_postCBVD_life
    mc <-cbind(c_chd, c_cbvd, c_fatalcvd_disc,c_fatalnoncvd_disc)
    
    
    ### 2.6 weighted cost (use probability for each event at each year (equation 1))
    cw<- prob_event_yr*mc
    
    ### 2.7 Results ---------->>  Total cost without intervention                             # This should be copmared to 'male parameters' - B5
    sum(cw)                                                                                   # if disc=0, this will give undiscounted cost. if disc=0.035, this will provide discounted cost.
    
    

    
    
    
    
    ###############################################################################
    ################################ Part III: QALYs ##############################
    ############################################################################### 

    
    ### 3.1 set up event utility decrements                                                   
    if (sex==1) {
      dec <- dec_m
    } else {
      dec <- dec_f
    }
    
    ### 3.2 Combined decrements for secondary events                                          # Set up a loop for each secondary event, and create a combined secondary event decrement matrix for each first event (CHD, CBVD) happended at age 60, 61, 63 etc.
    secondary_event_CHD <- c("SecondaryCHD_CHD","SecondaryCBVD_CHD","SecondaryHF_CHD","SecondaryPAD_CHD","Secondaryotherheart_CHD")
    secondary_event_CBVD <- c("SecondaryCHD_CBVD","SecondaryCBVD_CBVD","SecondaryHF_CBVD","SecondaryPAD_CBVD","Secondaryotherheart_CBVD")
    decrements_comb_CHD <- matrix(nrow = tm, ncol = tm)
    decrements_comb_CBVD <- matrix(nrow = tm, ncol = tm)
    
    for (j in 0: (tm-1)) {
          decrements_mat_CHD <- matrix(nrow = tm,ncol = length(dec),data = 0)
          decrements_mat_CBVD <- matrix(nrow = tm,ncol = length(dec),data = 0)
          
          # post CHD
          for (i in 1: length(dec)){
            eqn <- as.numeric(secondary_event_inputs[secondary_event_inputs[,1] == secondary_event_CHD[i],2:7])
            lp <- eqn[6]+eqn[5]*FH+eqn[4]*SIMD+eqn[3]*(age+j)+eqn[2]*coef_ct[,3]+eqn[1]*coef_ct[,2]
            prop <- pnorm(lp,0,1)
            decrement <- prop*dec[i]
            decrements_mat_CHD[,i]<- decrement
          }
      
          # post CBVD 
          for (i in 1: 4){
            eqn <- as.numeric(secondary_event_inputs[secondary_event_inputs[,1] == secondary_event_CBVD[i],2:7])
            lp <- eqn[6]+eqn[5]*FH+eqn[4]*SIMD+eqn[3]*(age+j)+eqn[2]*coef_ct[,4]+eqn[1]*coef_ct[,2]
            prop <- pnorm(lp,0,1)
            decrement <- prop*dec[i]
            decrements_mat_CBVD[,i]<- decrement
          }
      
          # sum decrement columns
          decrements_comb_CHD[,j+1] <- rowSums(decrements_mat_CHD)
          decrements_comb_CBVD[,j+1] <- rowSums(decrements_mat_CBVD)
    }
    

    ### 3.3 Sex specific utility values norm table                                            # This could be moved to the data input script.
    ut1 <- 1:(120+tm+tm)                                                                      # Because the model have a maximum of 115 cycles, and assuming the largest age input into the model is 120.
    if (sex <-1) {
      ut2_1 <- rep(0.831, times=24)
      ut2_2<-rep(c(0.8230,0.8200,0.8060,0.8010,0.7880,0.7740), each=10)
      ut2_3 <- rep(0.7740, times= 120+tm+tm-24-6*10) 
    } else {
      ut2_1 <- rep(0.809, times=24)
      ut2_2<-rep(c(0.8090,0.8020,0.7850,0.7870,0.7770, 0.7210), each=10)
      ut2_3 <- rep(0.7210, times= 120+tm+tm-24-6*10)
    }
    ut2 <- c(ut2_1, ut2_2, ut2_3)
    dict=data.frame(index=ut1, ut=ut2)                                                        # Look up dictionary
    
    
    ### 3.4 vlookup matching norm values
    agecyclem <- cyclem + m_agem                                                              # age matrix to look up
    m_utility_norm <- matrix(0,nrow = tm, ncol = tm)
    for (j in 1:tm) {
       a <- vlookup_df(agecyclem[,j], dict, result_column='ut')
       m_utility_norm[,j]<-as.vector(a[,1])                                                   # the result of vlookup_df is a dataframe so used 'as.vector'.
    }

    
    ### 3.5 remaining QALY after first event
    CHD_QALY_und <- Post_surv_AUC_CHD_und[1:tm,1:tm]*(m_utility_norm-(MI__Dec+decrements_comb_CHD))
    CBVD_QALY_und <- Post_surv_AUC_CBVD_und[1:tm,1:tm]*(m_utility_norm-(Stroke__Dec+decrements_comb_CBVD))
  
    CHD_QALY_dis <- CHD_QALY_und/((1+disc)^(cyclem[1:tm,1:tm]+discym[1:tm,1:tm]))             # can tidy from the beginning as the last column is not used in the discym and cyclem arrays.
    CBVD_QALY_dis <- CBVD_QALY_und/((1+disc)^(cyclem[1:tm,1:tm]+discym[1:tm,1:tm]))
    
    QALY_postCHD<- colSums(CHD_QALY_dis)
    QALY_postCBVD<-colSums(CBVD_QALY_dis)
    
    
    ### 3.6 QALY prior to events                                                              # i.e. QALY between starting age and age of event. this is also the total QALY for CVD fatal first event
    preeventQALY_und <- vlookup_df(age_cyc[1:tm], dict, result_column='ut')
    preeventQALY_und <- as.vector(preeventQALY_und[,1])
    preeventQALY_dis <- preeventQALY_und/((1+disc)^cycle)
    preeventQALY_dis_total <- NULL
    for (i in 1: length(cycle)) {
      preeventQALY_dis_total[i] <- sum(preeventQALY_dis[1:i])
    }
    
    
    ### 3.7 Total QALY
    QALY_CHD_total <- preeventQALY_dis_total + QALY_postCHD                                   # CHD
    QALY_CBVD_total <- preeventQALY_dis_total + QALY_postCBVD                                 # CBVD
    preeventQALY_dis_total                                                                    # fatal events
    mQALY <- cbind(QALY_CHD_total,QALY_CBVD_total, preeventQALY_dis_total, preeventQALY_dis_total)
    
    
    ### 3.8  Weighted QALY                                                                    # use a matrix of probability for each event at each year (equation 1), same for LY and COST * corresponding lifetime QALYs associated with that probability.
    QALYw<- prob_event_yr*mQALY
    
    
    ### 3.9 --------->> results                                                               # calculate the sum weighted QALY (undiscounted). this should be copmared to 'male parameters' - D3
    sum(QALYw)+age
    
    
    #<<<<>>>># Model output
    list <- list("LY"= sum(Total_LY)+age, "Cost" = sum(cw), "QALY" = sum(QALYw)+age, "prob_event_yr" = prob_event_yr)
    return(list)
    }
    
    
    ### I-III: ---------->> Base case results overall 

    # No intervention
    base_nointv <- model(age,SIMD,Diabetes,FH,CPD,SBP,TC,HDL,sex,first_event_coef,post_event_coef, secondary_event_inputs, coef_c1,coef_c2,disc, tm)   # discounted
    base_nointv$LY
    base_nointv$Cost
    base_nointv$QALY
    
    # Statin effect
    TC_statin <-TC-(TC-HDL)*statin_ldleffect
    HDL_statin <- HDL*statin_hdleffect
    
    base_intv <- model(age,SIMD,Diabetes,FH,CPD,SBP,TC=TC_statin,HDL=HDL_statin,sex,first_event_coef,post_event_coef, secondary_event_inputs, coef_c1,coef_c2,disc, tm)   # discounted
    base_intv$LY
    base_intv$Cost
    base_intv$QALY
    
    # Statin cost
    cycle <- seq(from = 1, to = tm,by = 1) 
    statin_anu <- statin_cost / (1+disc)^cycle
    statin_anucum <- cumsum(statin_anu)
    statin_anucum_sur <- statin_anucum * rowSums(base_intv$prob_event_yr)                               # take statin every year until the first event
    
    # STATIN cost results
    sum(statin_anucum_sur)
    
    # STATIN disutility
    statin_disu_anu <- statin_disu / (1+disc)^cycle
    statin_disu_anucum <- cumsum(statin_disu_anu)
    statin_disu_anucum_sur <- statin_disu_anucum * rowSums(base_intv$prob_event_yr)  
    
    # STATA disutility results
    sum(statin_disu_anucum_sur)
    
    # Intervention overall
    sum(statin_anucum_sur) + base_intv$Cost
    base_intv$QALY - sum(statin_disu_anucum_sur)
    
    # Incremental effect
    incre_cost <- sum(statin_anucum_sur) + base_intv$Cost - base_nointv$Cost
    incre_cost
    incre_QALY <- base_intv$QALY - sum(statin_disu_anucum_sur) - base_nointv$QALY
    incre_QALY
    
    
    
    ############################################################################
    ################################ Part IV: PSA ##############################  
    ############################################################################ 
    
    # Monte Carlo Simulation incorporating variance/covariance betweent the parameters
    # Another benefit of R - packages can help you so that you don't have to write the formula yourself, e.g. decomposition.
    # note: The matrix decomposition in the mvrnorm function is done via eigen. authors' note: "although a Choleski decomposition might be faster, the eigendecomposition is stabler."
    # https://www.rdocumentation.org/packages/MASS/versions/7.3-0/topics/mvrnorm
    # Here used 'mvrnorm' from the 'MASS' package. But it could also use the mvrnorm from the rockchalk package, which makes one small modification to 
    # (cont.) facilitate replication of random samples. The aim is to make sure that, after the seed is reset, the first rows of generated data are identical no matter what value is chosen for n. 
    
    
    ### 4.0 Check the data: IMPORTANT: make sure 'first_event_coef' and 'post_event_coef' stores the base case coefficients because the PSA (4.4) will replace these coefficients.
    if(sex == 1){                                                                           
      first_event_coef <- first_event_coef_m
      post_event_coef <- post_event_coef_m
    }else{
      first_event_coef <- first_event_coef_f
      post_event_coef <- post_event_coef_f
    }
    
    
    ### 4.1 Set the No. of simulations
    nrep <- 1000
    
    
    ### 4.2 Set the variance/covariance matrix according to sex
    if(sex == 1){
      Varcov_firstE_CHD <- Varcov_firstE_CHD_m
      Varcov_firstE_CBVD <- Varcov_firstE_CBVD_m
      Varcov_fatal_CVD <- Varcov_fatal_CVD_m
      Varcov_fatal_nonCVD <- Varcov_fatal_nonCVD_m
      Varcov_surv_postCHD <- Varcov_surv_postCHD_m
      Varcov_surv_postCBVD <- Varcov_surv_postCBVD_m
    }else{
      Varcov_firstE_CHD <- Varcov_firstE_CHD_f
      Varcov_firstE_CBVD <- Varcov_firstE_CBVD_f
      Varcov_fatal_CVD <- Varcov_fatal_CVD_f
      Varcov_fatal_nonCVD <- Varcov_fatal_nonCVD_f
      Varcov_surv_postCHD <- Varcov_surv_postCHD_f
      Varcov_surv_postCBVD <- Varcov_surv_postCBVD_f
    }

    
    ## 4.3 Generate draws  
    mu_firstE_CHD <- first_event_coef[,1]                                                   # get the mean coefficient from the base case data
    set.seed(12345)
    psa_pm_firstE_CHD <- MASS::mvrnorm(n=nrep, mu = mu_firstE_CHD, Sigma = Varcov_firstE_CHD)      
    
    mu_firstE_CBVD <- first_event_coef[,2]
    set.seed(12345)
    psa_pm_firstE_CBVD <- MASS::mvrnorm(n=nrep, mu = mu_firstE_CBVD, Sigma = Varcov_firstE_CBVD)
    
    mu_fatal_CVD <- first_event_coef[,3]
    set.seed(12345)
    psa_pm_fatal_CVD <- MASS::mvrnorm(n=nrep, mu = mu_fatal_CVD, Sigma = Varcov_fatal_CVD)
    
    mu_fatal_nonCVD <- first_event_coef[,4]
    set.seed(12345)
    psa_pm_fatal_nonCVD <- MASS::mvrnorm(n=nrep, mu = mu_fatal_nonCVD, Sigma = Varcov_fatal_nonCVD)
    
    mu_surv_postCHD <- post_event_coef[,1]
    set.seed(12345)
    psa_pm_surv_postCHD <- MASS::mvrnorm(n=nrep, mu = mu_surv_postCHD, Sigma = Varcov_surv_postCHD)
    
    mu_surv_postCBVD <- post_event_coef[,2]
    set.seed(12345)
    psa_pm_surv_postCBVD <- MASS::mvrnorm(n=nrep, mu = mu_surv_postCBVD, Sigma = Varcov_surv_postCBVD)
    
    
    ### 4.4 set up the PSA loop
    psa_result <- matrix(nrow = nrep, ncol = 10)
    colnames(psa_result) <- c("LY__nointv", "Cost__nointv", "QALY__nointv","LY_intv", "Cost_intv", "QALY_intv", 
                              "Cost_intv_add_statin", "disQALY_intv_statin", "Incre_cost", "Incre_QALY"  )
    for (i in 1:nrep) {
      first_event_coef[,1] <- psa_pm_firstE_CHD[i,]                 
      first_event_coef[,2] <- psa_pm_firstE_CBVD[i,]   
      first_event_coef[,3] <- psa_pm_fatal_CVD[i,]
      first_event_coef[,4] <- psa_pm_fatal_nonCVD[i,]
      post_event_coef[,1] <- psa_pm_surv_postCHD[i,]                           
      post_event_coef[,2] <- psa_pm_surv_postCBVD[i,]
      psa_nointv <- model(age,SIMD,Diabetes,FH,CPD,SBP,TC,HDL,sex,first_event_coef,post_event_coef, secondary_event_inputs, coef_c1,coef_c2,disc, tm)                                   # can set discount rate to be anything here.
      psa_result[i,1] <- psa_nointv$LY
      psa_result[i,2] <- psa_nointv$Cost
      psa_result[i,3] <- psa_nointv$QALY
      psa_intv <- model(age,SIMD,Diabetes,FH,CPD,SBP,TC=TC_statin,HDL=HDL_statin,sex,first_event_coef,post_event_coef, secondary_event_inputs, coef_c1,coef_c2,disc, tm)   # discounted
      psa_result[i,4] <- psa_intv$LY
      psa_result[i,5] <- psa_intv$Cost
      psa_result[i,6] <- psa_intv$QALY
      psa_result[i,7] <- sum(statin_anucum * rowSums(psa_intv$prob_event_yr))                # STATIN additional cost            
      psa_result[i,8] <- sum(statin_disu_anucum * rowSums(psa_intv$prob_event_yr))          # STATIN disutility due to taking pills
      }
    
  
    ### 4.5 plot the PSA results on CE plane
    psa_result[,9] <- psa_result[,5] + psa_result[,7] - psa_result[,2]
    psa_result[,10] <- psa_result[,6] - psa_result[,8]- psa_result[,3]
    psa_result <- as.data.frame(psa_result)
    
    # set up x y axes limit
    min_c <- min(psa_result[,9])
    max_c <- max(psa_result[,9])
    min_q <- min(psa_result[,10])
    max_q <- max(psa_result[,10])
    
    lim_fn <- function(minpara,maxpara) {
    if ((minpara<0 & maxpara<0) | (abs(minpara)>abs(maxpara))) {
      lim_min <- minpara
      lim_max <- -minpara
    } else  {
      lim_max <- maxpara
      lim_min <- -maxpara
    }
      list = list(min=lim_min, max=lim_max)
      return(list)
    }
    
    ylim <- lim_fn(min_c,max_c)
    xlim <- lim_fn(min_q,max_q)
    
    # plot
    ggplot(psa_result,aes(x=Incre_QALY, y=Incre_cost)) +
      geom_point(colour="violetred4") + 
      ylab("Cost (£)") +
      xlab("QALY") +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      xlim(xlim$min, xlim$max) +
      ylim(ylim$min, ylim$max) +
      geom_abline(aes (intercept = 0, slope = 20000, color= "blue"),linetype="dashed") +
      scale_color_identity(labels=c("£20k"), guide="legend")

    
    ### 4.6 plot CEAC

    prob_nmb <- NA
    for (i in 0:50000) {
        threshold <- i
        nmb <- psa_result[,10]* threshold - psa_result[,9]
        prob_nmb[i+1] <- sum(nmb > 0) / length(nmb)
    }
    
   ceac <- data.frame(x=seq(0, 50000, by=1), y=prob_nmb)
    
    ggplot(ceac,aes(x=x, y=y)) +
      geom_point(colour="violetred4") + 
      ylab("Probability for the intervention to be cost-effective") +
      xlab("ICER threshold (£)") +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      xlim(0, 50000) +
      ylim(0, 1) +
      geom_vline(xintercept = 20000, color="blue", linetype="dashed")
    
    
    
    
    ### I-IV : ---------->> Final results: base case & PSA
    base_nointv$LY                                                                  # No intervention
    quantile(psa_result[,1], c(0.025, 0.975), type=1)
    base_nointv$Cost
    quantile(psa_result[,2], c(0.025, 0.975), type=1)
    base_nointv$QALY
    quantile(psa_result[,3], c(0.025, 0.975), type=1)
    
    base_intv$LY                                                                    # With intervention, cost including statin cost, QALY including statin disutility
    quantile(psa_result[,4], c(0.025, 0.975), type=1)
    base_intv$Cost
    quantile(psa_result[,5] + psa_result[,7], c(0.025, 0.975), type=1)
    base_intv$QALY
    quantile(psa_result[,6] - psa_result[,8], c(0.025, 0.975), type=1)
    
 
