  
  # This Scottish Cardiovascular disease policy model is published in the following two papers:
  # Lewsey JD, et al. A cardiovascular disease policy model that predicts life expectancy taking into account socioeconomic deprivation. Heart 2015;101:201-208. doi:10.1136/heartjnl-2014-305637.
  # Lawson KD, et al. A cardiovascular disease policy model: part 2-preparing for economic evaluation and to assess health inequalities. Open Heart 2016;3:e000140. doi:10.1136/openhrt-2014-000140.
  # Here is the R shiny app for it, which is developed based on the R version of the model.



  library(shiny)
  library(shinydashboard)
  library(shinydashboardPlus)
  library(shinyWidgets)
  library(shinydashboard)
  
  library(dplyr)
  library(tidyr)
  library(expss)                                                                            
  library(MASS)
  library(ggplot2)
  
  
  source("CVD_fns.R",local = TRUE)
  
  
  ###############################################################################
  ################################### UI function ###############################
  ###############################################################################
  
  
  ui <- fluidPage(
      list(tags$head(
        tags$style(HTML("
        .navbar .navbar-nav {float: right; 
                              } ")))),  
      navbarPage(id="cvd",theme="custom-navbar.css",
                       "Health economics modelling in R",
                   header = tagList(
                     useShinydashboard()
                   ),
          tabPanel(id="Model", "Model",
                        tags$head(tags$script('
                          var dimension = [0, 0];
                          $(document).on("shiny:connected", function(e) {
                              dimension[0] = window.innerWidth;
                              dimension[1] = window.innerHeight;
                              Shiny.onInputChange("dimension", dimension);
                          });
                          $(window).resize(function(e) {
                              dimension[0] = window.innerWidth;
                              dimension[1] = window.innerHeight;
                              Shiny.onInputChange("dimension", dimension);
                          });
                      ')),
               fluidRow(
                   column(1),
                   column(10,
                       h2(tags$strong("Scottish cardiovascular disease policy model")),
                       br(),
                       p(tags$strong("A cardiovascular disease (CVD) policy model which can be used to model remaining life expectancy, quality-adjusted life-years (QALYs) and health care costs, including a measure of socioeconomic deprivation as an independent risk factor for CVD.")),
                       br(),
                       p(tags$a(href="https://heart.bmj.com/content/101/3/201", icon("mouse-pointer"), 
                                tags$strong("Lewsey JD, et al. A cardiovascular disease policy model that predicts life expectancy taking into account socioeconomic deprivation. Heart 2015;101:201-208. doi:10.1136/heartjnl-2014-305637.", 
                                style="color:#2196c4; background-color: #FFFFFF"),target="_blank")),
                       p(tags$a(href="https://openheart.bmj.com/content/3/1/e000140", icon("mouse-pointer"), 
                                tags$strong("Lawson KD, et al. A cardiovascular disease policy model: part 2-preparing for economic evaluation and to assess health inequalities. Open Heart 2016;3:e000140. doi:10.1136/openhrt-2014-000140.", 
                                            style="color:#2196c4; background-color: #FFFFFF"),target="_blank")),
                       p(tags$strong("Code for this model can be downloaded from Github repository by clicking"),
                         tags$a(href="https://github.com/yiqiaoxin/CVDmodel/", 
                                tags$strong("here",style="color:#2196c4; background-color: #FFFFFF"),target="_blank"), tags$strong(".")),
                       br(),
                       sidebarLayout(
                         sidebarPanel(
                                               h4(tags$strong("Please enter the patient characteristics below:")),
                                               br(),
                                               numericInput("age", "Age (Please enter a value between 0 - 120)", value = 60),
                                               radioButtons("sex", "Sex", c("Male" = 1,"Female" = 2), inline = TRUE),
                                               numericInput("simd", label = p("Scottish Index of Multiple Deprivation", tags$a(href="https://simd.scot/", 
                                                                                   tags$strong("(https://simd.scot)"),target="_blank"),
                                                                                   "(please enter a value between 0 (least deprived) and 100 (most deprived))")
                                                                                   , value = 60.8),
                                               radioButtons("dia", "Diabetes", c("With diabetes" = 1,"No diabetes" = 0),  selected = "0", inline = TRUE),
                                               radioButtons("fh", "Family history", c("With family history of CVD" = 1,"No family history of CVD" = 0), selected = "0", inline = TRUE),
                                               numericInput("cpd", "Cigarette per day", value = 20),                          
                                               numericInput("sbp", "Systolic blood pressure", value = 160),       
                                               numericInput("tc", "Total cholestorol level (mmol/L) (below 5.2 mmol/L is desirable)", value = 7),
                                               numericInput("hdl", "High density lipoproteins (mmol/L) (above 1.5 mmol/L is desirable)", value = 1),
                                               br(),
                                               h4(tags$strong("Treatment information:")),
                                               numericInput("statin_hdle", "Statin improves high density lipoproteins levels by (%):", value = 4),
                                               numericInput("statin_ldle", "Statin lowers low density lipoproteins levels by (%)", value = 26),
                                               numericInput("statin_cost", "Annual cost of statin (GBP)", value = 13),
                                               numericInput("statin_disu", "Disutility due to taking statin pills", value = 0.001),
                                               br(),
                                               h4(tags$strong("Model parameters:")),
                                               numericInput("ds", "Please enter a discount rate:", value = 0.035),
                                               numericInput("tm", "Please enter a time horizon", value = 100)
                                        ),
                         mainPanel(
                           tabsetPanel(
                             tabPanel("1. Model intro",
                                             br(),
                                             br(),
                                             p("This is a state transition model developed using the Scottish Heart Health Extended Cohort (SHHEC) linked to Scottish morbidity and 
                                               death records. Individuals start in a CVD-free state and can transit to three CVD event states plus a non-CVD death state. 
                                               Individuals who have a non-fatal first event are then followed up until death. Taking a competing risk approach, 
                                               the cause-specific hazards of a first event are modelled using parametric survival analysis. Survival following a first non-fatal 
                                               event is also modelled parametrically."),
                                             br(),
                                             p("To generate quality-adjusted life expectancy (QALE), the Scottish Health Survey was used to estimate background health utilities and the 
                                               impact of CVD events (utility decrements). The SF-6D algorithm generated utilities and decrements were modelled using ordinary least squares (OLS). 
                                               To generate lifetime hospital costs, the Scottish Heart Health Extended Cohort (SHHEC) was linked to the Scottish morbidity and death records (SMR) 
                                               to cost each continuous inpatient stay (CIS). OLS and restricted cubic splines estimated annual costs before and after each of the first four events. 
                                               A Kaplan-Meier sample average (KMSA) estimator was then used to weight expected health-related quality of life and costs by the probability of survival."),
                                             br(),
                                      column(1),
                                      column(10,
                                             img(src='CVD.jpg', width=600, height=400)),
                                      br(),
                                      br(),
                                      hr())
                                      ,
                             tabPanel("2. Base case results",
                                      br(),
                                      p("The base case analysis results are shown in the colored tiles below."),

                                      htmlOutput("discshow"),
                                      br(),
                                      br(),
                                      h4(tags$strong("Usual care (without statin)")),
                                      fluidRow(box(width = 12,
                                                   tags$head(tags$style(HTML(".small-box {height: 200px}"))),
                                                   valueBoxOutput("ly"),
                                                   valueBoxOutput("qaly"),
                                                   valueBoxOutput("cost")
                                                   
                                      )),
                                      h4(tags$strong("Intervention (with statin)")),
                                      fluidRow(
                                      box(width = 12,
                                          tags$head(tags$style(HTML(".small-box {height: 200px}"))),
                                          valueBoxOutput("intvly"),
                                          valueBoxOutput("intvqaly"),
                                          valueBoxOutput("intvcost")
                                          
                                      )),
                                      h4(tags$strong("Incremental outcomes"), style = "cursive;
                                      font-weight: 500; line-height: 1.1; 
                                      color: #FFFFFF;background-color: #00022e"),
                                      fluidRow(
                                        box(width=12,
                                            tags$head(tags$style(HTML(".small-box {height: 200px}"))),
                                            valueBoxOutput("increly"),
                                            valueBoxOutput("increqaly"),
                                            valueBoxOutput("increcost")
                                      )
                                      )),
                             tabPanel("3. PSA",
                                      br(),
                                      numericInput("ite", "Please enter number of iterations for the PSA:", value=100),
                                      actionButton("runpsa", "Please click here to run PSA"),
                                      br(),
                                      br(),
                                      h4(tags$strong("Incremental cost-effectiveness plane")),
                                      plotOutput("plotplane"),
                                      br(),
                                      h4(tags$strong("Cost-effectiveness acceptability curve")),
                                      plotOutput("plotceac")
                                      )
                             )))),
                   column(1)
              ),
              br(),
              hr(),
              p(em("Developed by"),br("Yiqiao Xin (yiqiao.xin@glasgow.ac.uk), Ewan Gray, Jose Antonio Robles-Zurita, Houra Haghpanahan, Robert Heggie, Ciaran Kohli-Lynch, Andrew Briggs, David McAllister,Kenny Lawson,  Jim Lewsey (Jim.Lewsey@glasgow.ac.uk)
                "),style="text-align:center; font-family: times")
              ),
          tabPanel(id="Intro", "Project introduction",
               column(1),
               column(10,
                      h2(tags$strong("Project Introduction")),
                      br(),
                      h4("Background",style = "cursive;
                        font-weight: 500; line-height: 1.1;
                        color: #FFFFFF;background-color: #2196c4"),
                      br(),
                      p("Many previously developed models are often adapted to be used in multiple projects rather than 
                        being superseded by de-novo models developed from scratch. Given the easy adaptability of models 
                        in R as well as the known benefits such as running efficiency and transparency, there is a growing 
                        interest to convert the existing Excel based models into R script for model re-use. The Scottish 
                        cardiovascular disease (CVD) policy model (", 
                        tags$a(href="https://heart.bmj.com/content/101/3/201",tags$strong("Lewsey JD 2015",style="color:#2196c4; background-color: #FFFFFF"),target="_blank"),
                      ",",
                        tags$a(href="https://openheart.bmj.com/content/3/1/e000140",tags$strong("Lawson KD 2016",style="color:#2196c4; background-color: #FFFFFF"),target="_blank"),
                        ") was initially developed in 
                        Excel and has now been converted into R. This Shiny app uses the R version of the model."),
                      br(),
                      h4("Contributors",style = "cursive;
                        font-weight: 500; line-height: 1.1;
                        color: #FFFFFF;background-color: #2196c4"),
                      br(),
                      p("Yiqiao Xin",tags$sup("1"), ", Ewan Gray",tags$sup("2"),", Jose Antonio Robles-Zurita", tags$sup("1"), ", Houra Haghpanahan", tags$sup("1"),", Robert Heggie", tags$sup("1"),", Ciaran Kohli-Lynch", tags$sup("1"),
                        ", Andrew Briggs", tags$sup("3"), ", David McAllister", tags$sup("4"),", Kenny Lawson", tags$sup("5"), ", and Jim Lewsey", tags$sup("1"), "."),
                      p(tags$strong("Affiliations")),
                      p("1 Health Technology Assessment and Health Economics (HEHTA), Institute of Health and Wellbeing, University of Glasgow, UK. "),
                      p("2 Edinburgh Clinical Trials Unit (ECTU), University of Edinburgh."),
                      p("3 Department of Health Services Research and Policy, London School of Hygiene and Tropical Medicine, London, UK."),
                      p("4 Public Health, Institute of Health and Wellbeing, University of Glasgow, Glasgow, UK."),
                      p("5 Translational Health Research Institute. Western Sydney University."),
                      br(),
                      h4("Feedback",style = "cursive;
                        font-weight: 500; line-height: 1.1;
                        color: #FFFFFF;background-color: #2196c4"),
                      br(),
                      p("We hope similar tools can be made for other health economics models in the future.
                         We would be very grateful to anyone who sends constructive feedback and suggestions. Please send any feedback to",
                      tags$a(href="mailto:yiqiao.xin@glasgow.ac.uk", "Dr Yiqiao Xin yiqiao.xin@glasgow.ac.uk",target="_blank"), "or",
                      tags$a(href="mailto:jim.lewsey@glasgow.ac.uk", "Prof Jim Lewsey jim.lewsey@glasgow.ac.uk",target="_blank"),"."),
                      br(),
                      br(),
                      br(),
                      p("THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
                         NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
                         IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
                         WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
                         OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."),
                      br()
                      ),
              column(1)
                  )
          )
)
  
  
  
  
  ###############################################################################
  ################################ Server function ##############################
  ###############################################################################

  server <- function(input, output, session) {
    
    output$discshow <- renderText({
      paste("Please note: you have set discount rate as", "<font color=\"#2196c4\"><b>" , input$ds,"</b></font>", ".")
    })
    
    ### base case result
    
    baseresult <- function(){
                    age <- input$age
                    SIMD <- input$simd                                                                          
                    Diabetes <- as.integer(input$dia)                                                              # diabetes: 0 - no diabetes; 1 - with diabetes
                    FH <- as.integer(input$fh)                                                                     # family history: 0 - no history ; 1 - with history
                    CPD <- input$cpd                                                                               # No. of cigarette per day (continuous)
                    SBP <- input$sbp                                                                               # Systolic blood pressure
                    TC <- input$tc                                                                                 # Total cholestorol 
                    HDL <- input$hdl
                    sex <- as.integer(input$sex)                                                                   # male. sex <-2 female
                    tm <- input$tm                                                                                 # As consistent with the Excel model - time_horizon set as 115 cycles. Alternatively, if setting tm to be bounded by age 100: tm <- 100 - age
                    disc <-input$ds  
                    
                    statin_ldleffect <- input$statin_ldle / 100                                                    # treatment effects on ldl and hdl
                    statin_hdleffect <- input$statin_hdle / 100+1
                    statin_cost <- input$statin_cost
                    statin_disu <- input$statin_disu
                    
                    source("CVD_parameters.R",local = TRUE)
                    
                    if(sex == 1){                                                                                  # set a switch to choose to use whether male or female parameters.
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
  
                    base_nointv <- model(age,SIMD,Diabetes,FH,CPD,SBP,TC,HDL,sex,first_event_coef,post_event_coef, secondary_event_inputs, coef_c1,coef_c2,disc, tm)   

                    
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
                    statin_anucum_sur <- statin_anucum * rowSums(base_intv$prob_event_yr)                         # take statin every year until the first event
                    sum(statin_anucum_sur)                                                                        # Result
                    
                    # STATIN disutility
                    statin_disu_anu <- statin_disu / (1+disc)^cycle
                    statin_disu_anucum <- cumsum(statin_disu_anu)
                    statin_disu_anucum_sur <- statin_disu_anucum * rowSums(base_intv$prob_event_yr)  
                    sum(statin_disu_anucum_sur)                                                                   # Result
                    
                    # Intervention overall
                    sum(statin_anucum_sur) + base_intv$Cost
                    base_intv$QALY - sum(statin_disu_anucum_sur)
                    
                    # Incremental effect
                    incre_cost <- sum(statin_anucum_sur) + base_intv$Cost - base_nointv$Cost
                    incre_cost
                    incre_QALY <- base_intv$QALY - sum(statin_disu_anucum_sur) - base_nointv$QALY
                    incre_QALY
                    
                    list(base_nointvLY = base_nointv$LY, base_nointv_Cost = base_nointv$Cost, base_nointvQALY =base_nointv$QALY, 
                         base_intvLY = base_intv$LY, base_intv_Cost = base_intv$Cost, base_intvQALY =base_intv$QALY,                    
                         statincost = sum(statin_anucum_sur), statin_disu = sum(statin_disu_anucum_sur),
                         base_intv_allCost = sum(statin_anucum_sur) + base_intv$Cost, base_intv_allQALY = base_intv$QALY - sum(statin_disu_anucum_sur),
                         incre_cost = incre_cost, incre_QALY= incre_QALY
                         )
                    
                    }

    output$ly <- renderValueBox({
      valueBox(round(baseresult()$base_nointvLY- input$age,2),"Remaining life years", 
               icon = icon("tree"), color = "green")
    })
    output$qaly <- renderValueBox({
      valueBox(round(baseresult()$base_nointvQALY -input$age,2),"Remaining quality adjusted life years", 
               icon = icon("laugh-beam"), color = "yellow")
    })
    output$cost <- renderValueBox({
      valueBox(format(baseresult()$base_nointv_Cost, big.mark=",", big.interval=3L,
                      digits=0, scientific=F), "Health care costs (GBP)", 
               icon = icon("coins"), color = "light-blue")
    })
    
    output$intvly <- renderValueBox({
      valueBox(round(baseresult()$base_intvLY- input$age,2),"Remaining life years", 
               icon = icon("tree"), color = "green")
    })
    output$intvqaly <- renderValueBox({
      valueBox(round(baseresult()$base_intv_allQALY -input$age,2),"Remaining quality adjusted life years", 
               icon = icon("laugh-beam"), color = "yellow")
    })
    output$intvcost <- renderValueBox({
      valueBox(format(baseresult()$base_intv_allCost, big.mark=",", big.interval=3L,
                      digits=0, scientific=F), "Health care costs (GBP)", 
               icon = icon("coins"), color = "light-blue")
    })
    
    output$increly <- renderValueBox({
      valueBox(round(baseresult()$base_intvLY-baseresult()$base_nointvLY,2),"Incremental life years", 
               icon = icon("tree"), color = "navy", width=6)
    })
    
    output$increqaly <- renderValueBox({
      valueBox(round(baseresult()$incre_QALY,2),"Incremental quality adjusted life years", 
               icon = icon("laugh-beam"), color = "navy", width=6)
    })
    output$increcost <- renderValueBox({
      valueBox(format(baseresult()$incre_cost, big.mark=",", big.interval=3L,
                      digits=0, scientific=F), "Incremental health care costs (GBP)", 
               icon = icon("coins"), color = "navy", width=6)
      
    }) 
    
    
    
    ### PSA

    psa <-  eventReactive(input$runpsa, {

      age <- input$age
      SIMD <- input$simd
      Diabetes <- as.integer(input$dia)
      FH <- as.integer(input$fh)
      CPD <- input$cpd
      SBP <- input$sbp
      TC <- input$tc
      HDL <- input$hdl
      sex <- as.integer(input$sex)
      tm <- input$tm
      disc <-input$ds
      statin_cost <- input$statin_cost
      statin_disu <- input$statin_disu
      statin_ldleffect <- input$statin_ldle / 100                                                                 # treatment effects on ldl and hdl
      statin_hdleffect <- input$statin_hdle / 100+1
      TC_statin <-TC-(TC-HDL)*statin_ldleffect
      HDL_statin <- HDL*statin_hdleffect

      source("CVD_parameters.R",local = TRUE)
      
      # Statin cost
      cycle <- seq(from = 1, to = tm,by = 1) 
      statin_anu <- statin_cost / (1+disc)^cycle
      statin_anucum <- cumsum(statin_anu)
                                                                  
      
      # STATIN disutility
      statin_disu_anu <- statin_disu / (1+disc)^cycle
      statin_disu_anucum <- cumsum(statin_disu_anu)
 

      if(sex == 1){
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

      ### Set the No. of simulations
      nrep <- input$ite

      ### Set the variance/covariance matrix according to sex
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


      ## Generate draws

      progress <- shiny::Progress$new()   # Adding progress bars
      on.exit(progress$close())
      progress$set(message="Conducting PSA. Please wait", value=0)

      mu_firstE_CHD <- first_event_coef[,1]                                                                       # get the mean coefficient from the base case
      set.seed(12345)
      psa_pm_firstE_CHD <- MASS::mvrnorm(n=nrep, mu = mu_firstE_CHD, Sigma = Varcov_firstE_CHD)

      progress$inc(0.1, detail="Completed 1/6 draws")

      mu_firstE_CBVD <- first_event_coef[,2]
      set.seed(12345)
      psa_pm_firstE_CBVD <- MASS::mvrnorm(n=nrep, mu = mu_firstE_CBVD, Sigma = Varcov_firstE_CBVD)

      progress$inc(0.1, detail="Completed 2/6 draws")

      mu_fatal_CVD <- first_event_coef[,3]
      set.seed(12345)
      psa_pm_fatal_CVD <- MASS::mvrnorm(n=nrep, mu = mu_fatal_CVD, Sigma = Varcov_fatal_CVD)

      progress$inc(0.1, detail="Completed 3/6 draws")

      mu_fatal_nonCVD <- first_event_coef[,4]
      set.seed(12345)
      psa_pm_fatal_nonCVD <- MASS::mvrnorm(n=nrep, mu = mu_fatal_nonCVD, Sigma = Varcov_fatal_nonCVD)

      progress$inc(0.1, detail="Completed 4/6 draws")

      mu_surv_postCHD <- post_event_coef[,1]
      set.seed(12345)
      psa_pm_surv_postCHD <- MASS::mvrnorm(n=nrep, mu = mu_surv_postCHD, Sigma = Varcov_surv_postCHD)

      progress$inc(0.1, detail="Completed 5/6 draws")

      mu_surv_postCBVD <- post_event_coef[,2]
      set.seed(12345)
      psa_pm_surv_postCBVD <- MASS::mvrnorm(n=nrep, mu = mu_surv_postCBVD, Sigma = Varcov_surv_postCBVD)

      progress$inc(0.1, detail="Completed 6/6 draws. Populating all draws into the model...")

      ### loop the model through the draws

      psa_result <- matrix(nrow = nrep, ncol = 10)
      colnames(psa_result) <- c("LY__nointv", "Cost__nointv", "QALY__nointv","LY_intv", "Cost_intv", "QALY_intv",
                                "Cost_intv_add_statin", "disQALY_intv_statin", "Incre_cost", "Incre_QALY"  )
      for (i in 1:nrep) {

        if (i==nrep/4) {
          progress$inc(0.1, detail="Completed 6/6 draws.Populating all draws into the model: 25% completed now ...")
        }

        if (i==nrep/2) {
          progress$inc(0.1, detail="Completed 6/6 draws.Populating all draws into the model: 50% completed now ...")
        }

        if (i==nrep/4*3) {
          progress$inc(0.1, detail="Completed 6/6 draws.Populating all draws into the model: 75% completed now ...")
        }

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
        psa_result[i,7] <- sum(statin_anucum * rowSums(psa_intv$prob_event_yr))                                  # STATIN additional cost
        psa_result[i,8] <- sum(statin_disu_anucum * rowSums(psa_intv$prob_event_yr))                             # STATIN disutility due to taking pills
        
        }
      progress$inc(0.1, detail="Completed. Rendering results now.")
      psa_result[,9] <- psa_result[,5] + psa_result[,7] - psa_result[,2]
      psa_result[,10] <- psa_result[,6] - psa_result[,8]- psa_result[,3]
      psa_result <- as.data.frame(psa_result)
      
    })

    #CE plane
    output$plotplane <- renderPlot({
            psa_result <- psa()
            min_c <- min(psa_result[,9])
            max_c <- max(psa_result[,9])
            min_q <- min(psa_result[,10])
            max_q <- max(psa_result[,10])
      
            ylim <- lim_fn(min_c,max_c)
            xlim <- lim_fn(min_q,max_q)
            
            ceplane_fn(psa_result,xlim$min, xlim$max,ylim$min, ylim$max)
    })
    
    # CEAC
    output$plotceac <- renderPlot({
      psa_result <- psa()
      ceac_fn(psa_result)
    })
      
    
  }
  shinyApp(ui = ui, server = server)