# Load packages ----
list.of.packages <- c("ggplot2", "rstan", "survival", "flexsurv", "bayesplot", "survminer", "asaur", 
                      "tidybayes", "dplyr", "tidyr", "gridExtra","ggfortify")

#Check to see if these are installed and install if not
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#load the packages
lapply(list.of.packages, require, character.only = TRUE)

# setwd('C:/Users/phili/Desktop/')
# Source helpers ----
source("helpers.R")

#Apply an rstan setting
rstan_options(auto_write = TRUE)
#Path for Stan files
# path <- paste0('C:/Users/phili/Desktop/', 'Homebrew Distributions/')
path <- 'Homebrew Distributions/'


# Load Dataframe

df <- gastricXelox
result.km <- survfit(Surv(timeWeeks, delta) ~ 1, data = df, conf.type="log-log")
km.data <- data.frame(cbind(result.km[[c("time")]],result.km[[c("surv")]], result.km[[c("upper")]],result.km[[c("lower")]]))
colnames(km.data) <- c("Time", "Survival", "upper", "lower")


#Organize Data for Bayesian analysis
y_cens <- df[which(df$delta == 0), ]$timeWeeks 
N_cens <- length(y_cens)
y_obs  <- df[which(df$delta == 1), ]$timeWeeks
N_obs  <- length(y_obs)

dist.exponential <- c("a","log_b", "Exponential", "Exponential.stan")
dist.lognormal   <- c("a","b", "LogNormal","LogNormal.stan")
dist.loglogistic <- c("a","b", "LogLogistic","LogLogistic.stan")
dist.weibull     <- c("a", "log_b", "Weibull","Weibull.stan")
dist.gompertz    <- c("LOG_SHAPE_PARAMETER", "LOG_SCALE_PARAMETER", "Gompertz" ,"Gompertz.stan")


# User interface ----
ui <- fluidPage(
  titlePanel("Survival Modelling using Expert Opinion"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Expert Opinion on Survival Probabilities"),
      sliderInput("timepoint", label = "Timepoint",
                  min = 0, max = 500, value = 200),
      sliderInput("probs", label = "Upper and Lower Survival Probabilities",
                  min = 0, max = 1, value = c(.25, .5)),
      plotOutput("plot_dist")),
    
    mainPanel(
      h3("Kaplan Meier Survival Plot"),
      plotOutput("plot"),
      selectInput("var", label = "Choose a survival model", 
                  choices = list("Exponential",
                                 "Weibull",
                                 "Gompertz",
                                 "Log-Logistic",
                                 "Log-normal"), 
                  selected = "Weibull"),
      actionButton("do", "Run Analysis"),
      br(),
      h3("Output of the Bayesian Analysis"),
      plotOutput("plot_result")
    
      )
  )
)

# Server logic
server <- function(input, output) {
  
  output$plot <- renderPlot({
    
    mean.surv <- ((input$probs[2] + input$probs[1])/2)
    sd.surv <- ((input$probs[2] - input$probs[1])/(2*1.96))
    
    df2 <- data.frame(y = seq(0, 1, 0.01),
                      x = input$timepoint -(100*sd.surv)*dnorm(seq(0, 1, 0.01),
                                        mean = mean.surv, sd = sd.surv,
                                        log = FALSE))
    ggplot() +
      geom_step(data = km.data, aes(x = Time, y = Survival), inherit.aes = F )+
      geom_step(data = km.data, aes(x = Time, y = upper),linetype = "dashed", inherit.aes = F )+
      geom_step(data = km.data, aes(x = Time, y = lower), linetype = "dashed", inherit.aes = F )+
      geom_path(data = df2, aes(x =x, y =y), colour = "salmon", lwd=1.1)+
      scale_x_continuous(limits = c(0, 500))+
      labs(x = "Time", y = "Survival" )
    
    })
  
  output$plot_dist <- renderPlot({
    mean.surv <- (input$probs[2] + input$probs[1])/2
    sd.surv <- (input$probs[2] - input$probs[1])/(2*1.96)
    ggdistribution(dnorm, seq(0, 1, 0.01), mean = mean.surv, sd = sd.surv)
    
  })
  
  observeEvent(input$do, {
    t_new <- seq(from = 0, to = 500, by = 10)
    N_pred <- length(t_new)

    data <- list(t = y_obs, c = y_cens,
                 N_cens = N_cens, N = N_obs,
                 N_pred = N_pred, t_new = t_new, St_upper = input$probs[2],
                 St_lower = input$probs[1], pred_time = input$timepoint)
    
    #Duplication
    mean.surv <- (input$probs[2] + input$probs[1])/2
    sd.surv <- (input$probs[2] - input$probs[1])/(2*1.96)
    
    df2 <- data.frame(y = seq(0, 1, 0.01),
                      x = input$timepoint -(100*sd.surv)*dnorm(seq(0, 1, 0.01),
                                                               mean = mean.surv, sd = sd.surv,
                                                               log = FALSE))
    
    distribution <- switch(input$var, 
                   "Exponential" = dist.exponential,
                   "Weibull" = dist.weibull,
                   "Gompertz" = dist.gompertz,
                   "Log-Logistic" = dist.loglogistic,
                   "Log-normal" = dist.lognormal)
    
    fit <-bayesian_surv(model = distribution, data = data, path = path, km.data =km.data, data2= df2 )
    output$plot_result <- renderPlot(fit[[2]])
    })
  

}

# Run the app
shinyApp(ui, server)
