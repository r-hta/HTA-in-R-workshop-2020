bayesian_surv <- function(model, data, path, km.data,data2){
  dist.chosen <- model
  dist.name   <- dist.chosen[3]
  par <- c(dist.chosen[1],dist.chosen[2]) 
  par1 <- dist.chosen[1]
  par2 <- dist.chosen[2]
  
  St_upper <- data[["St_upper"]]
  St_lower <- data[["St_lower"]]
  pred_time <- data[["pred_time"]]
  
  
  #Fit model
  model_fit <- stan(paste0(path, dist.chosen[4]), 
                    data = data, iter = 5000, chains = 2)
  
  # Check to see mixing of the chains
  
  mixingplot <- rstan::traceplot(model_fit, par = par)
  
  # Check to see autocorrelation in the simulations
  autocorrelationplot <- bayesplot::mcmc_acf(as.matrix(model_fit), pars = par)
  
  #Plot the posterior distributions
  posteriorplot <- bayesplot::mcmc_areas(as.matrix(model_fit), pars= par, prob = 0.95)
  
  #Extract out the draws
  survival_model_draws <- tidybayes::tidy_draws(model_fit)
  
  Surv_model <- function(par1, par2, time_vector, dist.name){
    Surv_vector <- rep(NA, length(time_vector))
    
    if(dist.name == "Weibull"){
      sapply(time_vector, function(x)exp(-(exp(par2)*x)^par1) ) #par1 is "a", par2 is "logb"
    }else if(dist.name == "Exponential"){
      sapply(time_vector, function(x)exp(-(exp(par2)*x)^par1) ) #   #par1 is "1", par2 is "b"
    }else if(dist.name == "LogNormal"){
      sapply(time_vector, function(x)plnorm(x,meanlog= par1, sdlog = par2,lower.tail = FALSE, log.p = FALSE)) #   #par1 is "a", par2 is "b"
    }else if(dist.name == "LogLogistic"){
      sapply(time_vector, function(x)(1/(1+exp((log(x)-par1)/par2)))) #   #par1 is "a", par2 is "b"
    }else{ #Gompertz
      sapply(time_vector, function(x)exp(-exp(par1)*(exp(exp(par2)*x)-1))) #par1 is "log_n", par2 is "log_b"
    }
  }
  
  
  extractSurv <- function(df,time_vector, par1, par2,dist.name){
    
    df_Surv <- matrix(NA, nrow = length(time_vector), ncol = nrow(df))
    
    for(i in 1:nrow(df)){
      
      df_Surv[,i] <- Surv_model(as.numeric(df[i,par1]), as.numeric(df[i,par2]), time_vector, dist.name)
    }
    colnames(df_Surv) <- 1:nrow(df)
    df_Surv <- data.frame(cbind(time_vector,df_Surv))
    df_Surv <- gather(df_Surv, key = iter, value = Survival, -time_vector)
    return(df_Surv)
  }
  
  Survivaldf <- extractSurv(df = survival_model_draws[sample(nrow(survival_model_draws), 250), ], time_vector = seq(from = 0, to = 500, by = 10),
                            par1, par2, dist.name) #We don't want to plot all the simulations so we take a random
  
  
  #Summary Stats
  model_summary <- data.frame(summary(model_fit, probs = c(0.025, 0.975))$summary)
  model_summary <- model_summary[grepl("^surv_pred", rownames(model_summary)), c("mean","X2.5.","X97.5.") ]
  model_summary$time <- data[["t_new"]]
  
  #https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html
  
  
  #Plot the data
  survplot <- ggplot(data = Survivaldf,
                     mapping = aes(x = time_vector, y = Survival, group = iter)) +
    geom_line(size = 0.1, alpha = 0.025, colour = "red")+
    scale_y_continuous(breaks=seq(0, 1, by = 0.05), limits = c(0,1))+
    scale_x_continuous(breaks=seq(0, 500, by = 50))+
    geom_line(data = model_summary,
              mapping = aes(x = time, y= X2.5.), colour = "blue",
              linetype = "dotted",inherit.aes = FALSE)+
    geom_line(data = model_summary,
              mapping = aes(x = time, y= X97.5.), colour = "blue",
              linetype = "dotted",inherit.aes = FALSE)+
    geom_line(data = model_summary ,
              mapping = aes(x = time, y= mean),colour = "blue",inherit.aes = FALSE)+
    geom_step(data = km.data, aes(x = Time, y = Survival), inherit.aes = F )+
    geom_step(data = km.data, aes(x = Time, y = upper),linetype = "dashed", inherit.aes = F )+
    geom_step(data = km.data, aes(x = Time, y = lower), linetype = "dashed", inherit.aes = F )+
    geom_path(data = data2, aes(x =x, y =y), colour = "salmon", lwd=1.1, inherit.aes = F)+
    ggtitle(label = paste0("Bayesian ",dist.name," Survial Model"),
            subtitle = paste0("Prior belief of ",St_lower*100," - ",St_upper*100,"% ",
                              "survival at ",pred_time," weeks"))+
    labs(x = "Time in weeks") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, face = "italic"))
  #print(survplot)
  #print(posteriorplot)
  return(list(model_fit, survplot))
}
