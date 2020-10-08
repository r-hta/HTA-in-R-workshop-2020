// Gompertz survival model

functions {
  // Defines the log hazard
  real gompertz_lpdf(real t, real SCALE, real SHAPE) {
    return log( SCALE * SHAPE) + (SCALE * t);
    }
  
  // Defines the log survival
  real gompertz_lccdf (real t, real SCALE, real SHAPE) {
    return  -SHAPE * (exp(SCALE * t) - 1);
    }
  
  // Defines the sampling distribution
  real gompertz_surv(real t, real SCALE, real SHAPE) {
    return exp(-SHAPE * (exp(SCALE * t) - 1));
    }
}

data {
  int<lower=0> N;
  int<lower=0> N_pred;
  int<lower=0> N_cens;
  real t[N];
  real c[N_cens];
  vector[N_pred] t_new;
  real St_upper;
  real St_lower;
  real pred_time;



  //int n;                  // number of observations
  //vector[n] t;            // observed times
  //vector[n] d;            // censoring indicator (1=observed, 0=censored)
  //real MEAN_LOG_SCALE_PARAMETER;
  //real<lower=0> SE_LOG_SCALE_PARAMETER;
  //vector[n] trt; 
  //real<lower=0> S_low_active;
  //real<lower=0> S_upp_active;
  //real<lower=0> time_ev;
    
}

parameters {
  real LOG_SCALE_PARAMETER;
  real LOG_SHAPE_PARAMETER;
 
}

model {
  LOG_SCALE_PARAMETER ~ normal(0, 5);
  
  LOG_SHAPE_PARAMETER ~ normal( 
      (log
        (-log(St_lower)/(exp(exp(LOG_SCALE_PARAMETER)*pred_time)-1)) +
      log(
        -log(St_upper)/(exp(exp(LOG_SCALE_PARAMETER)*pred_time)-1)))/2,
        
        (log(
          -log(St_lower)/(exp(exp(LOG_SCALE_PARAMETER)*pred_time)-1)) -
        log(
          -log(St_upper)/(exp(exp(LOG_SCALE_PARAMETER)*pred_time)-1)))/3.92
      );

    for(i in 1:N) {
    t[i] ~ gompertz(exp(LOG_SCALE_PARAMETER), exp(LOG_SHAPE_PARAMETER));
  }
  
 for(i in 1:N_cens){
    target += gompertz_lccdf(c[i]| exp(LOG_SCALE_PARAMETER), exp(LOG_SHAPE_PARAMETER));
  }
  
  
//t ~ surv_gompertz( d , exp( LOG_SCALE_PARAMETER ) , exp(LOG_SHAPE_PARAMETER));

}

generated quantities {
  vector[N] log_lik;
  vector[N_pred] surv_pred;
  
  for (i in 1:N) {
    log_lik[i] <- gompertz_lpdf(t[i]| exp(LOG_SCALE_PARAMETER), exp(LOG_SHAPE_PARAMETER));
  }
  
  for(j in 1:N_pred) {
    surv_pred[j] <- gompertz_surv(t_new[j], exp(LOG_SCALE_PARAMETER), exp(LOG_SHAPE_PARAMETER));
  } 
}
