// Gompertz survival model

functions {
  // Defines the log hazard
  vector log_h (vector t, vector SCALE, vector SHAPE) {
    vector[num_elements(t)] log_h;
        for (i in 1:num_elements(t)) {
          log_h[i] = log( SCALE[i] * SHAPE[i] ) + (SCALE[i] * t[i]);
        }
    return log_h;
  }
  
  // Defines the log survival
  vector log_S (vector t, vector SCALE, vector SHAPE) {
    vector[num_elements(t)] log_S;
    for (i in 1:num_elements(t)) {
      log_S[i] = -SHAPE[i] * (exp(SCALE[i] * t[i]) - 1);
    }
    return log_S;
  }
  
  // Defines the sampling distribution
  real surv_gompertz_lpdf (vector t, vector d, vector SCALE, vector SHAPE) {
    vector[num_elements(t)] log_lik;
    real prob;
    log_lik = d .* log_h(t,SCALE,SHAPE) + log_S(t,SCALE,SHAPE);
    prob = sum(log_lik);
    return prob;
  }
}

data {
  int n;                  // number of observations
  vector[n] t;            // observed times
  vector[n] d;            // censoring indicator (1=observed, 0=censored)
  real MEAN_LOG_SCALE_PARAMETER;
  real<lower=0> SE_LOG_SCALE_PARAMETER;
  real MEAN_LOG_SCALE_PARAMETER1;
  real<lower=0> SE_LOG_SCALE_PARAMETER1;
  real MEAN_LOG_SHAPE_PARAMETER1;
  real<lower=0> SE_LOG_SHAPE_PARAMETER1;
  vector[n] trt; 
  real<lower=0> S_low_plac;
  real<lower=0> S_upp_plac;
  real<lower=0> S_low_active;
  real<lower=0> S_upp_active;
  real<lower=0> time_ev;
    
}

parameters {
  real LOG_SCALE_PARAMETER;
  real LOG_SCALE_PARAMETER1;
  real LOG_SHAPE_PARAMETER; 
  real LOG_SHAPE_PARAMETER1; 
  
}

transformed parameters {
  vector[n] linpred;
  vector[n] mu;
  vector[n] linpred2;
  vector[n] scale_vector;
 
  linpred = (1-trt)*LOG_SHAPE_PARAMETER + trt*LOG_SHAPE_PARAMETER1;
  for (i in 1:n) {
    mu[i] = exp(linpred[i]);
  }

  linpred2 = (1-trt)*LOG_SCALE_PARAMETER + trt*LOG_SCALE_PARAMETER1;
  for (i in 1:n) {
    scale_vector[i] = exp(linpred2[i]);
  }

}

model {
  LOG_SCALE_PARAMETER ~ normal( MEAN_LOG_SCALE_PARAMETER , SE_LOG_SCALE_PARAMETER );
  LOG_SCALE_PARAMETER1 ~ normal( MEAN_LOG_SCALE_PARAMETER1 , SE_LOG_SCALE_PARAMETER1 );

  LOG_SHAPE_PARAMETER ~ normal( 
      (log
        (-log(S_low_plac)/(exp(exp(LOG_SCALE_PARAMETER)*time_ev)-1)) +
      log(
        -log(S_upp_plac)/(exp(exp(LOG_SCALE_PARAMETER)*time_ev)-1)))/2,
        
        (log(
          -log(S_low_plac)/(exp(exp(LOG_SCALE_PARAMETER)*time_ev)-1)) -
        log(
          -log(S_upp_plac)/(exp(exp(LOG_SCALE_PARAMETER)*time_ev)-1)))/3.92
      );


  LOG_SHAPE_PARAMETER1 ~ normal( 
      (log
        (-log(S_low_active)/(exp(exp(LOG_SCALE_PARAMETER1)*time_ev)-1)) +
      log(
        -log(S_upp_active)/(exp(exp(LOG_SCALE_PARAMETER1)*time_ev)-1)))/2,
        
        (log(
          -log(S_low_active)/(exp(exp(LOG_SCALE_PARAMETER1)*time_ev)-1)) -
        log(
          -log(S_upp_active)/(exp(exp(LOG_SCALE_PARAMETER1)*time_ev)-1)))/3.92
      );


  t ~ surv_gompertz( d , exp( LOG_SCALE_PARAMETER ) , exp(LOG_SHAPE_PARAMETER));
//  t ~ surv_gompertz( d , scale_vector , mu);
}

generated quantities {
  real perc; real perc1;               // Checking whether really through right percentage
  perc =exp(-exp(LOG_SHAPE_PARAMETER) * (exp(exp(LOG_SCALE_PARAMETER) * time_ev) - 1));
    perc1 =exp(-exp(LOG_SHAPE_PARAMETER1) * (exp(exp(LOG_SCALE_PARAMETER1) * time_ev) - 1));
}

