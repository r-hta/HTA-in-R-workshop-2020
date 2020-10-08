functions {

  real weib_lpdf(real t, real a, real b){
   return log(b) + log(a) + (a-1)*log(b*t) -(b*t)^a;
     }

  real weib_lccdf(real t, real a, real b) {
    return  -(b*t)^a ;
  }
  
  real weib_surv(real t, real a, real b) {
    return exp(-(b*t)^a);
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
}

parameters {
  real<lower=0> a;
  real log_b;
}

model {
  a ~ uniform(0, 1000);
  log_b ~ normal(((log(-log(St_lower))/a - 2*log(pred_time)) +(log(-log(St_upper))/a ))/2,
					(log(-log(St_lower))/a - log(-log(St_upper))/a)/(2*1.96));
  
  for(i in 1:N) {
    t[i] ~ weib(a, exp(log_b));
  }
  
 for(i in 1:N_cens){
    target += weib_lccdf(c[i]| a, exp(log_b));
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N_pred] surv_pred;
  
  for (i in 1:N) {
    log_lik[i] <- weib_lpdf(t[i]| a, exp(log_b));
  }
  
  for(j in 1:N_pred) {
    surv_pred[j] <- weib_surv(t_new[j], a, exp(log_b));
  } 
}
