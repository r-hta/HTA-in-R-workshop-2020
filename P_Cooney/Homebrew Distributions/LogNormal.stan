functions{
real lognormal_surv(real t, real a, real b){
  return 1 - lognormal_cdf(t, a, b);
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
  real a;
  real b;
}

model {

  a ~ normal((2*log(pred_time) - b*(inv_Phi(1-St_lower) + inv_Phi(1-St_upper)))/2,
				(b*(inv_Phi(1-St_lower) - inv_Phi(1-St_upper)))/(2*1.96));
  b ~ cauchy(0, 7);
  
  for(i in 1:N) {
    t[i] ~ lognormal(a, b);
  }
  
 for(i in 1:N_cens){
    target += lognormal_lccdf(c[i]| a, b);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N_pred] surv_pred;
  
  for (i in 1:N) {
    log_lik[i] = lognormal_lpdf(t[i]| a, b);
  }
  
  for(j in 1:N_pred) {
    surv_pred[j] =  lognormal_surv(t_new[j], a, b);
  }
}
