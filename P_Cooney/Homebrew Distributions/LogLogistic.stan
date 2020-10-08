//http://reliawiki.org/index.php/The_Loglogistic_Distribution

functions{
real logistic_surv(real t, real a, real b){
  return 1/(1+exp((log(t)-a)/b));
}

real logistic_surv_lpdf(real t, real a, real b){
  return log(exp((log(t)-a)/b)/(b*t*((1+exp((log(t)-a)/b))^2)));

}

real logistic_surv_lccdf(real t, real a, real b){
 return log(1/(1+exp((log(t)-a)/b)));
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
  real<lower=0> b;
}

model {
  a ~ normal((2*log(pred_time) - b*log((1/St_upper)-1) - b*log((1/St_lower)-1))/2,
			(-b*log((1/St_upper)-1) + b*log((1/St_lower)-1))/(2*1.96));
  b ~ cauchy(0,7);

  for(i in 1:N) {
    t[i] ~ logistic_surv(a, b);
  }
  
 for(i in 1:N_cens){
    target += logistic_surv_lccdf(c[i]| a, b);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N_pred] surv_pred;
  
  for (i in 1:N) {
    log_lik[i] = logistic_surv_lpdf(t[i]| a, b);
  }
  
  for(j in 1:N_pred) {
    surv_pred[j] =  logistic_surv(t_new[j], a, b);
  }
}
