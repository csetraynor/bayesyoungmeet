/*Null Model with parametric Weibull 

  Variable naming: 
  obs       = observed 
  cen       = (right) censored 
  N         = number of samples 
  tau       = scale parameter 
*/ 
data { 
 int<lower=0> Nobs; //number of observed data points
 int<lower=0> Ncen; //number of censored data points
 vector[Nobs] yobs; //time to observed events
 vector[Ncen] ycen; //time to censored events
} 
  
transformed data { 
 real<lower=0> tau_mu; 
 real<lower=0> tau_al; //constant scaling term
  
 tau_mu = 10.0; 
 tau_al = 10.0; 
} 
  
parameters { 
 real alpha_raw; //parameter with a normal prior 
 real mu; 
} 
  
transformed parameters { 
 real alpha; 
 alpha = exp(tau_al * alpha_raw); //non-centered parameterization
} 
  
model { 
 yobs ~ weibull(alpha, exp(-(mu)/alpha)); 
 target += weibull_lccdf(ycen | alpha, exp(-(mu)/alpha)); 
  
 alpha_raw ~ normal(0.0, 1.0); 
 mu ~ normal(0.0, tau_mu); 
}
