/*  Variable naming:
  // dimensions
N          = total number of observations (length of data)
S          = number of sample ids
T          = max timepoint (number of timepoint ids)
M          = number of covariates

// data
s          = sample id for each obs
t          = timepoint id for each obs
event      = integer indicating if there was an event at time t for sample s
x          = matrix of real-valued covariates at time t for sample n [N, X]
obs_t      = observed end time for interval for timepoint for that obs

*/
  
data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=0> M;
  int<lower=1, upper=N> s[N];     // sample id
  int<lower=1> t[N];     // timepoint id
  int<lower=0, upper=1> status[N]; // 1: event, 0:censor
  matrix[N, M] x;                 // explanatory vars
  real t_dur[T];
}
parameters {
  real sigma_baseline;
  vector[T] log_baseline; // unstructured log baseline hazard for each timepoint t
  real alpha;
  real<lower=0> sigma_beta;
  vector[M] beta; // beta for each covariate
}
transformed parameters {
  vector[N] hazard;
  real<lower=0> tau_baseline;
  //Hyperprior on historical variance parameters
  tau_baseline = 1/(sigma_baseline*sigma_baseline);
  //hazard calculation
  for (n in 1:N) {
    hazard[n] = t_dur[t[n]] * exp((x[n,]*beta) + log_baseline[t[n]]);
  }
}
model {
  beta ~ normal(alpha, sigma_beta);
  alpha ~ normal(0,1);
  sigma_beta ~ cauchy(0,1);
  status ~ poisson(hazard);
  //Prior on hazard parameters
  log_baseline[1] ~ normal(0, 0.0001);
  for (i in 2:T) {
    log_baseline[i] ~ normal(log_baseline[i-1], tau_baseline);
  }
  //Hyperprior on historical variance parameters
  sigma_baseline ~ uniform(0.01, 100);
}

