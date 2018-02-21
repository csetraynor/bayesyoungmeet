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
  
functions {
  vector sqrt_vec(vector x) {
    vector[dims(x)[1]] res;
      
    for (m in 1:dims(x)[1]){
      res[m] = sqrt(x[m]);
    }
    return res;
  }
  //Laplace prior for beta parameters
  vector bg_prior_lp(real r_global, vector r_local) {
    r_global ~ normal(0.0, 10.0);
    r_local ~ inv_chi_square(1.0);
    
    return r_global * sqrt_vec(r_local);
  }
}
data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=0> M;
  int<lower=1, upper=N> s[N];     // sample id
  int<lower=1> t[N];     // timepoint id
  int<lower=0, upper=1> status[N]; // 1: event, 0:censor
  matrix[N, M] x;                 // explanatory vars
  vector[T] log_t_dur;
}
parameters {
  vector[T] log_baseline_raw; // unstructured log baseline hazard for each timepoint t
  real<lower=0> sigma_baseline;
  real log_baseline_mu;
  vector[M] beta_bg_raw; // beta for each covariate
  real<lower=0> tau_s_bg_raw;
  vector<lower=0>[M] tau_bg_raw;
}
transformed parameters {
  vector[N] log_hazard;
  vector[T] log_baseline;
  vector[M] beta_bg;
  real<lower=0> tau_baseline;
  //Laplace prior for beta parameter
  beta_bg = bg_prior_lp(tau_s_bg_raw, tau_bg_raw) .* beta_bg_raw;
  //Hyperprior on hazard
  for (i in 1:T){
    log_baseline[i] = log_baseline_raw[i]  + log_baseline_mu;
  } 
  tau_baseline = 1/(sigma_baseline*sigma_baseline);
  //hazard calculation
  for (n in 1:N) {
    log_hazard[n] =  (x[n,]*beta_bg) + log_baseline[t[n]] + log_t_dur[t[n]];
  }
}
model {
  beta_bg_raw ~ normal(0.0, 1.0);
  status ~ poisson_log(log_hazard);
  //Prior on hazard parameters
  log_baseline_raw[1] ~ normal(0.0, 1.0);
  for (i in 2:T) {
    log_baseline_raw[i] ~ normal(log_baseline[i-1], tau_baseline);
  }
  //Hyperprior on historical variance parameters
  log_baseline_mu ~ normal(0.0, 1.0);
  sigma_baseline ~ uniform(0.01, 100);
}

