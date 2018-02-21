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
  //Gaussian prior for baseline
    vector gau_prior_lp(real r1_global, real r2_global, vector ones_baseline) {
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);

    return (r1_global * sqrt(r2_global)) * ones_baseline;
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
transformed data{
  vector[T] ones_baseline;
  
    for (i in 1:T) {
    ones_baseline[i] = 1.0;
  }
}
parameters {
  // unstructured log baseline hazard for each timepoint t
  vector[T] log_baseline_raw; 
  real<lower=0> tau_s1_baseline_raw;
  real<lower=0> tau_s2_baseline_raw;
  real log_baseline_mu;
  real<lower=0> sigma_baseline;
  
  
  vector[M] beta_bg_raw; // beta for each covariate
  real<lower=0> tau_s_bg_raw;
  vector<lower=0>[M] tau_bg_raw;
}
transformed parameters {
  vector[N] log_hazard;
  vector[T] log_baseline;
  real<lower=0> tau_baseline;
  
  vector[M] beta_bg;

  //Laplace prior for beta parameter
  beta_bg = bg_prior_lp(tau_s_bg_raw, tau_bg_raw) .* beta_bg_raw;
  //Hyperprior on hazard
  for (i in 1:T){
    log_baseline[i] = log_baseline_raw[i] + log_baseline_mu + log_t_dur[i];
  }
  log_baseline = gau_prior_lp(tau_s1_baseline_raw, tau_s2_baseline_raw, ones_baseline) .* log_baseline;
  tau_baseline = 1/(sigma_baseline*sigma_baseline);
  //hazard calculation
  for (n in 1:N) {
    log_hazard[n] =  (x[n,]*beta_bg) + log_baseline[t[n]];
  }
}
model {
  beta_bg_raw ~ normal(0.0, 1.0);
  status ~ poisson_log(log_hazard);
  //Hyperprior on hazard parameters
  log_baseline_raw[1] ~ normal(0.0, 1.0);
  for (i in 2:T) {
    log_baseline_raw[i] ~ normal(log_baseline[i-1], tau_baseline);
  }
  log_baseline_mu ~ normal(0.0, 1.0);
  sigma_baseline ~ uniform(0.01, 100.0);
}

