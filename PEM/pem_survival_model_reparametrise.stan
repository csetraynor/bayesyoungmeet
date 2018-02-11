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
// Jacqueline Buros Novik <jackinovik@gmail.com>

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=0> M;
  int<lower=1, upper=N> s[N];     // sample id
  int<lower=1, upper=T> t[N];     // timepoint id
  int<lower=0, upper=1> event[N]; // 1: event, 0:censor
  matrix[N, M] x;                 // explanatory vars
  real<lower=0> t_obs[N];         // observed end time for each obs
  real<lower=0> t_dur[N];         // duration for each timepoint
}
transformed data {
  real c;
  real r;

  // baseline hazard params (fixed)
  c = 0.001;
  r = 0.1;
}
parameters {
  vector<lower=0>[T] baseline; // unstructured baseline hazard for each timepoint t
  vector[M_bg] beta_bg_raw; // beta for each covariate
}
transformed parameters {
  vector<lower=0>[N] hazard;
  vector[M_bg] beta_bg;
  
  beta_bg = bg_prior_lp(tau_s_bg_raw, tau_bg_raw) .* beta_bg_raw;

  for (n in 1:N) {
    hazard[n] = exp(x[n,]*beta_bg_raw)*baseline[t[n]];
  }
}
model {
  for (i in 1:T) {
      baseline[i] ~ gamma(r * t_dur[i] * c, c);
  }
  beta_bg_raw ~ cauchy(0, 2);
  event ~ poisson(hazard);
}
generated quantities {
  real log_lik[N];

  for (i in 1:N) {
      log_lik[i] = poisson_lcdf(event[i] | hazard[i]);
  }
}
