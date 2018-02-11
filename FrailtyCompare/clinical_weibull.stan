/*
  N = number of observations (individuals)
  M = number of covariates
  x = covariates matrix
  y = obsrved survival time
  event = censor indicator (1:event, 0:censor)
  */
functions {
  vector sqrt_vec(vector x) {
    vector[dims(x)[1]] res;
    
    for (m in 1:dims(x)[1]){
      res[m] = sqrt(x[m]);
    }
    
    return res;
  }
  
  vector bg_prior_lp(real r_global, vector r_local) {
    r_global ~ normal(0.0, 10.0);
    r_local ~ inv_chi_square(1.0);
    
    return r_global * sqrt_vec(r_local);
  }
}
data {
  // dimensions
  int<lower=0> N;             // number of observations
  int<lower=1> M;             // number of predictors
  
  // observations
  matrix[N, M] x;             // predictors for observation n
  vector[N] y;                // time for observation n
  vector[N] event;            // event status (1:event, 0:censor) for obs n
}
transformed data {
  real<lower=0> tau_mu;
  real<lower=0> tau_al;
  tau_mu = 10.0;
  tau_al = 10.0;
}
parameters {
  real<lower=0> tau_s_raw;
  vector<lower=0>[M] tau_raw;
  
  real alpha_raw;
  vector[M] beta_raw;
  
  real mu;
}

transformed parameters {
  vector[M] beta;
  real<lower=0> alpha;
  vector[N] lp;
  beta = bg_prior_lp(tau_s_raw, tau_raw) .* beta_raw;
  alpha = exp(tau_al * alpha_raw);
  for (n in 1:N) {
    lp[n] = mu + dot_product(x[n], beta);
  }
}
model {
  // priors
  target += normal_lpdf(beta_raw | 0.0, 1.0);
  target += normal_lpdf(alpha | 0.0, 1.0);
  target += normal_lpdf(mu | 0.0, tau_mu)
  print("A");
  // likelihood
  for (n in 1:N) {
    if (event[n]==1)
      target += weibull_lpdf(y[n] | alpha, exp(-(lp[n])/alpha));
    else
      target += weibull_lccdf(y[n] | alpha, exp(-(lp[n])/alpha));
  }
}
generated quantities {
  vector<lower=0, upper=1>[N] yhat_uncens;
  vector<upper=0>[N] log_lik;
  
  for (n in 1:N) {
    yhat_uncens[n] = weibull_rng(alpha, exp(-(lp[n])/alpha));
    if (event[n]==1) {
      log_lik[n] = weibull_lpdf(y[n] | alpha, exp(-(lp[n])/alpha));
    } else {
      log_lik[n] = weibull_lccdf(y[n] | alpha, exp(-(lp[n])/alpha));
    }
  }
}
