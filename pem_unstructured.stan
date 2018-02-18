/*  Variable naming:
 // dimensions
 N          = total number of observations (length of data)
 S          = number of sample ids
 T          = max timepoint (number of timepoint ids)
 M          = number of covariates

 // main data matrix (per observed timepoint*record)
 s          = sample id for each obs
 t          = timepoint id for each obs
 event      = integer indicating if there was an event at time t for sample s
 x          = matrix of real-valued covariates at time t for sample n [N, X]

 // timepoint-specific data (per timepoint, ordered by timepoint id)
 t_obs      = observed time since origin for each timepoint id (end of period)
 t_dur      = duration of each timepoint period (first diff of t_obs)

*/
// Carlos Serra Traynor <carlos.serra91@gmail.com>

  
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
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> T;
  int<lower=0> M;

  // data matrix
  int<lower=1, upper=N> s[N];     // sample id
  int<lower=1, upper=T> t[N];     // timepoint id
  int<lower=0, upper=1> event[N]; // 1: event, 0:censor
  matrix[N, M] x;                 // explanatory vars
  real<lower=0> obs_t[N];         // observed end time for each obs

}
transformed data {
    // timepoint data
  vector<lower=0>[T] t_obs;
  vector<lower=0>[T] t_dur;
  vector[T] log_t_dur;  // log-duration for each timepoint
  int n_trans[S, T];
  
    // capture observation time for each timepoint id t
  for (i in 1:N) {
      // assume these are constant per id across samples
      t_obs[t[i]] = obs_t[i];
  }
  
  // duration of each timepoint
  // duration at first timepoint = t_obs[1] ( implicit t0 = 0 )
  t_dur[1] = t_obs[1];
  for (i in 2:T) {
      t_dur[i] = t_obs[i] - t_obs[i-1];
  }

  log_t_dur = log(t_dur);
  
    // n_trans used to map each sample*timepoint to n (used in gen quantities)
  // map each patient/timepoint combination to n values
  for (n in 1:N) {
      n_trans[s[n], t[n]] = n;
  }

  // fill in missing values with n for max t for that patient
  // ie assume "last observed" state applies forward (may be problematic for TVC)
  // this allows us to predict failure times >= observed survival times
  for (samp in 1:S) {
      int last_value;
      last_value = 0;
      for (tp in 1:T) {
          // manual says ints are initialized to neg values
          // so <=0 is a shorthand for "unassigned"
          if (n_trans[samp, tp] <= 0 && last_value != 0) {
              n_trans[samp, tp] = last_value;
          } else {
              last_value = n_trans[samp, tp];
          }
      }
  }

}
parameters {
    vector[T] log_baseline_raw; // unstructured baseline hazard for each timepoint t
  real<lower=0> baseline_sigma;
  real log_baseline_mu;
  real<lower=0> tau_s_bg_raw;
  vector<lower=0>[M] tau_bg_raw; // beta for each covariate
  vector[M] beta_bg_raw;
}
transformed parameters {
  vector[M] beta_bg;
  vector[N] log_hazard;
  vector[T] log_baseline;//unstructured baseline hazard for each timepoint t

  beta_bg = bg_prior_lp(tau_s_bg_raw, tau_bg_raw) .* beta_bg_raw;
  log_baseline = log_baseline_mu + log_baseline_raw + log_t_dur;

  for (n in 1:N) {
    log_hazard[n] = log_baseline[t[n]] + x[n,]*beta_bg;
  }
}
model {
  beta_bg_raw ~ normal(0.0, 1.0);                                             event ~ poisson_log(log_hazard);
  // Prior on hazard parameters
  log_baseline_mu ~ normal(0, 1);
  baseline_sigma ~ lognormal(0, 2);
  log_baseline_raw[1] ~ normal(0, baseline_sigma);
  	for(n in 2:T) {
		log_baseline_raw[n] ~ normal(log_baseline_raw[n-1],baseline_sigma);
	 }
}
generated quantities {
 real log_lik[N];
  vector[T] baseline;
  real y_hat_time[S];      // predicted failure time for each sample
  int y_hat_event[S];      // predicted event (0:censor, 1:event)

  // compute raw baseline hazard, for summary/plotting
  baseline = exp(log_baseline_mu + log_baseline_raw);

  // prepare log_lik for loo-psis
  for (n in 1:N) {
      log_lik[n] = poisson_log_log(event[n], log_hazard[n]);
  }

  // posterior predicted values
  for (samp in 1:S) {
      int sample_alive;
      sample_alive = 1;
      for (tp in 1:T) {
        if (sample_alive == 1) {
              int n;
              int pred_y;
              real log_haz;

              // determine predicted value of this sample's hazard
              n = n_trans[samp, tp];
              log_haz = log_baseline[tp] + x[n,] * beta_bg;

              // now, make posterior prediction of an event at this tp
              if (log_haz < log(pow(2, 30)))
                  pred_y = poisson_log_rng(log_haz);
              else
                  pred_y = 9;

              // summarize survival time (observed) for this pt
              if (pred_y >= 1) {
                  // mark this patient as ineligible for future tps
                  // note: deliberately treat 9s as events
                  sample_alive = 0;
                  y_hat_time[samp] = t_obs[tp];
                  y_hat_event[samp] = 1;
              }

          }
      } // end per-timepoint loop

      // if patient still alive at max
      if (sample_alive == 1) {
          y_hat_time[samp] = t_obs[T];
          y_hat_event[samp] = 0;
      }
  } // end per-sample loop
}
