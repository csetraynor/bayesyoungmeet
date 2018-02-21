#---Set initial values---#
gen_inits <- function(M) {
  function() 
    list(
      log_baseline_raw = rnorm(n= length(t_dur)),
      log_baseline_mu = rnorm(1),
      sigma_baseline = runif(1, 0.01, 100),
      beta_bg_raw = rnorm(M),
      tau_s_bg_raw = 0.1*abs(rnorm(1)),
      tau_bg_raw = abs(rnorm(M)),
      tau_s1_baseline_raw = 0.1*abs(rnorm(1)),
      tau_s2_baseline_raw = 0.1*abs(rnorm(1))
    )
}
#-----Run Stan-------#
nChain <- 2
stanfile <- 'pem_bg_gaus.stan'
rstan_options(auto_write = TRUE)
test_simulated <- stan(stanfile,
                       data = gen_stan_data(longdata),
                       init = gen_inits(M=2),
                       iter = 1000,
                       cores = min(nChain, parallel::detectCores()),
                       seed = 7327,
                       chains = nChain,
                       #  control = list(adapt_delta = 0.95),
                       pars = c("beta_bg", "log_baseline", "lp__")
)
