## create initial estimates
init_wei <- function() list(
  tau_s_raw = abs(rnorm(1)),
  tau_raw = abs(rnorm(M)),
  alpha_raw = rnorm( 1, sd = 0.1),
  beta_raw = rnorm(M),
  mu = rnorm(1)
)

## Specify the variables for which you want history and density plots
parametersToPlot <- c("beta_raw", "alpha_raw", "mu")

## Additional variables to monitor
otherRVs <- c( "log_lik", "yhat_uncens")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

library(rstan)
nchains <- 1
niter <- 1000
nwarmup <- 500

fit_wei_clinical_bg <- stan(file = "clinical_weibull.stan", 
                            data = data, 
                            pars = parameters,
                            init = init_wei, 
                            chains = nchains, 
                            iter = niter,
                            warmup = nwarmup,
                            control = list(stepsize = 0.01, adapt_delta = 0.99),
                            cores = min(nchains, parallel::detectCores())) 