pem_sim_data <- function(n, tau, estimates, X){
  
str(estimates)

beta <- estimates[[2]]
lambda <- estimates[[1]]

print("beta =");
str(beta)
print("lambda")
str(lambda)

#prognostic index
mu = exp (X[,1:2] %*% beta )
#extract first interval baseline hazard
lambda0 <- lambda[1]
#compute relative hazard for each interval respect to the first
rel_base_risk <- lambda/lambda0
rel_risk = lapply(mu, "*" , rel_base_risk)
#caculate duration
dt = diff(c(0,tau, Inf))
assertthat::assert_that(length(dt) == length(lambda))
#create a helping matrix
LD <- matrix(0, nrow = length(tau), ncol = length(rel_base_risk))
LD[lower.tri(LD)] <- 1;
#compute log survival
logsurv <- log(1-runif(n))
#compute log survival for each interval tau
lsm = lapply(rel_risk, function(x) -lambda0 * as.vector(LD %*% (x*dt)))
t <- (rep(NA,100))
#find appropiate time interval
t = mapply(function(x, y, z) {
  for (i in seq_along(lambda)) {
    t = ifelse(x[i]>=z & z>x[i+1], tau[i] + (x[i] - z)/lambda0/y[i], t)
  }
  return(t)
} , x = lsm, y = rel_risk , z = as.list(logsurv)
)

sim.data <- data_frame(surv_months = t) %>%
  mutate(os_status = ifelse(is.na(surv_months), 'LIVING', 'DECEASED'),
         surv_months = ifelse(is.na(surv_months), tau[length(tau)], surv_months),
         id = seq(n), 
         censor_months = rexp(n = n, rate = 1/100))%>%
  dplyr::mutate(os_status = ifelse(surv_months < censor_months & os_status != 'LIVING',
                                   'DECEASED', 'LIVING'
  ),
  os_months = ifelse(surv_months < censor_months  & os_status != 'LIVING',
                     surv_months, censor_months
  )
  ) %>%   cbind(X) %>%
  rename("continuos" = "1", "discrete" = "2")

return(sim.data)
}


pp_beta <- rstan::extract(simulated_fit,pars = 'beta', permuted = TRUE)$beta
pp_lambda <- rstan::extract(simulated_fit,pars = 'baseline', permuted = TRUE)$baseline

#pp_lambda <- matrix(unlist(pp_lambda), ncol = length(tau), byrow = TRUE)

pp_estimates <- list(pp_lambda,pp_beta)

set.seed(342)
test_n = 100
X = matrix(c(rnorm(100), sample(c(0,1), 100, replace = T)), ncol=2)


pp_newdata <- mapply(function(theta){pem_sim_data(estimates = theta, 
                                              tau = tau,
                                              n = test_n,
                                              X = X)},
                     theta = pp_estimates)


pp_newdata <- 
  purrr::pmap(.l = pp_estimates,
              .f = ~ pem_sim_data(estimates = .l, 
                                  tau = tau,
                                  n = test_n,
                                  X = X
              )
  )
