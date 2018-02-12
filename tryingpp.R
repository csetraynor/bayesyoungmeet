pem_sim_data <- function(n, tau, lambda, beta, X, ...){
  
# str(estimates)
# 
# lambda <- estimates[[1]]
# beta <- estimates[[2]]

beta <- as.vector(as.numeric(beta))
lambda<- as.vector(as.numeric(lambda))

print("beta =");
str(beta)
print("lambda")
str(lambda)

#prognostic index
mu = exp (X %*% beta )
#extract first interval baseline hazard
lambda0 <- lambda[1]
#compute relative hazard for each interval respect to the first
rel_base_risk <- lambda/lambda0
rel_risk = lapply(mu, "*" , rel_base_risk)
#caculate duration
dt = diff(c(0,tau))
# or dt = diff(c(0,tau, Inf))
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


pp_beta <- as.data.frame.array(rstan::extract(simulated_fit,pars = 'beta', permuted = TRUE)$beta) 

pp_lambda <- as.data.frame.array(rstan::extract(simulated_fit,pars = 'baseline', permuted = TRUE)$baseline)

# create list
pp_beta <-  split(pp_beta, seq(nrow(pp_beta)))
pp_lambda <-  split(pp_lambda, seq(nrow(pp_lambda)))

pp_estimates <- list(pp_lambda,pp_beta)

pp_newdata <- 
  purrr::pmap(list(pp_lambda, pp_beta),
              function(a,b) {pem_sim_data(lambda = a, 
                                          beta = b,
                                          tau = tau,
                                          n = test_n,
                                          X = X)
                } )


#pp_lambda <- matrix(unlist(pp_lambda), ncol = length(tau), byrow = TRUE)










set.seed(342)
test_n = 100
X = matrix(c(rnorm(100), sample(c(0,1), 100, replace = T)), ncol=2)
tau <- sim_data %>% 
  # filter(os_status == "DECEASED") %>%
  select(os_months) %>% 
  unlist() %>%
  sort()



pp_newdata <- 
  purrr::map2(.x = pp_beta,
              .y = pp_lambda,
              .f = ~ pem_sim_data(beta = .x,
                                  lambda = .y,
                                  tau = tau,
                                  n = test_n,
                                  X = X
              )
  )



pp_newdata <- mapply(function(theta){pem_sim_data(estimates = theta, 
                                              tau = tau,
                                              n = test_n,
                                              X = X)},
                     theta = pp_estimates)



