pem_sim_data <- function(n, tau, lambda, beta, X){
  
  #prognostic index
  mu = exp (X[,1:2] %*% beta )
  #extract first interval baseline hazard
  lambda0 <- lambda[1]
  #compute relative hazard for each interval respect to the first
  rel_base_risk <- lambda/lambda0
  rel_risk = lapply(mu, "*" , rel_base_risk)
  #caculate duration
  dt = c(diff(tau))
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
           os_months = ifelse(is.na(surv_months), tau[length(tau)], surv_months),
           id = seq(n))   %>%   cbind(X) %>%
    rename("continuos var" = "1", "discrete var" = "2")
  
  return(sim.data)
}
set.seed(342)
test_n = 100
sim_data <-  pem_sim_data(n = test_n, lambda = exp(-3)*rev(seq(0.1,1, by = 0.1)), tau = seq(0, 200, by = 20), X = matrix(c(rnorm(100), sample(c(0,1), 100, replace = T)), ncol=2), beta = c(0.5,1))

