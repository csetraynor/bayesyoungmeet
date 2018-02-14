library(rstan)
funnel <- stan_demo("funnel", seed = 12345)   # has 5 divergent transitions
pairs(funnel, pars = c("y", "x[1]", "lp__"), las = 1) # below the diagonal
funnel_reparam <- stan_demo("funnel_reparam") # has no divergent transitions

sim.data <- data_frame(surv_months = t) %>%
  mutate(os_status = ifelse(is.na(surv_months), 'LIVING', 'DECEASED'),
         surv_months = ifelse(is.na(surv_months), tau[length(tau)], surv_months),
         id = seq(n), 
         censor_months = rexp(n = n, rate = 1/100))   %>%
  dplyr::mutate(os_status = ifelse(surv_months < censor_months & os_status != 'LIVING',
                                   'DECEASED', 'LIVING'
  ),
  os_months = ifelse(surv_months < censor_months  & os_status != 'LIVING',
                     surv_months, censor_months
  )
  
  
#compute log survival
  logsurv <- log(1-runif(n))
hazard_1 <- (1 - rexp(n = n, rate = 1/100))
hazard_2 <- rexp(n = n, rate = 1/100))
hazard_3 <- rexp(n = n, rate = 1/100))