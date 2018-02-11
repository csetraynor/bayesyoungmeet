library(dplyr)

pem_sim_data <- function(n, tau, bh){
  #extract first interval baseline hazard
  bh0 <- bh[1]
  #compute relative baseline hazard for each interval respect to the first
  rbh <- bh/bh0
  #caculate duration
  dt = diff(tau)
  #create a helping matrix
  LD <- matrix(0, nrow = length(tau), ncol = length(rbh)); LD[lower.tri(LD)] <- 1;
  #compute log survival
  logsurv <- log(1-runif(n))
  #compute log survival for each interval tau
  lsm <- -bh0 * as.vector(LD %*% (rbh*dt))
  #find appropiate time interval
  t <- rep(NA,n)
  for (i in 1:length(h0+1)) {
    t <- ifelse(lsm[i]>=logsurv & logsurv>lsm[i+1], tau[i] + (lsm[i] - logsurv)/h0/rbh[i], t)
  }
  
  sim.data <- data_frame(surv_months = t) %>%
    mutate(
           id = seq(n), 
           censor_months = rexp(n = n, rate = exp(-4)))%>%
    dplyr::mutate(os_status = ifelse(surv_months < censor_months,
                                     'DECEASED', 'LIVING'
    ),
    os_months = ifelse(surv_months < censor_months,
                       surv_months, censor_months
    )
    )
  return(sim.data)
}



pem_sim_data(n = 100, bh = exp(-0.9)* c(1, 1.01, 0.381, 0.150), tau = c(0, 22, 35, 74, 102))

