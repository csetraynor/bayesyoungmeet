##Rainer Walke
pcbhsim <- function(){
  # set the start population size for both groups
  number1 <- 108
  number2 <- 34
  # set the time points
  TAU <- c(0, 1, 4.33, 26, 52)
  DT = TAU[2:5] - TAU[1:4]
  # set the absolute risk
  h0 <- exp(-2.01)
  # set the relative risks
  G <- c(1, 1.01, 0.381, 0.150)
  # set the relative risks for the population groups
  P <- c(1, 1.54)
  # create a helping matrix
  LD <- matrix(0,nrow=5, ncol=4)
  LD[lower.tri(LD)]<-1
  
  
  # start with population group 1
  LS <- log(1-runif(number1))
  GP <- P[1]*G
  # determine the ln(S) for all TAU
  LSM <- -h0 * as.vector(LD %*% (GP * DT))
  # find the appropriate time interval
  t1 <- rep(NA,number1)
  for (i in 1:4) {
    t1 <- ifelse(LSM[i]>=LS & LS>LSM[i+1], TAU[i] + (LSM[i] - LS)/h0/GP[i], t1)
  }
  # end of population group1
  
  # start with population group 2
  LS <- log(1-runif(100))
  GP <- P[2]*G
  # determine the ln(S) for all TAU
  LSM <- -h0 * as.vector(LD %*% (GP * DT))
  # find the appropriate time interval
  t2 <- rep(NA,number1)
  for (i in 1:4) {
    t2 <- ifelse(LSM[i]>=LS & LS>LSM[i+1], TAU[i] + (LSM[i] - LS)/h0/GP[i], t2)
  }
  # end of population group 2
  
  # combine both populations
  sim.data <- data.frame( rbind(cbind(t=t1,period=1),cbind(t2,2)))
  sim.data$occ <- ifelse(is.na(sim.data$t), 0, 1)
  sim.data$t <- ifelse(is.na(sim.data$t), TAU[5], sim.data$t)
  sim.data$id <- row(sim.data)[,1]
  
  return(sim.data)
}

#compute prognostic index
x <- as.matrix(rnorm(100), sample(c(0,1), 100, replace = T), ncol=2)
beta <- 1.35
pi = as.list( exp(cumsum(x*beta)))
gp = lapply(pi, "*" , rbh)

logsurv <- log(1-runif(100))

#compute log survival for each interval tau
lsm = lapply(gp, function(x) -bh0 * as.vector(LD %*% (x*dt)))

#find appropiate time interval
t <- (rep(NA,100))
t = mapply(function(x, y, z) {
  for (i in 1:length(rbh)) {
    t = ifelse(x[i]>=z & z>x[i+1], tau[i] + (x[i] - z)/bh0/y[i], t)
  }
  #t = do.call(rbind, t)
  return(t)
} , x = lsm, y = gp , z = as.list(logsurv)
)


t = do.call(rbind, t)


LSM <- -h0 * as.vector(LD %*% (GP * DT))
# find the appropriate time interval
t2 <- rep(NA,number1)
for (i in 1:4) {
  t2 <- ifelse(LSM[i]>=LS & LS>LSM[i+1], TAU[i] + (LSM[i] - LS)/h0/GP[i], t2)
}


pem_sim_data <- function(n, tau, bh, x, beta){
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
  #compute prognostic index
  
  
  #compute log survival for each interval tau
  lsm <- -bh0 * as.vector(LD %*% (rbh*dt))
  #find appropiate time interval
  t <- rep(NA,n)
  for (i in 1:length(bh)) {
    t <- ifelse(lsm[i]>=logsurv & logsurv>lsm[i+1], tau[i] + (lsm[i] - logsurv)/bh0/rbh[i], t)
  }
  
  sim.data <- data_frame(surv_months = t) %>%
    mutate(os_status = ifelse(is.na(surv_months), 'LIVING', 'DECEASED'),
           surv_months = ifelse(is.na(surv_months), tau[length(tau)], surv_months),
           id = seq(n), 
           censor_months = rexp(n = n, rate = exp(-2)))%>%
    dplyr::mutate(os_status = ifelse(surv_months < censor_months & os_status != 'LIVING',
                                     'DECEASED', 'LIVING'
    ),
    os_months = ifelse(surv_months < censor_months  & os_status != 'LIVING',
                       surv_months, censor_months
    )
    )
  return(sim.data)
}



sim_data <-  pem_sim_data(n = 100, bh = exp(-2.9)* c(1, 0.8, 0.2, 0.25, 0.125), tau = c(0, 10, 15, 30, 50, 100))


