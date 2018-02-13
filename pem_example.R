#----Load libraries---#

library(purrr)
library(httr)
library(readr)
library(survival)
library(rstan)
library(spBayesSurv)
library(pracma)
library(assertthat)
library(cgdsr)
suppressMessages(library(dplyr))
library(ggplot2)
require(ggfortify)
theme_set(theme_bw())
library(VIM)

#---Download Data----#

mycgds = CGDS("http://www.cbioportal.org/public-portal/")
study_list = getCancerStudies(mycgds)
id_sutdy = getCancerStudies(mycgds)[55,1]
case_list = getCaseLists(mycgds, id_sutdy)[2,1]
clinical_data <-  tbl_df(getClinicalData(mycgds, case_list)) 

#---- Data Cleaning ----#

clinical_data <- clinical_data %>% tibble::rownames_to_column("sample_id") 
names(clinical_data) <- tolower(names(clinical_data)) 
#convert missig values into NA
convert_blank_to_na <- function(x) {
  if(!purrr::is_character(x)){
    warning('input variate not character - return original')
    return(x)
  } else {
    ifelse(x == '', NA, x)
  }
}
clinical_data <- clinical_data %>%
  dplyr::mutate_all(funs(convert_blank_to_na))

#--Missing Data ---#

clinical_data %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

clinical_data %>% 
  filter(is.na(os_status) | os_status == "") %>%
  select(os_months, os_status) %>%
  glimpse

clinical_data %>% 
  filter(is.na(os_status) | os_status == "" |os_months < 0 | is.na(os_months)) %>%
  select(os_months, os_status) %>%
  glimpse

#--- Remove Missing Obs ---#

short_clin_dat <- 
  clinical_data %>% 
  filter(!is.na(os_status) & os_status != "" )

#confirm 44 fewer observations
assertthat::assert_that(nrow(short_clin_dat) == (nrow(clinical_data) - 44))
clinical_data <- tbl_df(short_clin_dat)
remove(short_clin_dat)

clinical_data %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

#--- Distribution event times ---#

clinical_data <- clinical_data %>%
  arrange(os_months)

clinical_data %>%
  ggplot(aes(x = os_months,
             group = os_status,
             colour = os_status,
             fill = os_status)) +
  geom_density(alpha = 0.5)

mle.surv <- survfit(Surv(os_months, os_deceased) ~ 1,
                    data = clinical_data %>%
                      mutate(os_deceased = (os_status == "DECEASED")))
require(ggfortify)
ggplot2::autoplot(mle.surv, conf.int = F) +
  ggtitle('KM survival for GGM Cohort')

#----PEM Model ----#

y0 <- c(1,1.1,0.8,0.6, 0.2)
sfun  <- stepfun(seq(20, 80, by = 20), y0, f = 0)
plot(sfun,  main = "Step function of PEM", xlab = "time (months)", ylab = "baseline hazard ratio", ylim= c(0, 1.2))

stanfile <- "'pem_biostan.stan'"
#biostan::print_stan_file(stanfile)
#or open stan file
if (interactive())
  file.edit(stanfile)

#---Testing the model on simulated data---#

pem_sim_data <- function(n, tau, beta, lambda, X, ...){
  
  #format check
  beta <- as.vector(as.numeric(beta))
  lambda<- as.vector(as.numeric(lambda))
  X <- array(matrix(as.numeric(X)), dim = c(n, length(beta)))

  #prognostic index
  mu = exp (X %*% beta )
  #extract first interval baseline hazard
  lambda0 <- lambda[1]
  #compute relative hazard for each interval respect to the first
  rel_base_risk <- lambda/lambda0
  rel_risk = lapply(mu, "*" , rel_base_risk)
  #caculate duration
  if(tau[1] != 0){
    tau <- c(0, tau)
  }
  dt = diff(tau)
  assertthat::assert_that(length(dt) == length(lambda))
  #create a helping matrix
  LD <- matrix(0, nrow = length(tau), ncol = length(rel_base_risk))
  LD[lower.tri(LD)] <- 1;
  #compute log survival
  logsurv <- log(1-runif(n))
  #compute log survival for each interval tau
  lsm = lapply(rel_risk, function(x) -lambda0 * as.vector(LD %*% (x*dt)))
  t <- (rep(NA,n))
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

set.seed(342)
test_n = 100
test_tau = c(seq(0, 900, length.out = test_n), 12*10^3)
test_baseline <- exp(-3)*runif(test_n , 0, 1)
# tau = c(seq(0, 300, by = 90), seq(300, 800, by = 100)) #time in months
# test_baseline <- exp(-4)*rev(seq(0.1, 1, by = 0.1))
X = matrix(c(rnorm(100), sample(c(0,1), 100, replace = T)), ncol=2)
test_beta = c(1.5, -2)
sim_data <-  pem_sim_data( beta = test_beta,
                           X = X,
                           tau = test_tau,
                           lambda = test_baseline,
                           n = test_n
)

## plot KM curve from simulated data
sim_data <- 
  sim_data %>%
  dplyr::mutate(os_deceased = os_status == 'DECEASED')

autoplot(survival::survfit(Surv(os_months, os_deceased) ~ 1,
                           data = sim_data
), conf.int = F) + 
  ggtitle('Simulated KM curve')
#------ long data format ----#

#set the tau interval times
tau <- sim_data %>% 
   filter(os_status == "DECEASED") %>%
  select(os_months) %>% 
  unlist() %>%
  sort()

longdata <- survival::survSplit(Surv(time = os_months, event = deceased) ~ . , 
                                cut = tau, data = (sim_data %>%
                                mutate(deceased = os_status == "DECEASED"))) %>%
                                arrange(id, os_months)
#create time point id
longdata <- longdata %>%
  group_by(id) %>%
  mutate(t = seq(n())) 

#----Generate stan data----#
M = length(beta)
gen_stan_data <- function(data){
  stan_data <- list(
    N = nrow(data),
    S = dplyr::n_distinct(data$id),
    "T" = dplyr::n_distinct(data$t),
    s = array(as.numeric(data$id)),
    t = data$t,
    M=M,
    event = data$deceased,
    obs_t = data$os_months,
    x = array(matrix(c(data$continuos, data$discrete), ncol=2), dim=c(nrow(data), 2))
  )
}

#---Set initial values---#

M = 2
gen_inits <- function() {
  list(
    beta = rcauchy(n = M, scale = 2),
    baseline = rgamma(n = n_distinct(tau), shape = 1, scale = 0.001)
    
  )
}


#-----Run Stan-------#
stanfile <- 'pem_survival_model.stan'
rstan_options(auto_write = TRUE)
simulated_fit <- stan(stanfile,
                      data = gen_stan_data(longdata),
                      init = gen_inits,
                      iter = 2000,
                      cores = min(4, parallel::detectCores()),
                      seed = 7327,
                      chains = 4
)

#----Convergence review -----#

print(simulated_fit)
rstan::traceplot(simulated_fit, 'lp__')
rstan::traceplot(simulated_fit, 'beta')
  
# if (interactive())
#   shinystan::launch_shinystan(simulated_fit)
  
  
#---Review posterior distribution of beta parameters--#


pp_beta1 <- rstan::extract(simulated_fit,'beta[1]')$beta
pp_beta2 <- rstan::extract(simulated_fit,'beta[2]')$beta

ggplot(data.frame(beta1 = pp_beta1, beta2 = pp_beta2)) + 
  geom_density(aes(x = beta1)) + 
  geom_vline(aes(xintercept = test_beta[1]), colour = 'red') +
  ggtitle('Posterior distribution of beta 1\nshowing true value in red')

ggplot(data.frame(beta1 = pp_beta1, beta2 = pp_beta2)) + 
  geom_density(aes(x = beta2)) + 
  geom_vline(aes(xintercept = test_beta[2]), colour = 'red') +
  ggtitle('Posterior distribution of beta 2\nshowing true value in red')

ggplot(data.frame(beta1 = pp_beta1, beta2 = pp_beta2)) + 
  geom_density2d(aes(x = beta1, y = beta2)) +
  geom_point(aes(x = test_beta[1], y = test_beta[2]), colour = 'red', size = 2) +
  ggtitle('Posterior distributions of beta\nshowing true parameter values in red')

#Compute probability of seeing a value beta1 >=1.5

mean(pp_beta1 >= test_beta[1])

#Compute probability of seeing a value beta1 >=-2

mean(pp_beta2 >= test_beta[2])
#Joint probability


mean(pp_beta1 >= test_beta[1] & pp_beta2 >= test_beta[2])

#---Posterior predictive checks---#

pp_beta <- as.data.frame.array(rstan::extract(simulated_fit,pars = 'beta', permuted = TRUE)$beta) 
pp_lambda <- as.data.frame.array(rstan::extract(simulated_fit,pars = 'baseline', permuted = TRUE)$baseline)

# create list
pp_beta <-  split(pp_beta, seq(nrow(pp_beta)))
pp_lambda <-  split(pp_lambda, seq(nrow(pp_lambda)))

pp_newdata <- 
  purrr::pmap(list(pp_beta),
              function(beta) {pem_sim_data(lambda = test_baseline, 
                                           beta = beta,
                                           tau = test_tau,
                                           n = test_n,
                                           X = X)
              } )

ggplot(pp_newdata %>%
         dplyr::bind_rows() %>%
         dplyr::mutate(type = 'posterior predicted values') %>%
         bind_rows(sim_data %>% dplyr::mutate(type = 'actual data'))
       , aes(x = os_months, group = os_status, colour = os_status, fill = os_status)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~type, ncol = 1)


#--Summarise posterior predictive values

## ----sim-pp-survdata-----------------------------------------------------
## cumulative survival rate for each posterior draw
pp_survdata <-
  pp_newdata %>%
  purrr::map(~ dplyr::mutate(., os_deceased = os_status == 'DECEASED')) %>%
  purrr::map(~ survival::survfit(Surv(os_months, os_deceased) ~ 1, data = .)) %>%
  purrr::map(fortify)

## summarize cum survival for each unit time (month), summarized at 95% confidence interval
pp_survdata_agg <- 
  pp_survdata %>%
  purrr::map(~ dplyr::mutate(., time_group = floor(time))) %>%
  dplyr::bind_rows() %>%
  dplyr::group_by(time_group) %>%
  dplyr::summarize(surv_mean = mean(surv)
                   , surv_p50 = median(surv)
                   , surv_lower = quantile(surv, probs = 0.025)
                   , surv_upper = quantile(surv, probs = 0.975)
  ) %>%
  dplyr::ungroup()

## km-curve for test data 
test_data_kmcurve <- 
  fortify(
    survival::survfit(
      Surv(os_months, os_deceased) ~ 1, 
      data = sim_data %>% 
        dplyr::mutate(os_deceased = os_status == 'DECEASED')
    )) %>%
  dplyr::mutate(lower = surv, upper = surv)

ggplot(pp_survdata_agg %>%
         dplyr::mutate(type = 'posterior predicted values') %>%
         dplyr::rename(surv = surv_p50, lower = surv_lower, upper = surv_upper, time = time_group) %>%
         bind_rows(test_data_kmcurve %>% dplyr::mutate(type = 'actual data')),
       aes(x = time, group = type, linetype = type)) + 
  geom_line(aes(y = surv, colour = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  xlim(c(0, 250))

#----Fitting Model to TCGA Glioblastome data----#

#--- Update gen stan data function to include covariates ---#

#--- This function will take a formula object as input --- #

gen_stan_data <- function(data, formula = as.formula(~1)) {
  if(!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  data <- tibble::rownames_to_column(data, var = "id")
  
  tau <- data %>% 
    filter(os_status == "DECEASED") %>%
    select(os_months) %>% 
    unlist() %>%
    unique() %>%
    sort()
  
  longdata <- tbl_df(survival::survSplit(Surv(time = os_months, event = deceased) ~ . , 
                                  cut = tau,
                                  data = (data %>%
                                          mutate(deceased = os_status == "DECEASED"))) %>%
                                          arrange(id, os_months))
  
  #create time point id
  longdata <- longdata %>%
    group_by(id) %>%
    mutate(t = seq(n())) 
  
  X_bg <- longdata %>%
    model.matrix(formula, data = .)
  
  M_bg <- ncol(X_bg)
  
  if (M_bg > 1){
    if("(Intercept)" %in% colnames(X_bg))
      X_bg <- array(X_bg[,-1], dim = c(nrow(longdata), M_bg - 1))
    M_bg <- ncol(X_bg)
  }
  
  
  stan_data <- list(
    N = nrow(longdata),
    S = dplyr::n_distinct(longdata$id),
    "T" = dplyr::n_distinct(longdata$t),
    s = as.numeric(longdata$id),
    t = longdata$t,
    M = M_bg,
    event = longdata$deceased,
    obs_t = longdata$os_months,
    x = X_bg
  )
}


#---Update initial values---#
tau <- clinical_data %>% 
  filter(os_status == "DECEASED") %>%
  select(os_months) %>% 
  unlist() %>%
  unique() %>%
  sort()


M = 2
gen_inits <- function() {
  list(
    beta = rcauchy(n = M, scale = 2),
    baseline = rgamma(n = n_distinct(tau), shape = 1, scale = 0.001)
    
  )
}


stanfile <- "pem_survival_model.stan"

stan_fit1 <- rstan::stan(stanfile,
                       data = gen_stan_data(clinical_data, '~ I(sex == "Male") + age'),
                       init = gen_inits,
                       iter = 2000,
                       cores = min(4, parallel::detectCores()),
                       seed = 7327,
                       chains = 4
)


###########




















pp_predict_surv <- function(pp_beta, pp_lambda, n, tau, X,
                            level = 0.9, 
                            plot = F, data = NULL,
                            sim_data_fun = pem_sim_data
) {
  pp_newdata <- 
    purrr::pmap(list(pp_beta, pp_lambda),
                function(a, b) {pem_sim_data(lambda = a, 
                                             beta = b,
                                             tau = tau,
                                             n = n,
                                             X = X)
                } )
  
  pp_survdata <-
    pp_newdata %>%
    purrr::map(~ dplyr::mutate(., os_deceased = os_status == 'DECEASED')) %>%
    purrr::map(~ survival::survfit(Surv(os_months, os_deceased) ~ 1, data = .)) %>%
    purrr::map(fortify)
  
  ## compute quantiles given level 
  lower_p <- 0 + ((1 - level)/2)
  upper_p <- 1 - ((1 - level)/2)
  
  pp_survdata_agg <- 
    pp_survdata %>%
    purrr::map(~ dplyr::mutate(.,
                               time_group = floor(time))) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(time_group) %>%
    dplyr::summarize(surv_mean = mean(surv)
                     , surv_p50 = median(surv)
                     , surv_lower = quantile(surv,
                                             probs = lower_p)
                     , surv_upper = quantile(surv,
                                             probs = upper_p)
    ) %>%
    dplyr::ungroup()
  
  if (plot == FALSE) {
    return(pp_survdata_agg)
  } 
  
  ggplot_data <- pp_survdata_agg %>%
    dplyr::mutate(type = 'posterior predicted values') %>%
    dplyr::rename(surv = surv_p50,
                  lower = surv_lower,
                  upper = surv_upper, time = time_group)
  
  if (!is.null(data))
    ggplot_data <- 
    ggplot_data %>% 
    bind_rows(
      fortify(
        survival::survfit(
          Surv(os_months, os_deceased) ~ 1, 
          data = data %>% 
            dplyr::mutate(
              os_deceased = os_status == 'DECEASED')
        )) %>%
        dplyr::mutate(lower = surv,
                      upper = surv, type = 'actual data')
    )
  
  pl <- ggplot(ggplot_data,
               aes(x = time, group = type, linetype = type)) + 
    geom_line(aes(y = surv, colour = type)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2)
  
  pl 
}



