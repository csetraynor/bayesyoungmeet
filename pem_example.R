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
library(caret)
library(mice)
library(missMDA)
#---Download Data----#
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
study_list = getCancerStudies(mycgds)
id_sutdy = getCancerStudies(mycgds)[56,1]
case_list = getCaseLists(mycgds, id_sutdy)[2,1]
clinical_data <-  tbl_df(getClinicalData(mycgds, case_list)) 

#---- Data Cleaning ----#
clinical_data <- clinical_data %>% tibble::rownames_to_column("sample_id") 
clinical_data <- clinical_data %>% tibble::rownames_to_column("num_id") 
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
  filter(is.na(os_status) | os_status == "") %>%
  select(os_months, os_status) %>%
  glimpse
clinical_data %>% 
  filter(is.na(os_status) | os_status == "" |os_months < 0 | is.na(os_months)) %>%
  select(os_months, os_status) %>%
  glimpse
clinical_data %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)
#--- Distribution event times ---#
clinical_data %>%
  ggplot(aes(x = os_months,
             group = os_status,
             colour = os_status,
             fill = os_status)) +
  geom_density(alpha = 0.5)
mle.surv <- survfit(Surv(os_months, os_deceased) ~ karnofsky_performance_score,
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
  #create a helping matrix for finding event time
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
  #create data set
  sim.data <- data_frame(surv_months = t) %>%
    mutate(os_status = ifelse(is.na(surv_months), 'LIVING', 'DECEASED'),
           surv_months = ifelse(is.na(surv_months), tau[length(tau)], surv_months),
           id = seq(n), 
           censor_months = rexp(n = n, rate = 1/100))   %>% #censoring rate
    dplyr::mutate(os_status = ifelse(surv_months < censor_months & os_status != 'LIVING',
                                     'DECEASED', 'LIVING'
    ),
    os_months = ifelse(surv_months < censor_months  & os_status != 'LIVING',
                       surv_months, censor_months
    )
    ) %>%   cbind(X) #joint covariates
  return(sim.data)
}
set.seed(342)
test_n = 100
test_tau = c(seq(0, 1200, length.out = test_n))
test_baseline <- exp(-3)*runif(test_n - 1 , 0, 1)
X_test = matrix(c(rnorm(100), sample(c(0,1), 100, replace = T)), ncol=2)
test_beta = c(0.5, 1)
sim_data <-  pem_sim_data( beta = test_beta,
                           X = X_test,
                           tau = test_tau,
                           lambda = test_baseline,
                           n = test_n)
## plot KM curve from simulated data
sim_data <- 
  sim_data %>%
  dplyr::mutate(os_deceased = os_status == 'DECEASED') %>%
  rename("continuos" = "1", "discrete" = "2")

autoplot(survival::survfit(Surv(os_months, os_deceased) ~ 1,
                           data = sim_data
), conf.int = F) + 
  ggtitle('Simulated KM curve')
#------ long data format ----#
#set the tau interval times
tau <- sim_data %>% select(os_months) %>% unlist %>% unique %>% sort()
longdata <- survival::survSplit(Surv(time = os_months, 
                                     event = deceased) ~ . , 
                                cut = tau, data = (sim_data %>%
                                mutate(deceased = os_status == "DECEASED")))

#create time point id
longdata <- longdata %>%
  arrange(id, os_months)  %>%
  group_by(id) %>%
  mutate(t_id = seq(n())) %>%
  ungroup()
t_obs <- sim_data %>% select(os_months) %>% unlist %>% unique %>% sort()
t_dur <- diff(tau)
#----Generate stan data----#
M = length(test_beta)
gen_stan_data <- function(data){
  stan_data <- list(
    N = nrow(data),
    S = length(unique(longdata$id)),
    "T" = dplyr::n_distinct(data$t_id),
    s = array(as.integer(data$id)),
    t_obs = t_obs,
    t_dur = t_dur,
    M=M,
    event = as.integer(data$deceased),
    t = data$t_id,
    x = array(matrix(c(data$continuos, data$discrete), ncol=M), dim=c(nrow(data), M))
  )
}
#---Set initial values---#
gen_inits <-  function(M) {
  function() 
  list(
    beta_bg_raw = rnorm(M),
    tau_s_bg_raw = 0.1*abs(rnorm(1)),
    tau_bg_raw = abs(rnorm(M)),
    c_raw = abs(rnorm(1)),
    r_raw = abs(rnorm(1)),
    baseline = rgamma(n = length(diff(tau)), shape =  mean(diff(tau)) * 0.1, scale = 0.01)
  )
}
#-----Run Stan-------#
nChain <- 1
stanfile <- 'pem_bg.stan'
rstan_options(auto_write = TRUE)
simulated_fit <- stan(stanfile,
                      data = gen_stan_data(longdata),
                      init = gen_inits(M=2),
                      iter = 1000,
                      cores = min(nChain, parallel::detectCores()),
                      seed = 7327,
                      chains = nChain,
                      pars = c("beta_bg", "baseline", "lp__")
                      #control = list(adapt_delta = 0.99)
)
#----Convergence review -----#
print(simulated_fit)
pairs(simulated_fit, pars = c("lp__", "beta_bg"), las = 1)
rstan::traceplot(simulated_fit, 'lp__')
rstan::traceplot(simulated_fit, 'beta_bg')
if (interactive())
  shinystan::launch_shinystan(simulated_fit)
#---Review posterior distribution of beta parameters--#
pp_beta1 <- rstan::extract(simulated_fit,'beta_bg[1]')$beta_bg
pp_beta2 <- rstan::extract(simulated_fit,'beta_bg[2]')$beta_bg
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
pp_beta_bg <- as.data.frame.array(rstan::extract(simulated_fit,pars = 'beta_bg', permuted = TRUE)$beta_bg) 
pp_lambda <- as.data.frame.array(rstan::extract(simulated_fit,pars = 'baseline', permuted = TRUE)$baseline)
# create list
pp_beta_bg <-  split(pp_beta_bg, seq(nrow(pp_beta_bg)))
pp_lambda <-  split(pp_lambda, seq(nrow(pp_lambda)))
pp_newdata <- 
  purrr::pmap(list(pp_beta_bg, pp_lambda),
              function(pp_beta, pp_lambda) {pem_sim_data(lambda = pp_lambda, 
                                           beta = pp_beta,
                                           tau = tau,
                                           n = test_n,
                                           X = X_test)} )
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
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + xlim(c(0,100))

#----Wrap as a function----#
pp_predict_surv <- function(pp_beta, pp_lambda, n, tau, X,
                            level = 0.9, 
                            plot = F, data = NULL,
                            sim_data_fun = pem_sim_data) {
  pp_newdata <- 
    purrr::pmap(list(pp_lambda, pp_beta),
                function(a, b) {pem_sim_data(lambda = a, 
                                             beta = b,
                                             tau = tau,
                                             n = n,
                                             X = X)} )
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
                                             probs = upper_p) ) %>%
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
              os_deceased = os_status == 'DECEASED') )) %>%
        dplyr::mutate(lower = surv,
                      upper = surv, type = 'actual data') )
  pl <- ggplot(ggplot_data,
               aes(x = time, group = type, linetype = type)) + 
    geom_line(aes(y = surv, colour = type)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2)
  
  pl 
}

#----Fitting Model to TCGA Glioblastome data----#
#-----Covariates-----#
#Create dummy vars
clinical_data$karnofsky_performance_score <- as.factor(clinical_data$karnofsky_performance_score)
Xdummies <- caret::dummyVars(Surv(os_months, os_deceased) ~ karnofsky_performance_score,
                             data =  clinical_data %>%
                               mutate(os_deceased = (os_status == "DECEASED")))
X <- tbl_df(predict(Xdummies, newdata =  clinical_data %>%
                      mutate(os_deceased = (os_status == "DECEASED"))))

#Think about imputation
pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(clinical_data %>% select(karnofsky_performance_score), 2, pMiss)
#we can discard the missing observations but is more than 5% thus by the rule of thumb is better to find an imputation strategy
clinical_data <- clinical_data %>%
  filter(!is.na(karnofsky_performance_score))  #clean data
#special boxplot
VIM::marginplot(clinical_data[c(1,2)])
#using imputation by Bayesian logistic regression
tempData <- mice(clinical_data,m=5,maxit=50,meth='pmm',seed=500)
summary(tempData)
tmp <- as.list(X)
lapply(tmp, function(x){
  mice.impute.logreg.boot( 
    y = as.vector(x),
    ry = !is.na(clinical_data$karnofsky_performance_score),
    x = model.matrix(~ prior_glioma + sex + treatment_status +
                       pretreatment_history + os_months + os_status,
                     data = clinical_data)[,-1])
})
#or using bagged trees
X <- X %>% cbind(clinical_data %>% select(prior_glioma, sex,
                                          treatment_status, 
                        pretreatment_history, os_months, os_status))
preProc <- caret::preProcess(X, method = c("bagImpute"))
X <- predict(preProc, X, na.action = na.pass)
X <- X %>% 
  select(dplyr::contains("karnofsky_performance_score")) %>% round()
# or multiple imputation
nb <- missMDA::estim_ncpMCA(as.data.frame(
  clinical_data %>%
  select(karnofsky_performance_score, os_status) %>%
    mutate_all(funs(as.factor))),
                            ncp.max=5) ## Time-consuming, nb = 4

# the bagged trees may be prefered for its computational efficiency

#Near Zero Variance Predictors
nzv <- caret::nearZeroVar(X, saveMetrics= TRUE)
nzv[nzv$nzv,]
#joint covariates lower than 60
X <- X %>% mutate(kpsless_or60 = (karnofsky_performance_score.40 | karnofsky_performance_score.60)) %>%
  rename(kps80 = karnofsky_performance_score.80, kps100 = karnofsky_performance_score.100) %>%
  select(-karnofsky_performance_score.40, -karnofsky_performance_score.60) %>%
  mutate_all(funs(as.integer))
#Double check near Zero Variance Predictors
nzv <- caret::nearZeroVar(X, saveMetrics= TRUE)
nzv[nzv$nzv,]

clinical_data <- clinical_data %>%
  select(sample_id, num_id, os_months, os_status) %>%
  cbind(X)
#--- Update gen stan data function to include covariates. This function will take a formula object as input --- #
gen_stan_data <- function(data, formula = as.formula(~1)) {
  if(!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  #set the tau interval times
  tau <- data %>% select(os_months) %>% unlist %>% unique %>% sort()
  if(tau[1] != 0){
    tau <- c(0, tau)
  }
  longdata <- survival::survSplit(Surv(time = os_months,
                                       event = deceased) ~ . , 
                                  cut = tau, data = (data %>%
                                  mutate(deceased = os_status == "DECEASED")))
  #create time point id
  longdata <- longdata %>%
    group_by(sample_id) %>%
    mutate(t_id = seq(n())) %>%
    ungroup()
  t_obs <- data %>% select(os_months) %>% unlist %>% unique %>% sort()
  t_dur <- diff(tau)
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
    S = length(unique(longdata$num_id)),
    "T" = dplyr::n_distinct(longdata$t_id),
    s = as.integer(longdata$num_id),
    t_obs = t_obs,
    t_dur = t_dur,
    M = M_bg,
    event = as.integer(longdata$deceased),
    t = longdata$t_id,
    x = X_bg
  )
}
#---Set initial values---#
gen_inits <-  function(M) {
  function() 
    list(
      beta_bg_raw = array(rnorm(M), dim =M),
      tau_s_bg_raw = 0.1*abs(rnorm(1)),
      tau_bg_raw = array(abs(rnorm(M)), dim = M),
      c_raw = abs(rnorm(1)),
      r_raw = abs(rnorm(1)),
      baseline = rgamma(n = length(diff(tau)), shape =  mean(diff(tau)) * 0.1, scale = 0.01)
    )
}

#-----Run Stan-------#
nChain <- 4
stanfile <- 'pem_bg.stan'
rstan_options(auto_write = TRUE)
pem_fit <- stan(stanfile,
                      data = gen_stan_data(clinical_data,
                                           '~ kpsless_or60 + 
                                           kps80 + kps100'),
                      init = gen_inits(M=3),
                      iter = 250,
                      cores = min(nChain, parallel::detectCores()),
                      seed = 7327,
                #control = list(adapt_delta = 0.99, max_treedepth = 15
                      chains = nChain,
                      pars = c("beta_bg", "baseline", "lp__"))
###########
#----Convergence review -----#
print(pem_fit, pars = c("baseline"))
pairs(pem_fit, pars = c("lp__", "beta_bg"), las = 1)
rstan::traceplot(pem_fit, 'lp__')
rstan::traceplot(pem_fit, 'beta_bg')
if (interactive())
  shinystan::launch_shinystan(simulated_fit)

###########
#-----Posteriro predictive check-----#
pp_baseline <- rstan::extract(pem_fit,'baseline')$baseline
pp_beta <- rstan::extract(pem_fit, 'beta_bg')$beta_bg
pp_baseline <-  split(pp_baseline, seq(nrow(pp_baseline)))
pp_beta <-  split(pp_beta, seq(nrow(pp_beta)))

#Create dummy vars
Xdummies <- dummyVars(Surv(os_months, os_deceased) ~ age +
                        g.cimp_methylation + idh1_status +
                        mgmt_status, data =  glio_clin_dat %>%
                        mutate(os_deceased = (os_status == "DECEASED")))
X <- tbl_df(predict(Xdummies, newdata =  glio_clin_dat %>%
                      mutate(os_deceased = (os_status == "DECEASED"))))
names(X)<- stringr::str_replace_all(names(X), "-", "")
names(X) <- tolower(names(X))

tau <- data %>% select(os_months) %>% unlist %>% unique %>% sort()

pl <- pp_predict_surv(pp_beta = pp_beta,
                      pp_lambda = pp_baseline,
                      n = nrow(clinical_data),
                      tau = clinical_data %>% select(os_months) %>% unlist %>% unique %>% sort(),
                      X = X,
                      level = 0.9, 
                      plot = T, 
                      data = clinical_data) 
pl + 
  xlim(NA, 250) +
  ggtitle('Posterior predictive checks \nfit to GBC 2008 historical cohort; showing 90% CI')
 

rcorr.cens(age, Surv(d.time, death))
r <- rcorrcens(Surv(d.time, death) ~ age + bp)
r
plot(r)


