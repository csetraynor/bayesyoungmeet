#---- Based on biostan package by Dr Jacqueline Buros --- #
library(purrr)
suppressMessages(library(tidyverse))
library(survival)
library(rstan)
library(assertthat)
library(corrplot)
library(cgdsr)

library(ggplot2)
theme_set(theme_bw())
###############################################
#Data obtantion
#------Obtain data by the cgdsr package from MSKCC CBioPortal ----# 

mycgds = CGDS("http://www.cbioportal.org/public-portal/")

study_list = getCancerStudies(mycgds)

id_sutdy = getCancerStudies(mycgds)[55,1]
case_list = getCaseLists(mycgds, id_sutdy)[2,1]
clinical_data <-  tbl_df(getClinicalData(mycgds, case_list)) 

#inspect dataframe
glimpse(clinical_data)


####################################################################
#Data Cleaning

clinical_data <- clinical_data %>% tibble::rownames_to_column("sample") 
#convert to lower case
names(clinical_data) <- tolower(names(clinical_data)) 

#convert missig values
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

#inspect resulting dataframe
glimpse(clinical_data)
clinical_data %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

######################################################################
#######--------------  Data Exploration  ----------------#################
#---------   Considering overall survival   ----------------------#

#filter unknown or negative survival times (os_monts < 0)
clinical_data %>%
  filter(is.na(os_status) | os_status == '') %>%
  filter(os_months <= 0 | is.na(os_months)) %>%
  select(os_status, os_months) %>%
  dplyr::glimpse()
  

clinical_data %>%
  filter(is.na(os_status) | os_status == '') %>%
  select(os_status, os_months) %>%
  str() 

#for now this observation will be remove from the analysis

clinical_data <- clinical_data %>%
  filter(!is.na(os_status) & os_status != '') %>%
  filter(os_months >= 0 & !is.na(os_months))

#Check 44 fewer observations than original
assertthat::assert_that(nrow(clinical_data) == nrow(glioblastome_2013_clinical_data) - 44)

########## Distribution of event times  ######################

clinical_data %>%
  ggplot(aes(x = os_months,
             group = os_status,
             colour = os_status,
             fill = os_status)) +
  geom_density(alpha = 0.5)

#KM curve
mle.surv <- survfit(Surv(os_months, os_deceased) ~ 1,
                    data = clinical_data %>%
                      mutate(os_deceased = (os_status == "DECEASED")))
ggplot2::autoplot(mle.surv, conf.int = F) +
  ggtitle('KM survival for GGM Cohort')

##############################################################
#######------ Considering disease free survival ------ #######

#filter unknown or negative survival times (os_monts < 0)
clinical_data %>%
  filter(is.na(dfs_status) | dfs_status == '') %>%
  filter(dfs_months <= 0 | is.na(dfs_months)) %>%
  dplyr::glimpse()

#for now this observation will be remove from the analysis

clinical_data <- clinical_data %>%
  filter(!is.na(dfs_status) | dfs_status != '') %>%
  filter(dfs_months > 0  | !is.na(dfs_months)) 

#Check 43 fewer observations than original
assertthat::assert_that(nrow(clinical_data) == nrow(glioblastome_2013_clinical_data) - 78)

##################################################################
##########------  Parametric Survival Model  --- #####################

#########----- generate dataset for Null Model ------########
gen_stan_data <- function(data) {
  observed_data <- data %>%
    dplyr::filter(os_status == 'DECEASED')
  
  censored_data <- data %>%
    dplyr::filter(os_status != 'DECEASED')
  
  stan_data <- list(
    Nobs = nrow(observed_data),
    Ncen = nrow(censored_data),
    yobs = observed_data$os_months,
    ycen = censored_data$os_months
  )
}

#---- Setting intial values for Weibull Distribution Stan http://mc-stan.org/users/documentation/ --#
gen_inits <- function() {
  list(
    alpha_raw = 0.01*rnorm(1),
    mu = rnorm(1)
  )
}

#---- Run Stan ----#

stanfile <- "ClinicoGenomicBayesianModels/ClinicalGliobastomeParametricWithWeibull.stan"
#open stan file
if (interactive())
  file.edit(stanfile)

weibull_null_model <-  stan(stanfile,
                            data = gen_stan_data(clinical_data),
                            chains = 4,
                            iter = 1000,
                            init = gen_inits
)

####################### Checking convergence ###################

print(weibull_null_model) #(Check Rhat close to 1)

rstan::traceplot(weibull_null_model, 'lp__') #Review traceplot for log-posterior

rstan::traceplot(weibull_null_model, c('alpha','mu'), ncol = 1)    #Review traceplot for parameters of interest

if(interactive())
  shinystan::launch_shinystan(weibull_null_model)        #Launch shiny stan


################# Posterior predictive check ###################
#Simulate time to event data
#weibull_sim_data function takes two parameters (alpha and mu) as inputs and a desired sample size (n). 
weibull_sim_data <- function(alpha, mu, n) {
  
  data <- data.frame(surv_months = rweibull(n = n, alpha, exp(-(mu)/alpha)),
                     censor_months = rexp(n = n, rate = 1/100),
                     stringsAsFactors = F
  ) %>%
    dplyr::mutate(os_status = ifelse(surv_months < censor_months,
                                     'DECEASED', 'LIVING'
    ),
    os_months = ifelse(surv_months < censor_months,
                       surv_months, censor_months
    )
    )
  
  return(data)
}

################# Posterior predictive check

pp_predict_surv <- function(pp_alpha, pp_mu, n,
                            level = 0.9,
                            plot = F, data = NULL,
                            sim_data_fun = weibull_sim_data
) {
  pp_newdata <- 
    purrr::map2(.x = pp_alpha,
                .y = pp_mu,
                .f = ~ sim_data_fun(alpha = .x, mu = .y, n = n)
    )
  
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
  
  if (!is.null(data)){
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
      )}
  
  pl <- ggplot(ggplot_data,
               aes(x = time, group = type, linetype = type)) + 
    geom_line(aes(y = surv, colour = type)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2)
  
  pl 
}

pl <- pp_predict_surv(pp_alpha = extract(weibull_null_model,'alpha')$alpha,
                      pp_mu = extract(weibull_null_model,'mu')$mu,
                      n = nrow(clinical_data),
                      data = clinical_data,
                      plot = T
) 
pl + 
  xlim(NA, 250) +
  ggtitle('Posterior predictive checks for NULL weibull model\nfit to GGC data; showing 90% CI')

# Proportion of times is the observed event rate within the 90 % Confidence Interval

## summarize 90% CI of predicted event rate for each interval
pp_agg <- pp_predict_surv(pp_alpha = extract(weibull_null_model,'alpha')$alpha,
                          pp_mu = extract(weibull_null_model,'mu')$mu,
                          n = nrow(clinical_data)
)


## summarize observed data into same time_groups
act_agg <- 
  survival::survfit(Surv(os_months, I(os_status == 'DECEASED')) ~ 1,
                    data = clinical_data
  ) %>%
  fortify() %>%
  dplyr::mutate(time_group = floor(time)) %>%
  dplyr::group_by(time_group) %>%
  dplyr::summarise(observed_surv = mean(surv)) %>%
  dplyr::ungroup()

## compute proportion of observed values within 90% ci
act_agg %>%
  dplyr::inner_join(pp_agg, by = 'time_group') %>%
  dplyr::mutate(within_interval = ifelse(observed_surv >= surv_lower & observed_surv <= surv_upper,
                                         1, 0),
                time_set = cut(time_group, breaks = c(0,100))
  ) %>%
  dplyr::group_by(time_set) %>%
  dplyr::summarize(mean(within_interval))


############# ---- Add Covariates ---- #######


#Correlation overview
stan_file <- system.file('stan', 'weibull_survival_model.stan', package =  'biostan')

## open stan file to review contents 
if (interactive())
  file.edit(stan_file)

biostan::print_stan_file(stan_file, section = 'parameters')


##########Select Clinical Variables

#Dummy Variables for Treatment

clinical_data <- clinical_data %>%
    tibble::rownames_to_column(var = "index")  %>%
  mutate(
    tmz_therapy = index %in% (dplyr::starts_with("TMZ", vars = therapy_class)),
    unspecified_therapy = index %in% (dplyr::starts_with("Unspecified", vars = therapy_class)),
    nonstandard_therapy = index %in% (dplyr::starts_with("Nonstandard", vars = therapy_class)),
    alkylating_therapy = index %in% (dplyr::starts_with("Alkylating", vars = therapy_class)),
    comb_tmz_rad_therapy = I(therapy_class == "TMZ Chemoradiation, TMZ Chemo")
  )


#--- Update gen stan data function to include covariates ---#

#--- This function will take a formula object as input --- #
gen_stan_data <- function(data, formula = as.formula(~1)) {
  if(!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  observed_data <- data %>%
    dplyr::filter(os_status == "DECEASED")
  
  censored_data <- data %>%
    dplyr::filter(os_status != "DECEASED")
  
  Xobs_bg <- observed_data %>%
    model.matrix(formula, data = .)
  
  Xcen_bg <- censored_data %>%
    model.matrix(formula, data = .)
  
  assertthat::assert_that(ncol(Xcen_bg) == ncol(Xobs_bg))
  M_bg <- ncol(Xcen_bg)
  
  if (M_bg > 1){
    if("(Intercept)" %in% colnames(Xobs_bg))
      Xobs_bg <- array(Xobs_bg[,-1], dim = c(nrow(observed_data), M_bg - 1))
    if("(Intercept)" %in% colnames(Xcen_bg))
      Xcen_bg <- array(Xcen_bg[,-1], dim = c(nrow(censored_data), M_bg - 1))
    assertthat::assert_that(ncol(Xcen_bg) == ncol(Xobs_bg))
    M_bg <- ncol(Xcen_bg)
  }
  
  stan_data <- list(
    Nobs = nrow(observed_data),
    Ncen = nrow(censored_data),
    yobs = observed_data$os_months,
    ycen = censored_data$os_months,
    M_bg = M_bg,
    Xcen_bg = array(Xcen_bg, dim = c(nrow(censored_data), M_bg)),
    Xobs_bg = array(Xobs_bg, dim = c(nrow(observed_data), M_bg))
  )
}

stan_input_data <- gen_stan_data(clinical_data, '~ comb_tmz_rad_therapy + I(sex == "Male") +
                                            alkylating_therapy + nonstandard_therapy + age') 

#---------- Update inits function ----------#
biostan::print_stan_file(stan_file, section = 'parameters')

gen_inits2 <- function(M_bg){
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1),
      tau_s_bg_raw = 0.1*abs(rnorm(1)),
      tau_bg_raw = array(abs(rnorm(M_bg)), dim = c(M_bg)),
      beta_bg_raw = array(rnorm(M_bg), dim = c(M_bg))
    )
}

###------ Run Stan --------##
nChains <- 4

testfit <- rstan::stan(stan_file,
                       data = gen_stan_data(clin_data, '~ I(sex == "Male")'),
                       init = gen_inits,
                       iter = 4,
                       chains = 1
)

fullfit <- rstan::stan(stan_file,
                       data = gen_stan_data(clinical_data, '~ comb_tmz_rad_therapy + I(sex == "Male") +
                                            alkylating_therapy + nonstandard_therapy + age'),
                       init = gen_inits2(M_bg = 1),
                       iter = 5000, 
                       control = list(stepsize = 0.01, adapt_delta = 0.99),
                       cores = min(nChains, parallel::detectCores()),
                       chains = nChains)

##--Review model convergence--#

#Fit object
print(fullfit)  #Rhat are close to 1?

#Traceplots
rstan::traceplot(fullfit, 'lp__')
rstan::traceplot(fullfit, c('alpha', 'mu'), ncol = 1)
rstan::traceplot(fullfit, 'beta_bg')

if(interactive())
  shinystan::launch_shinystan(weibull_null_model)        #Launch shiny stan. There are some divergent transitions but overall the sampling has gone well

##--- Review parameter estimates --#

pp_beta_bg <- rstan::extract(fullfit, 'beta_bg')$beta_bg
ggplot(data = data.frame(beta_bg = unlist(pp_beta_bg)),
       aes(x = beta_bg)) + 
  geom_density()

#How likely is the coefficient be greater than 0
mean(pp_beta_bg >= 0)

#How well does this model fit the data

#among not comb treated
pp_predict_surv(pp_alpha = rstan::extract(fullfit, 'alpha')$alpha,
                pp_mu = rstan::extract(fullfit, 'mu')$mu,
                n = nrow(clinical_data %>% dplyr::filter(comb_tmz_rad_therapy == 0)),
                data = clinical_data %>% dplyr::filter(comb_tmz_rad_therapy == 0),
                plot = TRUE)

#among comb treated
pp_predict_surv(pp_alpha = rstan::extract(fullfit, 'alpha')$alpha,
                pp_mu = rstan::extract(fullfit, 'mu')$mu,
                n = nrow(clinical_data %>% dplyr::filter(comb_tmz_rad_therapy == 1)),
                data = clinical_data %>% dplyr::filter(comb_tmz_rad_therapy == 1),
                plot = TRUE)

