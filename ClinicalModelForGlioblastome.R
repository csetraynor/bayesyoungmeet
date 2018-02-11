#----Based on biostan package and the excellent explanation by Dr Jacqueline Buros ---#

library(purrr)
library(dplyr)
library(ggplot2)
library(survival)
library(rstan)
library(scales)
###############################################
#Data obtantion
#get data from with MSKCC package 
require(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

glioblastome_2013_id_sutdy = getCancerStudies(mycgds)[55,1]
glioblastome_2013_case_list = getCaseLists(mycgds, glioblastome_2013_id_sutdy)[2,1]
glioblastome_2013_clinical_data <-  getClinicalData(mycgds, glioblastome_2013_case_list)

#inspect dataframe
str(glioblastome_2013_clinical_data, no.list = T, vec.len = 2)

####################################################################
#Data Cleaning

#convert to lower case
names(glioblastome_2013_clinical_data) <- tolower(names(glioblastome_2013_clinical_data)) 

#convert missig values
convert_blank_to_na <- function(x) {
  if(!purrr::is_character(x)){
    warning('input variate not character - return original')
    return(x)
  } else {
    ifelse(x == '', NA, x)
  }
}
glio_clin_dat <- glioblastome_2013_clinical_data %>%
  dplyr::mutate_all(funs(convert_blank_to_na))

#inspect resulting dataframe
str(glio_clin_dat)

######################################################################
#Data Exploration
#Considering overall survival#

glio_clin_dat %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

#filter unknown or negative survival times (os_monts < 0)

glio_clin_dat %>%
  filter(!is.na(os_status) & os_status != '') %>%
  filter(os_months < 0 | is.na(os_months)) %>%
  select(os_status, os_months) %>%
  head()

glio_clin_dat %>%
  filter(is.na(os_status) | os_status == '') %>%
  select(os_status, os_months) %>%
  str() 

#for now this observation will be remove from the analysis

glio_clin_dat <- glio_clin_dat %>%
  filter(!is.na(os_status) & os_status != '') %>%
  filter(os_months >= 0 & !is.na(os_months))
  
#Check 44 fewer obsrvations than original
assertthat::assert_that(nrow(glio_clin_dat) == nrow(glioblastome_2013_clinical_data) - 44)


########## Distribution of event times  ######################

glio_clin_dat %>%
  ggplot(aes(x = os_months,
             group = os_status,
             colour = os_status,
             fill = os_status)) +
  geom_density(alpha = 0.5)

#KM curve

mle.surv <- survfit(Surv(os_months, os_deceased) ~ 1,
                    data = glio_clin_dat %>%
                      mutate(os_deceased = (os_status == "DECEASED")))
ggplot2::autoplot(mle.surv, conf.int = F) +
  ggtitle('KM survival for GGM Cohort')



############ Parametric Survival Model #####################


observed_data <- glio_clin_dat %>%
  filter(os_status == "DECEASED")

censored_data <- glio_clin_dat %>%
  filter(os_status != "DECEASED")

stan_data <- list(
  Nobs = nrow(observed_data),
  Ncen = nrow(censored_data),
  yobs = observed_data$os_months,
  yceb = censored_data$os_months
)
rm(censored_data)
rm(observed_data)

str(stan_data)

#Wraped in a function

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

######### Setting intial values

gen_inits <- function() {
  list(
    alpha_raw = 0.01*rnorm(1),
    mu = rnorm(1)
  )
}

########################## Stan run ###########################

stanfile <- "ClinicoGenomicBayesianModels/ClinicalGliobastomeParametricWithWeibull.stan"
#open stan file
if (interactive())
  file.edit(stanfile)

weibull_null_model <-  stan(stanfile,
              data = gen_stan_data(glio_clin_dat),
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



######################### Posterior predicitive checks ###################################
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
#Censoring is "arbitrarily" rexp() , censoring is assumed to be noninformative.

######## Simulating data for each posterior draw #
test_n <- nrow(glio_clin_dat)
pp_newdata <- purrr::map2(.x = pp_alpha,
                          .y = pp_mu,
                          .f = ~weibull_sim_data(alpha = .x,
                                                 mu = .y,
                                                 n = test_n))

###### Plot time to event in the posterior draws compare to actual time in dataset
ggplot(pp_newdata %>%
         bind_rows() %>%
         mutate(type = 'posterior predicted values') %>%
         bind_rows(glio_clin_dat %>% mutate(type = 'actual data'))
       , aes(x = os_months, group = os_status, colour = os_status, fill = os_status))+ 
  geom_density(alpha = 0.5) +
  facet_wrap(~type, ncol = 1)

#### summarise posterior predictive draws

## cumulative survival rate at each draw from the posterior

pp_survdata <- 
  pp_newdata %>%
  map(~ mutate(., os_deceased = os_status == 'DECEASED')) %>%
  map(~ survfit(Surv(os_months, os_deceased) ~ 1, data = .)) %>%
  map(fortify)
  
## summarise cum survival for each unit time (month), summarised at 95% confidence interval
pp_survdata_agg <- 
  pp_survdata %>%
  map(~mutate(., time_group = floor(time))) %>%
  bind_rows() %>%
  group_by(time_group) %>%
  summarize(surv_mean = mean(surv),
            surv_p50 = median(surv),
            surv_lower = quantile(surv, probs = 0.025),
            surv_upper = quantile(surv, probs = 0.975)) %>%
  ungroup()

## km
kmcurve_data <-   fortify(
    survfit(
      Surv(os_months, os_deceased) ~ 1,
      data = glio_clin_dat %>%
        mutate(os_deceased = os_status == 'DECEASED')
    )) %>%
  mutate(lower =  surv,
         upper = surv)

ggplot(pp_survdata_agg %>%
         mutate(type = 'posterior predicted values') %>%
         rename(surv = surv_p50, lower = surv_lower, upper = surv_upper, time = time_group)
       %>%
         bind_rows(kmcurve_data %>% mutate(type = 'actual data')),
       aes(x = time, group = type, linetype = type)) +
  geom_line(aes(y = surv, colour = type)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  xlim(c(0, 200))


############# Save as a function  ######################

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


################# Posterior predictive check ###################
pl <- pp_predict_surv(pp_alpha = extract(weibull_null_model,'alpha')$alpha,
                      pp_mu = extract(weibull_null_model,'mu')$mu,
                      n = nrow(glio_clin_dat),
                      data = glio_clin_dat,
                      plot = T
) 
pl + 
  xlim(NA, 250) +
  ggtitle('Posterior predictive checks for NULL weibull model\nfit to GGC data; showing 90% CI')

# Proportion of times is the observed event rate within the 90 % Confidence Interval

## summarize 90% CI of predicted event rate for each interval
pp_agg <- pp_predict_surv(pp_alpha = extract(weibull_null_model,'alpha')$alpha,
                          pp_mu = extract(weibull_null_model,'mu')$mu,
                          n = nrow(glio_clin_dat)
)


## summarize observed data into same time_groups
act_agg <- 
  survival::survfit(Surv(os_months, I(os_status == 'DECEASED')) ~ 1,
                    data = glio_clin_dat
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
####### Parame