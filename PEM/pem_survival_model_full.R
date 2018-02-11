#---- PEM Survival Model --- #
library(purrr)
suppressMessages(library(tidyverse))
library(survival)
library(rstan)
library(assertthat)
library(corrplot)
library(cgdsr)
suppressMessages(library(dplyr))

library(ggplot2)
require(ggfortify)
theme_set(theme_bw())
###############################################
#Data obtantion
#------Obtain data by the cgdsr package from MSKCC CBioPortal ----# 

mycgds = CGDS("http://www.cbioportal.org/public-portal/")

study_list = getCancerStudies(mycgds)

id_sutdy = getCancerStudies(mycgds)[55,1]
case_list = getCaseLists(mycgds, id_sutdy)[2,1]
clinical_data <-  tbl_df(getClinicalData(mycgds, case_list)) 
clinical_data <- clinical_data %>% tibble::rownames_to_column("sample") 

#inspect dataframe
glimpse(clinical_data)

#separate two cohorts
id_2008_sutdy = getCancerStudies(mycgds)[56,1] #cohort 2008
case_list_2008 = getCaseLists(mycgds, id_2008_sutdy)[2,1]
clinical_data_2008 <-  tbl_df(getClinicalData(mycgds, case_list_2008))
clinical_data_2008 <- clinical_data_2008 %>% tibble::rownames_to_column("sample")
# 
clinical_data <- clinical_data %>% 
   filter(sample %in% setdiff(clinical_data$sample, clinical_data_2008$sample))
 
####################################################################
#Data Cleaning

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

clinical_data %>% 
  filter(is.na(os_status) | os_status == "") %>%
  select(os_months, os_status) %>%
  glimpse
  
clinical_data %>% 
  filter(is.na(os_status) | os_status == "" |os_months < 0 | is.na(os_months)) %>%
  select(os_months, os_status) %>%
  glimpse

#For the moment we will remove these observations from the analysis
clin_data <- tbl_df(clinical_data)
clinical_data <- 
  clin_data %>% 
  filter(!is.na(os_status) & os_status != "" )

assertthat::assert_that(nrow(clinical_data) == (nrow(clin_data) - 44))
remove(clin_data)

clinical_data %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

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
require(ggfortify)
ggplot2::autoplot(mle.surv, conf.int = F) +
  ggtitle('KM survival for GGM 2013 Cohort')


#------------ Prepare for fitting model--------------------------#

#Center continuos covariates

ggplot(clinical_data, aes(x = age))+
  geom_density()+
  geom_vline(xintercept = median(clinical_data$age))

centered <- function(x){
  x_centered <- x - mean (x)
  return(x_centered)
}

clinical_data <- clinical_data %>%
  dplyr::mutate( age_centered = centered(age)) 

#Impute or delete missing Covariate values

clinical_data %>%
  select(mgmt_status) %>%
  table(exclude = NULL)

clinical_data %>%
  select(idh1_status, g.cimp_methylation) %>%
  table(exclude = NULL)


#remove nas
clinical_data <- clinical_data %>%
  filter(!is.na(mgmt_status))

clinical_data %>%
  select(mgmt_status, g.cimp_methylation) %>%
  table(exclude = NULL)


#generate data we need A long data formating


gen_stan_data <- function(data, formula = as.formula(~1)) {
  

  if(!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
   #add sample id
  data <- data %>%
    mutate(s = seq(n()))
  
  #obtain vector of unique event times
  t_obs <- data %>% 
    filter(os_status == "DECEASED") %>%
    select(os_months) %>% unique %>% ungroup %>% arrange(os_months) %>% unlist
  
  t_dur <- diff(c(0,t_obs))
  
  #create long data
  longdata <- survival::survSplit(Surv(time = os_months, event = deceased) ~ . , 
                                  cut = t_obs, data = (data %>%
                                                         filter(!is.na(os_status) & os_status != '') %>%
                                                         filter(os_months > 0 & !is.na(os_months))  %>%
                                                         mutate(deceased = os_status == "DECEASED")))
  #create time point id
  longdata <- longdata %>%
    group_by(sample) %>%
    mutate(t = seq(n())) 
  
  #covariate matrix
  x <- longdata %>%
    model.matrix(formula, data = .)
  M <- ncol(x)
  
  if (M > 1){
    if("(Intercept)" %in% colnames(x))
      x <- array(x[,-1], dim = c(nrow(longdata), M -1))
    M <- ncol(x)
  }
  
  stan_data <- list(
    N = nrow(longdata),
    S = dplyr::n_distinct(longdata$sample),
    "T" = length(t_obs),
    M = M, 
    s = array(as.numeric(longdata$s)),
    t = array(longdata$t),
    event = array(longdata$deceased),
    t_obs = array(t_obs),
    t_dur = array(t_dur),
    x = array(x, dim = c(nrow(longdata), M))
  )
}

#---------- inits function ----------#
stanfile <- "pem.stan"
biostan::print_stan_file(stanfile, section = 'parameters')

#t length of t_obs

t_obs <- clinical_data %>% 
  filter(os_status == "DECEASED") %>%
  select(os_months) %>% unique %>% ungroup %>% arrange(os_months) %>% unlist
t <- length(t_obs)

M <- 3 #number of covariates

gen_inits <- function(M = M, t = t){
  function()
    list(
      beta = array(rnorm(M), dim = c(M)),
      baseline_sigma = abs(rnorm(1)),
      log_baseline_mu = rnorm(1)
    )
}



###------ Run Stan --------##
nChains <- 4

# #open stan file
if (interactive())
  file.edit(stanfile)

stan(stanfile,
     data = gen_stan_data(clinical_data, '~ age_centered + mgmt_status + g.cimp_methylation'),
     iter = 5,
     cores = min(1, parallel::detectCores()),
     init = gen_inits,
     chains = 1)

fit_pem_age <-  stan(stanfile,
                     data = gen_stan_data(clinical_data, '~ age_centered + mgmt_status + g.cimp_methylation'),
                     iter = 1000,
                     init = gen_inits(M = M, t = t),
                     cores = min(nChains, parallel::detectCores()),
                     chains = nChains
)

fit_pem_age <-  stan(stanfile,
                     data = gen_stan_data(clinical_data, '~ age_centered + mgmt_status + g.cimp_methylation'),
                     iter = 2000,
                     control = list(max_treedepth =15),
                     cores = min(nChains, parallel::detectCores()),
                     chains = nChains
)


fit_pem_age <-  stan(stanfile,
                            data = gen_stan_data(clinical_data, '~ age_centered + mgmt_status + g.cimp_methylation'),
                            iter = 2000,
                            init = gen_inits(M = M, t = t),
                            control = list(stepsize = 0.01, adapt_delta = 0.99),
                            cores = min(nChains, parallel::detectCores()),
                            chains = nChains
)


##--Review model convergence--#

#Fit object
print(fit_pem_age)  #Rhat are close to 1?

#Traceplots
rstan::traceplot(fit_pem_age, 'lp__')
rstan::traceplot(fit_pem_age, 'beta')

if(interactive())
  shinystan::launch_shinystan(fit_pem_age)        #Launch shiny stan. There are some divergent 

##--- Review parameter estimates --#

pp_beta <- rstan::extract(fit_pem_age, 'beta', permuted = TRUE)$beta
ggplot(data = data.frame(beta = unlist(pp_beta)),
       aes(x = beta))+
  geom_density()


#How likely is the coefficient be greater than 0
mean(pp_beta >= 0)




#------------------ Posterior predictive check ---------------------#
#Simulate time to event data
#weibull_sim_data function takes two parameters (alpha and mu) as inputs and a desired sample size (n). 
pem_sim_data <- function(alpha, mu, n) {
  
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

##pp function

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

