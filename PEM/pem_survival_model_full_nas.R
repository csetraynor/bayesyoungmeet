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
# id_2008_sutdy = getCancerStudies(mycgds)[56,1] #cohort 2008
# case_list_2008 = getCaseLists(mycgds, id_2008_sutdy)[2,1]
# clinical_data_2008 <-  tbl_df(getClinicalData(mycgds, case_list_2008))
# clinical_data_2008 <- clinical_data_2008 %>% tibble::rownames_to_column("sample")
# clinical_data_2008 <- clinical_data %>%
#   filter(sample %in% clinical_data_2008$sample)
# clinical_data_2013 <- clinical_data %>%
#   filter(!
#       (sample %in% clinical_data_2008$sample))

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
  ggtitle('KM survival for GGM 2008 Cohort')

clinical_data %>%
  filter(!is.na( idh1_status)) %>%
  ggplot(aes(x = os_months, group = os_status, colour = os_status, fill = os_status))+
  geom_density(alpha=0.5)


#Impute variables Clinical Variables of Importance are idh1_status, mgmt_status, g.cimp_methylation and age
glimpse(clinical_data)
clinical_data %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

clinical_data %>%
  filter(is.na(g.cimp_methylation))%>%
  ggplot(aes(x = os_months,
             group = os_status,
             colour = os_status,
             fill = os_status)) +
  geom_density(alpha = 0.5)







#generate data we need A long data formating





gen_stan_data <- function(data, formula = as.formula(~1)) {
  

  if(!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
   #add sample id
  data <- data %>%
    mutate(s = seq(n()))
  
  #obtain vector of unique event times
  times <- data %>% 
    filter(os_status == "DECEASED") %>% select(os_months) %>% unique %>% ungroup %>% arrange(os_months) %>% unlist
  
  #create long data
  longdata <- survival::survSplit(Surv(time = os_months, event = deceased) ~ . , 
                                  cut = times, data = (data %>%
                                                         filter(!is.na(os_status) & os_status != '') %>%
                                                         filter(os_months > 0 & !is.na(os_months))  %>%
                                                         mutate(deceased = os_status == "DECEASED")))
  #create time point id
  longdata <- longdata %>%
    group_by(sample) %>%
    mutate(t = seq(n())) 
  
  #calculate duration
  longdata <- longdata %>%
    group_by(sample) %>%
    mutate(t_dur = os_months - tstart) 
  
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
    "T" = length(times),
    M = M, 
    s = array(as.numeric(longdata$s)),
    t = array(longdata$t),
    event = array(longdata$deceased),
    t_obs = array(longdata$os_months),
    t_dur = array(longdata$t_dur),
    x = array(x, dim = c(nrow(longdata), M))
  )
}

#---------- inits function ----------#
stanfile <- "pem_survival_model.stan"
biostan::print_stan_file(stanfile, section = 'parameters')

gen_inits <- function(M){
  function()
    list(
      beta = array(rnorm(M), dim = c(M))
    )
}



###------ Run Stan --------##
nChains <- 4

# #open stan file
if (interactive())
  file.edit(stanfile)

stan(stanfile,
     data = gen_stan_data(clinical_data, '~ age'),
     iter = 5,
     init = gen_inits(M = 1),
     cores = min(1, parallel::detectCores()),
     chains = 1)

fit_pem_age <-  stan(stanfile,
                            data = gen_stan_data(clinical_data, '~ age'),
                            iter = 5000,
                            init = gen_inits(M = 1),
                            control = list(stepsize = 0.01, adapt_delta = 0.99),
                            cores = min(nChains, parallel::detectCores()),
                            chains = nChains
)


##--Review model convergence--#

#Fit object
print(pem_null_model)  #Rhat are close to 1?

#Traceplots
rstan::traceplot(pem_null_model, 'lp__')
rstan::traceplot(pem_null_model, 'baseline')

if(interactive())
  shinystan::launch_shinystan(pem_null_model)        #Launch shiny stan. There are some divergent 
