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

#Simulate data for arbitrary imput values
test_alpha <- 0.8
test_mu <- -3

## sample size from TCGA blca data
test_n <- 500

## test these inputs for arbitrary values of alpha & mu
simulated_data <- 
  weibull_sim_data(alpha = test_alpha,
                   mu = test_mu,
                   n = test_n
  ) 
head(simulated_data)


## plot KM curve from simulated data
require(survival)

simulated_data <- 
  simulated_data %>%
  dplyr::mutate(os_deceased = os_status == 'DECEASED')

ggfortify::autoplot(survfit(Surv(os_months, os_deceased) ~ 1,
                            data = simulated_data
), conf.int = F) + 
  ggtitle('Simulated KM curve')

########################## Stan run ###########################
library(rstan)

stanfile <- "ClinicoGenomicBayesianModels/ClinicalGliobastomeParametricWithWeibull.stan"
#open stan file
if (interactive())
  file.edit(stanfile)

weibull_null_model <-  stan(stanfile,
                            data = gen_stan_data(weibull_sim_data),
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


#Review posterior distributions of parameters

pp_alpha <- rstan::extract(weibull_null_model, 'alpha')$alpha
pp_mu <- rstan::extract(weibull_null_model, 'mu')$mu


ggplot(data.frame(alpha = pp_alpha, mu = pp_mu)) + 
  geom_density(aes(x = alpha)) + 
  geom_vline(aes(xintercept = test_alpha), colour = 'red') +
  ggtitle('Posterior distribution of alpha\nshowing true value in red')

ggplot(data.frame(alpha = pp_alpha, mu = pp_mu)) + 
  geom_density(aes(x = mu)) + 
  geom_vline(aes(xintercept = test_mu), colour = 'red') +
  ggtitle('Posterior distribution of mu\nshowing true value in red')

#Degree of correlation between alpha and mu
ggplot(data.frame(alpha = pp_alpha, mu = pp_mu)) + 
  geom_density2d(aes(x = alpha, y = mu)) +
  geom_point(aes(x = test_alpha, y = test_mu), colour = 'red', size = 2) +
  ggtitle('Posterior distributions of mu and alpha\nshowing true parameter values in red')

#