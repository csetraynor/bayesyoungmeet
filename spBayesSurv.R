library("coda")
library("survival")
library("spBayesSurv")
library("fields")
library("BayesX")
library("R2BayesX")
data("LeukSurv")
attach(LeukSurv)
d <- LeukSurv[order(district),]
n <- nrow(d); detach(LeukSurv)
head(d)

#PH Model
set.seed(1)
mcmc <- list(nburn = 5000, nsave = 2000, nskip = 4, ndisplay = 1000);
prior <- list(M = 20, nknots = 100, nblock = 1043);
ptm <- proc.time()

sim_data <- 
  sim_data %>%
  dplyr::mutate(os_deceased = os_status == 'DECEASED')


res <-survregbayes(formula = Surv(os_months, os_deceased) ~ continuos + discrete,data = sim_data, survmodel = "PH", mcmc = mcmc, prior = prior)

sfit <- summary(res); sfit

