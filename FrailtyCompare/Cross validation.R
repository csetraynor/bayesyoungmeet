#Cross validation and Monte Carlo Cross Validation

library(WilcoxCV)

#generate CV
?generate.cv
g_clin_tib_cv <- generate.cv(n = 541, m = 10)

#generate MCCV
?generate.split
g_clin_tib_mccv <- generate.split(niter = 50, n = 541, ntest= 180)

g_clin_tib_mccv[1,]

#test survival
require(survival)
library(tidyverse)
#get sample id
g_clin_tib <- as.numeric(rownames_to_column(g_clin_tib, var = 'sample'))

#write not in function
'!%in%' <- function(x,y)!('%in%'(x,y))

K <- 2
for(i in 1:K){
  fit <- coxph(Surv(OS_MONTHS, IsDeceased) ~ GCIMPorCL6 + CombTherapy + Agecentered , data = g_clin_tib[notin(g_clin_tib$sample , g_clin_tib_mccv[1,]),])
}

g_clin_tib[g_clin_tib$sample, g_clin_tib_mccv[1,])

notin <- function(x,y)!('%in%'(x,y))

fit <- do.call(rbind, lapply(split(g_clin_tib, g_clin_tib_mccv), function(x){
  fit <- coxph(Surv(OS_MONTHS, IsDeceased) ~ GCIMPorCL6 + CombTherapy + Agecentered , data = g_clin_tib[notin(g_clin_tib$sample , g_clin_tib_mccv[x,]),])
}))


res <- do.call(rbind,lapply(split(dat, dat$year),function(x){
  fit <- lm(value~exp(AgeR), data=x)
  res <- data.frame(year=unique(x$year),coeff=coef(fit)[2])
  res
}))


#in addition check rms package http://biostat.mc.vanderbilt.edu/wiki/Main/RmS
