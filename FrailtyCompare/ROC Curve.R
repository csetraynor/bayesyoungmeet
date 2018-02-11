library(tidyverse)
require(survivalROC)

#obtain a subset of Dataset to evaluate how well NPI predicts survival outcome
dataSubset <- metabric_clinicaldata %>%
  filter(!is.na(OS_STATUS)) %>%                    #filtrate available observation
  mutate(IsDeceased = (OS_STATUS == "DECEASED")) %>% # convert to boolean
  select(NPI, IsDeceased, OS_MONTHS) %>%
  na.omit() 

nobs <- NROW(dataSubset)
cutoff <- 12

#fit ROC curve
npi_fit <- survivalROC(Stime = dataSubset$OS_MONTHS,
                       status = dataSubset$IsDeceased,
                       marker = dataSubset$NPI,
                       predict.time = cutoff, span = 0.25*nobs^(-0.2))

#plot
plot(npi_fit$FP, npi_fit$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.1$AUC,3)),
     ylab="TP",main="NPI, Method = NNE \n Year = 1")
abline(0,1)