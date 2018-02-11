library(tidyverse)
library(Coxnet)

metabric_clinicaldata <- as.tibble(metabric_clinicaldata)
metabric_clinicaldata[metabric_clinicaldata == "" | metabric_clinicaldata == " "] <- NA 



# create a subset of the data with only our variables of interest (variables
# that aren't converted numbers won't work)
dataSubset <- metabric_clinicaldata %>%
  filter(!is.na(OS_STATUS)) %>%                    #filtrate available observation
  mutate(IsDeceased = (OS_STATUS == "DECEASED")) %>% # convert to boolean
  select(TUMOR_SIZE, TUMOR_STAGE, GRADE, AGE_AT_DIAGNOSIS, COHORT, IsDeceased, OS_MONTHS) %>%
  na.omit() # remove non-numeric values

# convert our input variables to a matrix
input <- dataSubset %>%
  select(TUMOR_SIZE, TUMOR_STAGE,GRADE) %>% # dont include the variable we're predicting!
  as.matrix()

# get a vector with our output variable
output <- Surv(dataSubset$OS_MONTHS, dataSubset$IsDeceased)

#model coxph
require(survival)
coxfit <- coxph(output ~ input)
summary(coxfit)

#output to fit a Coxnet

y <- dataSubset %>%
  select(OS_MONTHS, IsDeceased) %>% 
  rename(time = OS_MONTHS, status = IsDeceased) %>%
  as.matrix()

#fit a Elastic net
netfit <-  Coxnet(input, y, penalty = "Lasso")
print(netfit)
