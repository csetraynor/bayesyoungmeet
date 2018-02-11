library('ggplot2') # visualization
library('ggthemes') # visualization
library('ggridges') # visualization
library('ggforce') # visualization
library('ggExtra') # visualization
library('GGally') # visualisation
library('scales') # visualization
library('grid') # visualisation
library('gridExtra') # visualisation
library('corrplot') # visualisation
library('VIM') # missing values
#suppressPackageStartupMessages(library(heatmaply)) # visualisation
library('dplyr') # data manipulation
library('tidyr') # data manipulation
library('readr') # data input
library('stringr') # string manipulation
library('forcats') # factor manipulation
library('modelr') # factor manipulation
library('randomForest') # classification
library('xgboost') # classification
library('ROCR') # model validation

  
summary(g_clin_tib)
glimpse(g_clin_tib)

#Missing values
library(VIM)

aggr(g_clin_tib, prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

g_clin_tib <- g_clin_tib %>%
  select(-CANCER_TYPE, -CANCER_TYPE_DETAILED, -ONCOTREE_CODE) %>%   ##Non informative all are Gliobastome Multiforme and GBM
  filter(!is.na(DFS_MONTHS))
aggr(g_clin_tib, prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

##sample
set.seed(1573)
#g_clin_tib <- g_clin_tib[sample(nrow(g_clin_tib), 200),]

#########################################################
##########Select Clinical Variables
g_clin_tib <- rownames_to_column(g_clin_tib, var = "index")
x <- g_clin_tib %>%
  mutate(
    TMZTherapy = index %in% (starts_with("TMZ", vars = THERAPY_CLASS)),
    UnspecifiedTherapy = index %in% (starts_with("Unspecified", vars = THERAPY_CLASS)),
    NonstandardTherapy = index %in% (starts_with("Nonstandard", vars = THERAPY_CLASS)),
    AlkylatingTherapy = index %in% (starts_with("Alkylating", vars = THERAPY_CLASS))
  ) %>%
  select(TMZTherapy, UnspecifiedTherapy, NonstandardTherapy, AlkylatingTherapy, AGE, SEX)

x <- model.matrix(~., x)

#Select response variates

y <- g_clin_tib$DFS_MONTHS

event <- g_clin_tib %>% mutate(
  event = (DFS_STATUS == "Recurred/Progressed")
)  %>%
  select(event) %>% unlist()

############ create data set

N <- length(y)
M <- ncol(x)
data <-  list(x = x,
              y = y,
              event = event,
              N = N,
              M =M)