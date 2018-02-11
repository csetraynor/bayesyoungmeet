#Missing values
library(VIM)

aggr(g_clin_tib, prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

g_clin_tib <- g_clin_tib %>%
  select(-CANCER_TYPE, -CANCER_TYPE_DETAILED, -ONCOTREE_CODE) %>%   ##Non informative all are Gliobastome Multiforme and GBM
  filter(!is.na(DFS_MONTHS))
aggr(g_clin_tib, prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)
#########################################################

library(data.table)
data <- (g_clin_tib)

#Expand accord to observation time
###ADD DUMMY POINTS
data <- rownames_to_column(data, var = "ID")
data$DFS_STATUS <- (data$DFS_STATUS != "Recurred/Progressed")
data$DFS_MONTHS <- data$DFS_MONTHS * 10
for (i in 1:max(data$ID)){
  obst<-sort(unique(data$DFS_MONTHS))
  stpT<-obst[1:which(obst==data$DFS_MONTHS[i])]
  id<-rep(i,length(stpT))
  stat<-c(rep(0,length(stpT)-1),data$DFS_STATUS[i])                            
  strT<-lag(stpT,1);strT[1]=0 
  iln<-stpT-strT
  
  df<-data.frame(ID=id,Start=strT,Stop=stpT,Status=stat,ILen=iln)
  if(i==1){data_obs=df}
  else{data_obs=rbind(data_obs,df)}
}
data_obs<-merge(data_obs,data,by='ID')
summary(data_obs)

library(DT)
dim(data_obs)
datatable(data_obs)

data_obs$Agecentered <- NULL
g_clin_tib <- as.tibble(data_obs)

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


y <- g_clin_tib %>%
  mutate(
    
  )

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

