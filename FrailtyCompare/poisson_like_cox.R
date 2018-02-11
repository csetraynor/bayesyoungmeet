#Cox proportional Hazard
library(survival)
library(mgcv)
library(dplyr)
library(data.table)
library(mgcv)
library(stargazer)
library(DT)
library(texreg)

N <-  42 
NT <-  17
obs_t <- c(1, 1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 11, 11, 12, 12, 15, 
    17, 22, 23, 6, 6, 6, 6, 7, 9, 10, 10, 11, 13, 16, 17, 19, 20, 
    22, 23, 25, 32, 32, 34, 35) 
fail <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 
    0) 
Z <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, 
    -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 
    -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5) 
t <- c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 22, 23, 
    35) 
##Cox as Poisson

data<-data.frame(ID=1:length(obs_t),Time=obs_t,Status=fail,X=Z)
dim(data)

datatable(data)

#Expand dataset accordin to calendar time
for (i in 1:max(data$ID)){
  stpT<-1:data$Time[i]
  id<-rep(i,length(stpT))
  stat<-c(rep(0,length(stpT)-1),data$Status[i])                            
  strT<-lag(stpT,1);strT[1]=0 
  iln<-stpT-strT
  
  df<-data.frame(ID=id,Start=strT,Stop=stpT,Status=stat,ILen=iln)
  if(i==1){data_cal=df}
  else{data_cal=rbind(data_cal,df)}
}
data_cal<-merge(data_cal,data[,c('ID','X')],by='ID')
dim(data_cal)

summary(data_cal)

#Customise dataset by cutting all intervals where no event has occured
# create subset of cal-time dataset - i.e., remove all non-event times
D<-as.data.table(data_cal)
clear<-unlist(D[,sum(Status),by=Start][V1==0 & Start!=0,Start])
data_corcal<-subset(data_cal, !(Start %in% clear))

# redo start-time and interval length calculations
D<-as.data.table(data_corcal)
D[,strt:=lag(Stop,1),by=ID][is.na(strt),strt:=0][,Ilen2:=Stop-strt]

# get cleaner dataset
data_corcal<-data.frame(ID=D$ID,Start=D$strt,Stop=D$Stop,Status=D$Status,ILen=D$Ilen2)
data_corcal<-merge(data_corcal,data[,c('ID','X')],by='ID')

dim(data_corcal)
datatable(data_corcal)

#Expand accord to observation time
for (i in 1:max(data$ID)){
  obst<-sort(unique(data$Time))
  stpT<-obst[1:which(obst==data$Time[i])]
  id<-rep(i,length(stpT))
  stat<-c(rep(0,length(stpT)-1),data$Status[i])                            
  strT<-lag(stpT,1);strT[1]=0 
  iln<-stpT-strT
  
  df<-data.frame(ID=id,Start=strT,Stop=stpT,Status=stat,ILen=iln)
  if(i==1){data_obs=df}
  else{data_obs=rbind(data_obs,df)}
}
data_obs<-merge(data_obs,data[,c('ID','X')],by='ID')
dim(data_obs)

###Cox Proportional Hazard
###Standard Cox proportional Hazard
m0<-coxph(Surv(Time,Status)~1+X,data)
m1<-coxph(Surv(Start,Stop,Status)~1+X,data_obs)
m2<-coxph(Surv(Start,Stop,Status)~1+X,data_cal)
m3<-coxph(Surv(Start,Stop,Status)~1+X,data_corcal)

######Poisson trick

#Important aspects of the modelling approach: * Dependent variable is the binary event indicator. * Adding an intercept is recommendable to provide a decent and interpretable reference category for factor variables.
#Use the logarithm of the interval length as an offset. * Add a smoothing-spline to model the baseline hazard.

m1p<-gam(Status~1+offset(log(ILen))+s(Stop)+X,data_obs,family='poisson')
m2p<-gam(Status~1+offset(log(ILen))+s(Stop)+X,data_cal,family='poisson')
m3p<-gam(Status~1+offset(log(ILen))+s(Stop)+X,data_corcal,family='poisson')

m1p2<-gam(Status~1+offset(log(ILen))+s(Start)+X,data_obs,family='poisson')
m2p2<-gam(Status~1+offset(log(ILen))+s(Start)+X,data_cal,family='poisson')
m3p2<-gam(Status~1+offset(log(ILen))+s(Start)+X,data_corcal,family='poisson')
texreg(list(Obs_Time=m1p2,Obs_Time=m1p,Calendar_Time=m2p2,Calendar_Time=m2p,Custom_Time=m3p2,Custom_Time=m3p))

#** Deciding on which model to use: looking at performance for easy models (baseline hazards) **
T<-c(sort(unique(data$Time)))

p1<-predict(gam(Status~1+offset(log(ILen))+s(Stop),data_obs,family='poisson'),data.frame(ILen=1,Stop=T))
S1<-exp(-cumsum(exp(p1)))

p2<-predict(gam(Status~1+offset(log(ILen))+s(Stop),data_cal,family='poisson'),data.frame(ILen=1,Stop=T))
S2<-exp(-cumsum(exp(p2)))

p3<-predict(gam(Status~1+offset(log(ILen))+s(Stop),data_corcal,family='poisson'),data.frame(ILen=1,Stop=T))
S3<-exp(-cumsum(exp(p3)))
plotframe<-data.frame(Time=c(0,T),S1=c(1,S1),S2=c(1,S2),S3=c(1,S3))

#plot that show the baseline survival estimates close to kaplan and meier estimates
plot(survfit(coxph(Surv(Time,Status)~1,data)),ylab='Survival Probability',xlab='Time',
     main='Baseline Survival Curves')
lines(plotframe$Time,plotframe$S1,col='green')
lines(plotframe$Time,plotframe$S2,col='red')
lines(plotframe$Time,plotframe$S3,col='blue')
legend('topright',c('Observation Time','Calendar Time','Customized Calendar Time'),bty='n',lty=1,col=c('green','red','blue'))
plotframe<-data.frame(Time=c(0,T),S1=c(1,S1),S2=c(1,S2),S3=c(1,S3))