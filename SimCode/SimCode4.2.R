setwd("C:/Users/TMurray/Documents/PwL_Surv_Model/RCode/SimResults/SimRes4.2")
library(R2jags);library(survival)

###################
## JAGS Programs ##
###################
#Piecewise Exponential Proportional Hazards Model
write("data{
	## Data Construction for the Current Likelihood ##
	for (i in 1:N) {
		for (k in 1:K) {
		# indicates event-time in interval k
		 zed[i,k] <- delta[i]*step(t[i] - tilde.t[k] - 0.00001)*step(tilde.t[k+1] - t[i]);   
		# length of overlap of y[i] with interval k
		 theta[i,k] <- (min(t[i], tilde.t[k+1]) - tilde.t[k])*step(t[i] - tilde.t[k]);      
		}
	} 
}
model{
	 ## Current Likelihoods ##
	for (i in 1:N) {
		for (k in 1:K) {
		# define the likelihood
		 mu[i,k] <- theta[i,k]*exp(gamma[k]+beta*trt[i])
		 zed[i,k] ~ dpois(mu[i,k]);
		}
	} 
	# Prior on hazard parameters
		 beta ~ dnorm(0,0.0001)
		 gamma[1] ~ dnorm(0,0.0001)
		for(k in 2:K){
			gamma[k] ~ dnorm(gamma[k-1],tau.gamma)
	 	}
	# Hyperpriors on historical variance parameters
		 tau.gamma <- 1/(sigma.gamma*sigma.gamma)
		 sigma.gamma ~ dunif(0.01,100)
}","PEPHMod.txt")

#Piecewise Exponential Time-dependent Model
write("data{
	## Data Construction for the Current Likelihood ##
	for (i in 1:N) {
		for (k in 1:K) {
		# indicates event-time in interval k
		 zed[i,k] <- delta[i]*step(t[i] - tilde.t[k] - 0.00001)*step(tilde.t[k+1] - t[i]);   
		# length of overlap of y[i] with interval k
		 theta[i,k] <- (min(t[i], tilde.t[k+1]) - tilde.t[k])*step(t[i] - tilde.t[k]);      
		}
	} 
}
model{
	 ## Current Likelihoods ##
	for (i in 1:N) {
		for (k in 1:K) {
		# define the likelihood
		 mu[i,k] <- theta[i,k]*exp(gamma[1,k]+gamma[2,k]*trt[i])
		 zed[i,k] ~ dpois(mu[i,k]);
		}
	} 
	# Prior on hazard parameters
	for(q in 1:2){
		 gamma[q,1] ~ dnorm(0,0.0001)
		for(k in 2:K){
			gamma[q,k] ~ dnorm(gamma[q,k-1],tau.gamma[q])
	 	}
	# Hyperpriors on historical variance parameters
		 tau.gamma[q] <- 1/(sigma.gamma[q]*sigma.gamma[q])
		 sigma.gamma[q] ~ dunif(0.01,100)
	}
}","PETDMod.txt")

#Piecewise Linear Proportional Hazards Model
write("model{
	for (i in 1:N){
		 lhaz[i] <- inprod(alpha[],T[i,]) + beta*trt[i]
		for(k in 1:K){
		 Haz[i,k] <-  exp(inprod(alpha[],TT[i,k,]))*(1-exp(-(TT[i,k,2]-tilde.t[k])*inprod(alpha[],U[k,])))/inprod(alpha[],U[k,])
		}
		 neg.LL[i] <- -delta[i]*lhaz[i] + sum(Haz[i,])*exp(beta*trt[i]) + C
		 zeros[i] ~ dpois(neg.LL[i]);
	}
	# Prior on hazard parameters
		 beta ~ dnorm(0,0.0001)
		 alpha[1] ~ dnorm(0,0.0001)
		 alpha[2] ~ dnorm(0,0.0001)
		for(k in 2:K){
			alpha[k+1] ~ dnorm(0,tau.alpha)
	 	}
	# Hyperpriors on historical variance parameters
		 tau.alpha <- 1/(sigma.alpha*sigma.alpha)
		 sigma.alpha ~ dunif(0.01,100)
}","PLPHMod.txt")

#Piecewise Linear Time-dependent Model
write("model{
	for (i in 1:N){
		 tlp[i,1:(K+1)] <- alpha[1,] + alpha[2,]*trt[i]
		 lhaz[i] <- inprod(tlp[i,],T[i,])
		for(k in 1:K){
		 Haz[i,k] <-  exp(inprod(tlp[i,],TT[i,k,]))*(1-exp(-(TT[i,k,2]-tilde.t[k])*inprod(tlp[i,],U[k,])))/inprod(tlp[i,],U[k,])
		}
		 neg.LL[i] <- -delta[i]*lhaz[i] + sum(Haz[i,]) + C
		 zeros[i] ~ dpois(neg.LL[i]);
	}
	# Prior on hazard parameters
	for(q in 1:2){
		 alpha[q,1] ~ dnorm(0,0.0001)
		 alpha[q,2] ~ dnorm(0,0.0001)
		for(k in 2:K){
			alpha[q,k+1] ~ dnorm(0,tau.alpha[q])
	 	}
	# Hyperpriors on historical variance parameters
		 tau.alpha[q] <- 1/(sigma.alpha[q]*sigma.alpha[q])
		 sigma.alpha[q] ~ dunif(0.01,100)
	}
}","PLTDMod.txt")

#########################
## Simulation Function ##
#########################
sim = function(N,K,scen){
inverse = function (f, lower = -100, upper = 100){
   function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}

#True Hazard in Control Group
hc = function(t) 1.5*sin(pi*t) + 1
Hc = function(t) -1.5*cos(pi*t)/pi + 1.5/pi + t
Hc.inv = inverse(function(t) ifelse(t < 1, Hc(t) , 2*(t-1) + Hc(1)) , 0, 100)

if(scen==1){
 #Time-Dependent Scenario Hazard and Cumulative Hazard Specifications
 ht = function(t) 4*(.5*cos(pi*t+pi/2) + 1.5/(t+1) - .2) #Time-Dependent Scenario
 Ht = function(t) 4*(.5*cos(pi*t)/pi + 1.5*log(t+1) - .2*t - 1/(2*pi)) #Time-Dependent Scenario
}
if(scen==2){
 #Time-Independent Scenario Hazard and Cumulative Hazard Specifications
 ht = function(t) hc(t)*exp(1) #Time-Independent Scenario
 Ht = function(t) Hc(t)*exp(1) #Time-Independent Scenario
}

Ht.inv = inverse(function (t) ifelse(t < 1, Ht(t), 2*(t-1) + Ht(1)) , 0, 100)
hr = function(x){  ht(x)/hc(x) }

#Generate Data
trt = rbinom(N,1,.5)
#trt = rbinom(N,1,.5)-.5
u = runif(N)
tc = as.vector(unlist(sapply(-log(u),function(z) Hc.inv(z))))
tt = as.vector(unlist(sapply(-log(u),function(z) Ht.inv(z))))
c = runif(N,0,5)
y = ifelse(trt==0,tc,tt)
t = pmin(y,c,rep(1,N))
delta = ifelse(y==t,1,0)
#Time Axis Partition
tilde.t = 0:K/K

####################
## Cox's PH Model ##
####################
cox.fit = coxph(Surv(t,delta)~trt)

##################################
## Piecewise Exponential Models ##
##################################
##PE-PH Model
dta <- list("N", "K", "tilde.t", "t", "delta", "trt")
ints <- function() {list(sigma.gamma=runif(1,0.01,100), gamma=rnorm(K), beta=rnorm(1))}
pars <- c("gamma","beta")

PEPH.fit <- jags(data=dta,inits=ints,model.file="PEPHMod.txt",parameters=pars,
n.chains=2,n.iter=12000,n.burnin=2000,n.thin=1)

## PE-TD Model
ints <- function() {list(sigma.gamma=runif(2,0.01,100), gamma=matrix(rnorm(2*K),2,K))}
pars <- c("gamma")

PETD.fit <- jags(data=dta,inits=ints,model.file="PETDMod.txt",parameters=pars,
n.chains=2,n.iter=12000,n.burnin=2000,n.thin=1)

#############################
## Piecewise Linear Models ##
#############################
#Time Transformation Matrix
OMEGA_alpha <- abs(outer(tilde.t[-c(1,K+1)],tilde.t[-c(1,K+1)],"-"))
svd.OMEGA_alpha<-svd(OMEGA_alpha)
sqrt.OMEGA_alpha<-t(svd.OMEGA_alpha$v%*%(t(svd.OMEGA_alpha$u)*sqrt(svd.OMEGA_alpha$d)))
inv.D <- solve(cbind(c(1,rep(0,K)),c(0,1,rep(0,K-1)),rbind(rep(0,K-1),rep(0,K-1),sqrt.OMEGA_alpha)))

#Construct Time Design Matrices
T_K = cbind(1,t(sapply(t,function(z) abs(z - tilde.t[-c(K+1)]) - tilde.t[-c(K+1)])))
T = T_K%*%inv.D
sk = t(sapply(t,function(z) pmax(pmin(z,tilde.t[-1]),tilde.t[-(K+1)])))
TT_K = TT = array(NA,c(N,K,K+1))
for(k in 1:K) TT_K[,k,] <- cbind(1,t(sapply(sk[,k],function(z) abs(z - tilde.t[-c(K+1)]) - tilde.t[-c(K+1)])))
for(i in 1:N) TT[i,,] = TT_K[i,,]%*%inv.D
U_K <- matrix(1, K, K); U_K[upper.tri(U_K)] <- -1; U <- U_K%*%inv.D[-1,]

## PL-PH Model
zeros = rep(0,N); C = 1000
dta <- list("N", "K", "tilde.t", "T", "TT", "U", "delta", "trt", "zeros", "C")
ints <- function() {list(sigma.alpha=runif(1,0.01,100), alpha=rnorm(K+1), beta=rnorm(1))}
pars <- c("alpha", "beta")

PLPH.fit <- jags(data=dta,inits=ints,model.file="PLPHMod.txt",parameters=pars,
n.chains=2,n.iter=12000,n.burnin=2000,n.thin=1)

## PL-TD Model 
ints <- function() {list(sigma.alpha=runif(2,0.01,100), alpha=matrix(rnorm(2*(K+1)),2,K+1))}
pars <- c("alpha")

PLTD.fit <- jags(data=dta,inits=ints,model.file="PLTDMod.txt",parameters=pars,
n.chains=2,n.iter=12000,n.burnin=2000,n.thin=1)

#Write parameter estimtates to files
write.table(rbind(round(cox.fit$coeff,4)),paste("CoxPHres",scen,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep=',',append=TRUE)
write.table(rbind(round(c(PEPH.fit$BUGSoutput$mean$gamma,PEPH.fit$BUGSoutput$mean$beta),4)),paste("PEPHres",scen,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep=',',append=TRUE)
write.table(rbind(round(c(PETD.fit$BUGSoutput$mean$gamma[1,],PETD.fit$BUGSoutput$mean$gamma[2,]),4)),paste("PETDres",scen,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep=',',append=TRUE)
write.table(rbind(round(c(PLPH.fit$BUGSoutput$mean$alpha,PLPH.fit$BUGSoutput$mean$beta),4)),paste("PLPHres",scen,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep=',',append=TRUE)
write.table(rbind(round(c(PLTD.fit$BUGSoutput$mean$alpha[1,],PLTD.fit$BUGSoutput$mean$alpha[2,]),4)),paste("PLTDres",scen,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep=',',append=TRUE)
}

###############################
## Parallel Computation Code ##
###############################
library(snowfall)
sfInit(socketHosts=rep("localhost",8),cpus=8,type='SOCK',parallel=TRUE)
sfLibrary(R2jags);sfLibrary(survival)
sfExport('sim')
wrapper = function(x){sim(N=200)}
res = sfLapply(1:200,wrapper)
sfStop()

######################
##  Compile Results ##
######################
#Scenario
scen = 1
##True Hazard in Control Group
hc = function(t) 1.5*sin(pi*t) + 1
Hc = function(t) -1.5*cos(pi*t)/pi + 1.5/pi + t
if(scen==1){
 ##Time-Dependent Scenario 
 ht = function(t) 4*(.5*cos(pi*t+pi/2) + 1.5/(t+1) - .2) #Time-Dependent Scenario
 Ht = function(t) 4*(.5*cos(pi*t)/pi + 1.5*log(t+1) - .2*t - 1/(2*pi)) #Time-Dependent Scenario
}
if(scen==2){
 ##Proportional Hazards Scenario
 ht = function(t) hc(t)*exp(1) #Time-Independent Scenario
 Ht = function(t) Hc(t)*exp(1) #Time-Independent Scenario
}
##Treatment versus Control Hazard Ratio
hr = function(x) ht(x)/hc(x)

#Read in parameter estimates for all simulation runs
CoxPHres = as.matrix(read.csv(paste("CoxPHres",scen,".txt",sep=""),header=FALSE))
PEPHres = as.matrix(read.csv(paste("PEPHres",scen,".txt",sep=""),header=FALSE))
PETDres = as.matrix(read.csv(paste("PETDres",scen,".txt",sep=""),header=FALSE))
PLPHres = as.matrix(read.csv(paste("PLPHres",scen,".txt",sep=""),header=FALSE))
PLTDres = as.matrix(read.csv(paste("PLTDres",scen,".txt",sep=""),header=FALSE))

#Partition
K = dim(PEPHres)[2]-1; tilde.t = c(0:K/K)

##control log-hazard estimates
t = seq(0,1,length.out=10000)
#PE models
S = t(sapply(t,function(t) ifelse(t<c(tilde.t[-c(1,K+1)],1.1),1,0)*ifelse(t>=tilde.t[-(K+1)],1,0)))
PEPH.lhaz = t(S%*%t(PEPHres[,1:K]))
PETD.lhaz = t(S%*%t(PETDres[,1:K]))

#PL models
#Time Transformation Matrix
OMEGA_alpha <- abs(outer(tilde.t[-c(1,K+1)],tilde.t[-c(1,K+1)],"-"))
svd.OMEGA_alpha<-svd(OMEGA_alpha)
sqrt.OMEGA_alpha<-t(svd.OMEGA_alpha$v%*%(t(svd.OMEGA_alpha$u)*sqrt(svd.OMEGA_alpha$d)))
inv.D <- solve(cbind(c(1,rep(0,K)),c(0,1,rep(0,K-1)),rbind(rep(0,K-1),rep(0,K-1),sqrt.OMEGA_alpha)))
#Time Design Matrix
T_K = cbind(1,t(sapply(t,function(z) abs(z - tilde.t[-c(K+1)]) - tilde.t[-c(K+1)])))
T = T_K%*%inv.D
PLPH.lhaz = t(T%*%t(PLPHres[,1:(K+1)]))
PLTD.lhaz = t(T%*%t(PLTDres[,1:(K+1)]))
#Average and Quantile Estimates
PEPH.lhaz.hat = apply(PEPH.lhaz,2,function(x) c(mean(x),quantile(x,c(0.025,0.975))))
PETD.lhaz.hat = apply(PETD.lhaz,2,function(x) c(mean(x),quantile(x,c(0.025,0.975))))
PLPH.lhaz.hat = apply(PLPH.lhaz,2,function(x) c(mean(x),quantile(x,c(0.025,0.975))))
PLTD.lhaz.hat = apply(PLTD.lhaz,2,function(x) c(mean(x),quantile(x,c(0.025,0.975))))
#RISE of average estimate and best individual fit
RISE.lhaz = c(sqrt(mean((PEPH.lhaz.hat[1,]-log(hc(t)))**2)),sqrt(mean((PETD.lhaz.hat[1,]-log(hc(t)))**2)),sqrt(mean((PLPH.lhaz.hat[1,]-log(hc(t)))**2)),sqrt(mean((PLTD.lhaz.hat[1,]-log(hc(t)))**2)))
best.fit.lhaz = c(which.min(apply(PEPH.lhaz,1,function(x) sqrt(mean((x-log(hc(t)))**2)))),which.min(apply(PETD.lhaz,1,function(x) sqrt(mean((x-log(hc(t)))**2)))),which.min(apply(PLPH.lhaz,1,function(x) sqrt(mean((x-log(hc(t)))**2)))),which.min(apply(PLTD.lhaz,1,function(x) sqrt(mean((x-log(hc(t)))**2)))))

##log-hazard ratio estimates
#Cox PH model
CoxPH.lhr = as.vector(CoxPHres)
#PE models
PEPH.lhr = PEPHres[,K+1]
PETD.lhr = t(S%*%t(PETDres[,-c(1:K)]))
#PL models
PLPH.lhr = PLPHres[,K+2]
PLTD.lhr = t(T%*%t(PLTDres[,-c(1:(K+1))]))
#Average and Quantile Estimates
CoxPH.lhr.hat = c(mean(CoxPH.lhr),quantile(CoxPH.lhr,c(0.025,0.975)))
PEPH.lhr.hat = c(mean(PEPH.lhr),quantile(PEPH.lhr,c(0.025,0.975)))
PETD.lhr.hat = apply(PETD.lhr,2,function(x) c(mean(x),quantile(x,c(0.025,0.975))))
PLPH.lhr.hat = c(mean(PLPH.lhr),quantile(PLPH.lhr,c(0.025,0.975)))
PLTD.lhr.hat = apply(PLTD.lhr,2,function(x) c(mean(x),quantile(x,c(0.025,0.975))))
#RISE of average estimate and best individual fit
RISE.lhr = c(sqrt(mean((CoxPH.lhr.hat[1]-log(hr(t)))**2)),sqrt(mean((PEPH.lhr.hat[1]-log(hr(t)))**2)),sqrt(mean((PETD.lhr.hat[1,]-log(hr(t)))**2)),sqrt(mean((PLPH.lhr.hat[1]-log(hr(t)))**2)),sqrt(mean((PLTD.lhr.hat[1,]-log(hr(t)))**2)))
best.fit.lhr = c(which.min(sapply(CoxPH.lhr,function(x) sqrt(mean((x-log(hr(t)))**2)))),which.min(sapply(PEPH.lhr,function(x) sqrt(mean((x-log(hr(t)))**2)))),which.min(apply(PETD.lhr,1,function(x) sqrt(mean((x-log(hr(t)))**2)))),which.min(sapply(PLPH.lhr,function(x) sqrt(mean((x-log(hr(t)))**2)))),which.min(apply(PLTD.lhr,1,function(x) sqrt(mean((x-log(hr(t)))**2)))))

###############################
# Generate panels in Figure 2 #
###############################
#PETD control log-hazard
pdf(paste("PETDlhaz",scen,".pdf",sep=""),family='serif')
par(mfrow=c(1,1),mar=c(2.7,3.2,.2,.2),las=1,mgp=c(1.5,.5,0),cex.lab=1.5,cex=2.5)
curve(log(hc(x)),main=" ",xlab=expression(t),ylab=expression(paste("log{h(t|0)}")),ylim=c(-2,2),col='white')
polygon(c(c(t[1],rep(t[-1],each=2)),rev(c(rep(t[-length(t)],each=2),t[length(t)]))),
 c(c(rep(PETD.lhaz.hat[3,-length(t)],each=2),PETD.lhaz.hat[3,length(t)]),c(rev(rep(PETD.lhaz.hat[2,-1],each=2)),PETD.lhaz.hat[2,1])),col='grey60',border='grey60')
lines(t,PETD.lhaz.hat[1,],lwd=3,type='s')
lines(t,PETD.lhaz[best.fit.lhaz[2],],lwd=3,type='s',col='grey40')
curve(log(hc(x)),from=0,to=1,add=TRUE,lwd=3,lty=2)
text(.28, 1.95, paste("RISE =",formatC(round(RISE.lhaz[2],2),format='f',digits=2)),cex=1.4)
if(scen==1) legend("bottomleft",c("Truth","Average Estimate","Best Estimate"),lwd=3,col=c("black","black","grey40"),lty=c(2,1,1),bty='n',cex=1.2,seg.len=1.3)
dev.off()

#PLTD control log-hazard
pdf(paste("PLTDlhaz",scen,".pdf",sep=""),family='serif')
par(mfrow=c(1,1),mar=c(2.7,3.2,.2,.2),las=1,mgp=c(1.5,.5,0),cex.lab=1.5,cex=2.5)
curve(log(hc(x)),main=" ",xlab=expression(t),ylab=expression(paste("log{h(t|0)}")),ylim=c(-2,2),col='white')
polygon(c(c(t[1],rep(t[-1],each=2)),rev(c(rep(t[-length(t)],each=2),t[length(t)]))),
 c(c(rep(PLTD.lhaz.hat[3,-length(t)],each=2),PLTD.lhaz.hat[3,length(t)]),c(rev(rep(PLTD.lhaz.hat[2,-1],each=2)),PLTD.lhaz.hat[2,1])),col='grey60',border='grey60')
lines(t,PLTD.lhaz.hat[1,],lwd=3)
lines(t,PLTD.lhaz[best.fit.lhaz[4],],lwd=3,col='grey40')
curve(log(hc(x)),from=0,to=1,add=TRUE,lwd=3,lty=2)
text(.28, 1.95, paste("RISE =",formatC(round(RISE.lhaz[4],2),format='f',digits=2)),cex=1.4)
dev.off()

#PETD log-hazard ratio
pdf(paste("PETDlhr",scen,".pdf",sep=""),family='serif')
par(mfrow=c(1,1),mar=c(2.7,3.2,.2,.2),las=1,mgp=c(1.5,.5,0),cex.lab=1.5,cex=2.5)
curve(log(hr(x)),main=" ",xlab=expression(t),ylab=expression(paste("log{h(t|1)/h(t|0)}")),ylim=c(-2.5,4),col='white')
polygon(c(c(t[1],rep(t[-1],each=2)),rev(c(rep(t[-length(t)],each=2),t[length(t)]))),
 c(c(rep(PETD.lhr.hat[3,-length(t)],each=2),PETD.lhr.hat[3,length(t)]),c(rev(rep(PETD.lhr.hat[2,-1],each=2)),PETD.lhr.hat[2,1])),col='grey60',border='grey60')
lines(t,PETD.lhr.hat[1,],lwd=3,type='s')
lines(t,PETD.lhr[best.fit.lhr[3],],lwd=3,type='s',col='grey40')
curve(log(hr(x)),from=0,to=1,add=TRUE,lwd=3,lty=2)
text(.28, 3.9, paste("RISE =",formatC(round(RISE.lhr[3],2),format='f',digits=2)),cex=1.4)
dev.off()

#PLTD log-hazard ratio
pdf(paste("PLTDlhr",scen,".pdf",sep=""),family='serif')
par(mfrow=c(1,1),mar=c(2.7,3.2,.2,.2),las=1,mgp=c(1.5,.5,0),cex.lab=1.5,cex=2.5)
curve(log(hr(x)),main=" ",xlab=expression(t),ylab=expression(paste("log{h(t|1)/h(t|0)}")),ylim=c(-2.5,4),col='white')
polygon(c(c(t[1],rep(t[-1],each=2)),rev(c(rep(t[-length(t)],each=2),t[length(t)]))),
 c(c(rep(PLTD.lhr.hat[3,-length(t)],each=2),PLTD.lhr.hat[3,length(t)]),c(rev(rep(PLTD.lhr.hat[2,-1],each=2)),PLTD.lhr.hat[2,1])),col='grey60',border='grey60')
lines(t,PLTD.lhr.hat[1,],lwd=3)
lines(t,PLTD.lhr[best.fit.lhr[5],],lwd=3,col='grey40')
curve(log(hr(x)),from=0,to=1,add=TRUE,lwd=3,lty=2)
text(.28, 3.9, paste("RISE =",formatC(round(RISE.lhr[5],2),format='f',digits=2)),cex=1.4)
dev.off()

########################################################################
# Table 1: Distributional summaries of individual estimate RISE values #
########################################################################
CoxPH.RISEs = sapply(CoxPH.lhr,function(x) sqrt(mean((x-log(hr(t)))**2)))
PEPH.RISEs = sapply(PEPH.lhr,function(x) sqrt(mean((x-log(hr(t)))**2)))
PETD.RISEs = apply(PETD.lhr,1,function(x) sqrt(mean((x-log(hr(t)))**2)))
PLPH.RISEs = sapply(PLPH.lhr,function(x) sqrt(mean((x-log(hr(t)))**2)))
PLTD.RISEs = apply(PLTD.lhr,1,function(x) sqrt(mean((x-log(hr(t)))**2)))

#Show Results
round(t(apply(rbind(CoxPH.RISEs,PEPH.RISEs,PETD.RISEs,PLPH.RISEs,PLTD.RISEs),1,function(x) c(mean(x),sd(x),quantile(x,c(.025,.25,.5,.75,.975))) )),2)

#Limit to t in [0,0.50]
CoxPH.RISEs = sapply(CoxPH.lhr,function(x) sqrt(mean((x-log(hr(t[1:5000])))**2)))
PEPH.RISEs = sapply(PEPH.lhr,function(x) sqrt(mean((x-log(hr(t[1:5000])))**2)))
PETD.RISEs = apply(PETD.lhr[,1:5000],1,function(x) sqrt(mean((x-log(hr(t[1:5000])))**2)))
PLPH.RISEs = sapply(PLPH.lhr,function(x) sqrt(mean((x-log(hr(t[1:5000])))**2)))
PLTD.RISEs = apply(PLTD.lhr[,1:5000],1,function(x) sqrt(mean((x-log(hr(t[1:5000])))**2)))
round(t(apply(rbind(CoxPH.RISEs,PEPH.RISEs,PETD.RISEs,PLPH.RISEs,PLTD.RISEs),1,function(x) c(mean(x),sd(x),quantile(x,c(.025,.25,.5,.75,.975))) )),2)
