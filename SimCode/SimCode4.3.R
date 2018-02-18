setwd("C:/Users/TMurray/Documents/PwL_Surv_Model/RCode/SimResults/SimRes4.3")
library(R2jags)

###################
## JAGS Programs ##
###################
#Piecewise Linear Proportional Hazards Model with a Linear effect
write("model{
	for (i in 1:N){
		 lhaz[i] <- inprod(alpha[],T[i,]) + beta*x[i]
		for(k in 1:K){
		 Haz[i,k] <-  exp(inprod(alpha[],TT[i,k,]))*(1-exp(-(TT[i,k,2]-tilde.t[k])*inprod(alpha[],U[k,])))/inprod(alpha[],U[k,])
		}
		 neg.LL[i] <- -delta[i]*lhaz[i] + sum(Haz[i,])*exp(beta*x[i]) + C
		 zeros[i] ~ dpois(neg.LL[i]);
	}
	# Prior on hazard parameters
	 alpha[1] ~ dnorm(0,0.0001)
	 alpha[2] ~ dnorm(0,0.0001)
	for(k in 2:K){
	 alpha[k+1] ~ dnorm(0,tau.alpha)
	}
	# Prior on regression parameters
	 beta ~ dnorm(0,0.0001)
	# Hyperpriors on historical variance parameters
	 tau.alpha <- 1/(sigma.alpha*sigma.alpha)
	 sigma.alpha ~ dunif(0.01,100)
}","LinearMod.txt")

#Piecewise Linear Proportional Hazards Model with an Smooth Effect
write("model{
	for (i in 1:N){
		 lhaz[i] <- inprod(alpha[],T[i,]) + inprod(beta[],X[i,])
		for(k in 1:K){
		 Haz[i,k] <-  exp(inprod(alpha[],TT[i,k,]))*(1-exp(-(TT[i,k,2]-tilde.t[k])*inprod(alpha[],U[k,])))/inprod(alpha[],U[k,])
		}
		 neg.LL[i] <- -delta[i]*lhaz[i] + sum(Haz[i,])*exp(inprod(beta[],X[i,])) + C
		 zeros[i] ~ dpois(neg.LL[i]);
	}
	# Prior on hazard parameters
	 alpha[1] ~ dnorm(0,0.0001)
	 alpha[2] ~ dnorm(0,0.0001)
	for(k in 2:K){
	 alpha[k+1] ~ dnorm(0,tau.alpha)
	}
	# Prior on regression parameters
	 beta[1] ~ dnorm(0,0.0001)
	for(j in 2:J){
	 beta[j] ~ dnorm(0,tau.beta)
	}
	# Hyperpriors on historical variance parameters
	 tau.alpha <- 1/(sigma.alpha*sigma.alpha)
	 sigma.alpha ~ dunif(0.01,100)
	 tau.beta <- 1/(sigma.beta*sigma.beta)
	 sigma.beta ~ dunif(0.01,100)
}","SmoothMod.txt")

#Piecewise Linear Proportional Hazards Model with a Shape-restricted Effect
write("model{
	for (i in 1:N){
		 lhaz[i] <- inprod(alpha[],T[i,]) + inprod(beta[],X2[i,])
		for(k in 1:K){
		 Haz[i,k] <-  exp(inprod(alpha[],TT[i,k,]))*(1-exp(-(TT[i,k,2]-tilde.t[k])*inprod(alpha[],U[k,])))/inprod(alpha[],U[k,])
		}
		 neg.LL[i] <- -delta[i]*lhaz[i] + sum(Haz[i,])*exp(inprod(beta[],X2[i,])) + C
		 zeros[i] ~ dpois(neg.LL[i]);
	}
	# Prior on hazard parameters
	 alpha[1] ~ dnorm(0,0.0001)
	 alpha[2] ~ dnorm(0,0.0001)
	for(k in 2:K){
	 alpha[k+1] ~ dnorm(0,tau.alpha)
	}
	# Prior on regression parameters
	for(j in 1:(J+1)){
	 beta[j] <- iota[j]*beta.star[j]
	 beta.star[j] ~ dnorm(0,0.0001)T(0,)
	 iota[j] ~ dbern(0.5)
	}
	# Hyperpriors on historical variance parameters
	 tau.alpha <- 1/(sigma.alpha*sigma.alpha)
	 sigma.alpha ~ dunif(0.01,100)
}","RestrictMod.txt")

#########################
## Simulation Function ##
#########################
sim = function(N,K,J,scen){
inverse = function (f, lower = -100, upper = 100){
   function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}

 #Flat Nonlinear Scenario
 if(scen==1) f = function(x) 5*(exp(20*(x-0.5))/(1+exp(20*(x-0.5)))-0.5)
 #Smooth Nonlinear Scenario
 if(scen==2) f = function(x) 20*(x-.5)**3  
 #Linear Scenario
 if(scen==3) f = function(x) 4*(x-.5)
 
 #True Baseline Hazard Function 
 h = function(t) 4*(.5*cos(pi*t+pi/2) + 1.5/(t+1) - log(1.125))
 H = function(t) 4*(.5*cos(pi*t)/pi + 1.5*log(t+1) - log(1.125)*t - 1/(2*pi))
 H.inv = inverse(function (t) ifelse(t < 1, H(t), 2*(t-1) + H(1)) , 0, 100)

 #Generate Data
 x = seq(0,1,length.out=N)
 u = runif(N)
 y = as.vector(unlist(sapply(-log(u)*exp(-f(x)),function(z) H.inv(z))))
 c = runif(N,0,5)
 t = pmin(y,c,rep(1,N))
 delta = ifelse(t==y,1,0)

###############
## PL Models ##
###############
#Time Axis Partition
tilde.t = 0:K/K

#Time Transformation Matrix
OMEGA_alpha <- abs(outer(tilde.t[-c(1,K+1)],tilde.t[-c(1,K+1)],"-"))
svd.OMEGA_alpha<-svd(OMEGA_alpha)
sqrt.OMEGA_alpha<-t(svd.OMEGA_alpha$v%*%(t(svd.OMEGA_alpha$u)*sqrt(svd.OMEGA_alpha$d)))
inv.D <- solve(cbind(c(1,rep(0,K)),c(0,1,rep(0,K-1)),rbind(rep(0,K-1),rep(0,K-1),sqrt.OMEGA_alpha)))

#Construct Time Design Matrices
T_K = cbind(1,t(sapply(t,function(z) abs(z - tilde.t[-c(K+1)]) - tilde.t[-c(K+1)])))
T = T_K%*%inv.D
tk = t(sapply(t,function(z) pmax(pmin(z,tilde.t[-1]),tilde.t[-(K+1)])))
TT_K = TT = array(NA,c(N,K,K+1))
for(k in 1:K) TT_K[,k,] <- cbind(1,t(sapply(tk[,k],function(z) abs(z - tilde.t[-c(K+1)]) - tilde.t[-c(K+1)])))
for(i in 1:N) TT[i,,] = TT_K[i,,]%*%inv.D
U_K <- matrix(1, K, K); U_K[upper.tri(U_K)] <- -1; U <- U_K%*%inv.D[-1,]

###################
## Linear Effect ##
###################
C = 1000; zeros = rep(0,N)
dta <- list("N", "K", "tilde.t", "T", "TT", "U", "delta", "x", "zeros", "C")
ints <- function() {list(sigma.alpha=runif(1,0.01,100), alpha=rnorm(K+1), beta=rnorm(1))}
pars <- c("beta")

Linear.fit <- jags(data=dta,inits=ints,model.file="LinearMod.txt",parameters=pars,
n.chains=2,n.iter=12000,n.burnin=2000,n.thin=1)

###################
## Smooth Effect ##
###################
#Covariate Range Partition
tilde.x = 0:J/J

#Transformation/Penalty Matrix
OMEGA_beta<-rbind(c(1,rep(0,J-1)),cbind(0,abs(outer(tilde.x[-c(1,J+1)],tilde.x[-c(1,J+1)],"-"))^3))
svd_OMEGA_beta<-svd(OMEGA_beta)
inv.D_beta<-solve(t(svd_OMEGA_beta$v%*%(t(svd_OMEGA_beta$u)*sqrt(svd_OMEGA_beta$d))))

#Design Matrix
X = cbind(x-mean(x),t(sapply(x,function(z) abs(z - tilde.x[-c(1,J+1)])^3-abs(mean(x)-tilde.x[-c(1,J+1)])^3)))%*%inv.D_beta

#JAGS Setup
dta <- list("N", "K", "tilde.t", "T", "TT", "U", "delta", "J", "X", "zeros", "C")
ints <- function() {list(sigma.alpha=runif(1,0.01,100), alpha=rnorm(K+1),sigma.beta=runif(1,0.01,100), beta=rnorm(J))}

#Call JAGS
Smooth.fit <- jags(data=dta,inits=ints,model.file="SmoothMod.txt",parameters=pars,
n.chains=2,n.iter=22000,n.burnin=2000,n.thin=1)

#############################
## Shape-Restricted Effect ##
#############################
#Construct Inverse Transformation Matrix
L.inv = matrix(NA,ncol=J+1,nrow=J+1)
L.inv[1,] = c(1,rep(0,J))
L.inv[2,] = c(-1/2/tilde.x[2],1/2/tilde.x[2],rep(0,J-1))
for(j in 2:J){ L.inv[j+1,] = c(rep(0,j-2),1/(2*(tilde.x[j]-tilde.x[j-1])),-(tilde.x[j+1]-tilde.x[j-1])/(2*(tilde.x[j+1]-tilde.x[j])*(tilde.x[j]-tilde.x[j-1])),1/(2*(tilde.x[j+1]-tilde.x[j])),rep(0,J-j)) }
X2 = cbind(x,x**2,t(sapply(x,function(z) ifelse(z>tilde.x[-c(1,J+1)],1,0)*(z - tilde.x[-c(1,J+1)])**2)))%*%L.inv

#JAGS Setup
dta <- list("N", "K", "tilde.t", "T", "TT", "U", "delta", "J", "X2", "zeros", "C")
ints <- function() {list(sigma.alpha=runif(1,0.01,100), alpha=rnorm(K+1), beta.star=rexp(J+1))}

Restrict.fit <- jags(data=dta,inits=ints,model.file="RestrictMod.txt",parameters=pars,
n.chains=2,n.iter=22000,n.burnin=2000,n.thin=1)

#Write parameter estimates to files
write.table(rbind(round(Linear.fit$BUGSoutput$mean$beta,4)),paste("Linres",scen,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep=',',append=TRUE)
write.table(rbind(round(Smooth.fit$BUGSoutput$mean$beta,4)),paste("Smthres",scen,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep=',',append=TRUE)
write.table(rbind(round(Restrict.fit$BUGSoutput$mean$beta,4)),paste("Rstres",scen,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep=',',append=TRUE)
}

###############################
## Parallel Computation Code ##
###############################
library(snowfall)
sfInit(socketHosts=rep("localhost",8),cpus=8,type='SOCK',parallel=TRUE)
sfLibrary(R2jags)
sfExport('sim')
wrapper = function(x){sim(N=200,K=20,J=20,scen=1)}
res = sfLapply(1:200,wrapper)
sfStop()

#####################
## Compile Results ##
#####################
#Scenario
scen = 1

#True Baseline Hazard Function 
h = function(t) 4*(.5*cos(pi*t+pi/2) + 1.5/(t+1) - log(1.125))
H = function(t) 4*(.5*cos(pi*t)/pi + 1.5*log(t+1) - log(1.125)*t - 1/(2*pi))

#Flat Nonlinear Effect Scenario
if(scen==1) f = function(x) 5*(exp(20*(x-0.5))/(1+exp(20*(x-0.5)))-0.5)
#Smooth Nonlinear Effect Scenario
if(scen==2) f = function(x) 20*(x-.5)**3  
#Linear Effect Scenario
if(scen==3) f = function(x) 4*(x-.5)

#Read in parameter estimates for all simulation runs
Linres = as.matrix(read.csv(paste("Linres",scen,".txt",sep=""),header=FALSE))
Smthres = as.matrix(read.csv(paste("Smthres",scen,".txt",sep=""),header=FALSE))
Rstres = as.matrix(read.csv(paste("Rstres",scen,".txt",sep=""),header=FALSE))

#Covariate Range Partition (fixed for all simulation runs)
J = dim(Smthres)[2]; tilde.x = 0:J/J

#Grid to calculate estimates and RISE
x = seq(0,1,length.out=10000)

#Transformation/Penalty Matrix for Smooth Effect
OMEGA_beta<-rbind(c(1,rep(0,J-1)),cbind(0,abs(outer(tilde.x[-c(1,J+1)],tilde.x[-c(1,J+1)],"-"))^3))
svd_OMEGA_beta<-svd(OMEGA_beta)
inv.D_beta<-solve(t(svd_OMEGA_beta$v%*%(t(svd_OMEGA_beta$u)*sqrt(svd_OMEGA_beta$d))))

#Smooth Effect Design Matrix
X = cbind(x-0.5,t(sapply(x,function(x) abs(x - tilde.x[-c(1,J+1)])^3-abs(0.5-tilde.x[-c(1,J+1)])^3)))%*%inv.D_beta

#Transformation Matrix for Restricted Effect
L.inv = matrix(NA,ncol=J+1,nrow=J+1)
L.inv[1,] = c(1,rep(0,J))
L.inv[2,] = c(-1/2/tilde.x[2],1/2/tilde.x[2],rep(0,J-1))
for(j in 2:J){ L.inv[j+1,] = c(rep(0,j-2),1/(2*(tilde.x[j]-tilde.x[j-1])),-(tilde.x[j+1]-tilde.x[j-1])/(2*(tilde.x[j+1]-tilde.x[j])*(tilde.x[j]-tilde.x[j-1])),1/(2*(tilde.x[j+1]-tilde.x[j])),rep(0,J-j)) }

#Restricted Effect Design Matrix
X2 = cbind(x,x**2,t(sapply(x,function(z) ifelse(z>tilde.x[-c(1,J+1)],1,0)*(z - tilde.x[-c(1,J+1)])**2)))%*%L.inv
#Will need to subtract prediction at mean(x) = 0.5 
X2.bar = rbind(c(0.5,0.5**2,ifelse(0.5>tilde.x[-c(1,J+1)],1,0)*(0.5 - tilde.x[-c(1,J+1)])**2)%*%L.inv)

#HR Estimates over domain of x, relative to x = mean(x) (= 0.5)
Lin.lhr = t(apply(Linres,1,function(z) x*z-0.5*z))
Smth.lhr = t(X%*%t(Smthres))
Rst.lhr = t(apply(Rstres,1,function(x) X2%*%x-as.vector(X2.bar%*%x)))
#Average and Quantile Estimates
Lin.lhr.hat = apply(Lin.lhr,2,function(x) c(mean(x),quantile(x,c(0.025,0.975))))
Smth.lhr.hat = apply(Smth.lhr,2,function(x) c(mean(x),quantile(x,c(0.025,0.975))))
Rst.lhr.hat = apply(Rst.lhr,2,function(x) c(mean(x),quantile(x,c(0.025,0.975))))
#RISE of average estimate and best individual fit
RISE = c(sqrt(mean((Lin.lhr.hat[1,]-(f(x)-f(0.5)))**2)),sqrt(mean((Smth.lhr.hat[1,]-(f(x)-f(0.5)))**2)),sqrt(mean((Rst.lhr.hat[1,]-(f(x)-f(0.5)))**2)))
best.fit = c(which.min(apply(Lin.lhr,1,function(z) sqrt(mean((z-(f(x)-f(0.5)))**2)))),which.min(apply(Smth.lhr,1,function(z) sqrt(mean((z-(f(x)-f(0.5)))**2)))),which.min(apply(Rst.lhr,1,function(z) sqrt(mean((z-(f(x)-f(0.5)))**2)))))

#################################
## Generate panels in Figure 3 ##
#################################
pdf(paste("Lin",scen,".pdf",sep=""),family='serif')
par(mfrow=c(1,1),mar=c(2.7,3.2,.2,.2),las=1,mgp=c(1.5,.5,0),cex.lab=1.5,cex=2.5)
curve(f(x)-f(mean(x)),main=" ",xlab=expression(x),ylab=expression(f(x)-f(bar(x))),ylim=c(-3,3),col='white')
polygon(c(c(x[1],rep(x[-1],each=2)),rev(c(rep(x[-length(x)],each=2),x[length(x)]))),
 c(c(rep(Lin.lhr.hat[3,-length(x)],each=2),Lin.lhr.hat[3,length(x)]),c(rev(rep(Lin.lhr.hat[2,-1],each=2)),Lin.lhr.hat[2,1])),col='grey60',border='grey60')
lines(x,Lin.lhr.hat[1,],lwd=3)
lines(x,Lin.lhr[best.fit[1],],lwd=3,col='grey40')
curve(f(x)-f(mean(x)),from=0,to=1,add=TRUE,lwd=3,lty=2)
text(0.25, 2.9, paste("RISE =",formatC(round(RISE[1],2),format='f',digits=2)),cex=1.2)
if(scen==1) legend("bottomright",c("Truth","Average Estimate","Best Estimate"),lwd=3,col=c("black","black","grey40"),lty=c(2,1,1),bty='n',seg.len=1.8,cex=.9)
dev.off()

pdf(paste("Smth",scen,".pdf",sep=""),family='serif')
par(mfrow=c(1,1),mar=c(2.7,3.2,.2,.2),las=1,mgp=c(1.5,.5,0),cex.lab=1.5,cex=2.5)
curve(f(x)-f(mean(x)),main=" ",xlab=expression(x),ylab=expression(f(x)-f(bar(x))),ylim=c(-3,3),col='white')
polygon(c(c(x[1],rep(x[-1],each=2)),rev(c(rep(x[-length(x)],each=2),x[length(x)]))),
 c(c(rep(Smth.lhr.hat[3,-length(x)],each=2),Smth.lhr.hat[3,length(x)]),c(rev(rep(Smth.lhr.hat[2,-1],each=2)),Smth.lhr.hat[2,1])),col='grey60',border='grey60')
lines(x,Smth.lhr.hat[1,],lwd=3)
lines(x,Smth.lhr[best.fit[2],],lwd=3,col='grey40')
curve(f(x)-f(mean(x)),from=0,to=1,add=TRUE,lwd=3,lty=2)
text(0.25, 2.9, paste("RISE =",formatC(round(RISE[2],2),format='f',digits=2)),cex=1.2)
dev.off()

pdf(paste("Rst",scen,".pdf",sep=""),family='serif')
par(mfrow=c(1,1),mar=c(2.7,3.2,.2,.2),las=1,mgp=c(1.5,.5,0),cex.lab=1.5,cex=2.5)
curve(f(x)-f(mean(x)),main=" ",xlab=expression(x),ylab=expression(f(x)-f(bar(x))),ylim=c(-3,3),col='white')
polygon(c(c(x[1],rep(x[-1],each=2)),rev(c(rep(x[-length(x)],each=2),x[length(x)]))),
 c(c(rep(Rst.lhr.hat[3,-length(x)],each=2),Rst.lhr.hat[3,length(x)]),c(rev(rep(Rst.lhr.hat[2,-1],each=2)),Rst.lhr.hat[2,1])),col='grey60',border='grey60')
lines(x,Rst.lhr.hat[1,],lwd=3)
lines(x,Rst.lhr[best.fit[3],],lwd=3,col='grey40')
curve(f(x)-f(mean(x)),from=0,to=1,add=TRUE,lwd=3,lty=2)
text(0.25, 2.9, paste("RISE =",formatC(round(RISE[3],2),format='f',digits=2)),cex=1.2)
dev.off()