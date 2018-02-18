setwd("C:/Users/TMurray/Documents/PwL_Surv_Model/RCode/SimResults/SimRes4.1")
library(R2jags)

###################
## JAGS Programs ##
###################
#Piecewise Exponential Model
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
		 mu[i,k] <- theta[i,k]*exp(gamma[k])
		 zed[i,k] ~ dpois(mu[i,k]);
		}
	} 
	# Prior on hazard parameters
		gamma[1] ~ dnorm(0,0.0001)
		for(k in 2:K){
			gamma[k] ~ dnorm(gamma[k-1],tau.gamma)
	 	}
	# Hyperpriors on historical variance parameters
		 tau.gamma <- 1/(sigma.gamma*sigma.gamma)
		 sigma.gamma ~ dunif(0.01,100)
}","PEMod.txt")

#Piecewise Linear log-Hazard Model
write("model{
	for (i in 1:N){
		 lhaz[i] <- inprod(alpha[],T[i,])
		for(k in 1:K){
		 Haz[i,k] <-  exp(inprod(alpha[],TT[i,k,]))*(1-exp(-(TT[i,k,2]-tilde.t[k])*inprod(alpha[],U[k,])))/inprod(alpha[],U[k,])
		}
		 neg.LL[i] <- -delta[i]*lhaz[i] + sum(Haz[i,]) + C
		 zeros[i] ~ dpois(neg.LL[i]);
	}
	# Prior on hazard parameters
		 alpha[1] ~ dnorm(0,0.0001)
		 alpha[2] ~ dnorm(0,0.0001)
		for(k in 2:K){
			alpha[k+1] ~ dnorm(0,tau.alpha)
	 	}
	# Hyperpriors on historical variance parameters
		 tau.alpha <- 1/(sigma.alpha*sigma.alpha)
		 sigma.alpha ~ dunif(0.01,100)
}","PLMod.txt")

#########################
## Simulation Function ##
#########################
#Simulation Function
sim = function(N,K,scen){
 #Inverse Function
 inverse = function (f, lower = -100, upper = 100){
   function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
 }
 if(scen==1){
  ##Scenario 1 (Valley)
  h = function(x) 4*(.5*cos(pi*x+pi/2) + 1.5/(x+1) - log(1.125))
  H = function(x) 4*(.5*cos(pi*x)/pi + 1.5*log(x+1) - log(1.125)*x - 1/(2*pi))
  H.inv = inverse(function (x) ifelse(x < 1, H(x), 2*(x-1) + H(1)) , 0, 100)
 }
 if(scen==2){
  ##Scenario 2 (Peak)
  h = function(x) ifelse(x<.5,exp(8*x)/5,exp(-8*(x-1))/5)
  H = function(x) ifelse(x<.5,exp(8*x)/40-1/40,exp(4)/20-1/40-exp(-8*(x-1))/40)
  H.inv = inverse(function (x) ifelse(x<1, H(x), 2*(x-1) + H(1)) , 0, 100)
 }
 if(scen==3){
  ##Scenario 3 (Double-Hill)
  h = function(x) 1.4*(sin(pi*3*x) + 1.5)
  H = function(x) 1.4*(-cos(pi*3*x)/(3*pi) + 1/(3*pi) + 1.5*x)
  H.inv = inverse(function (x) H(x) , 0, 100)
 }

 #Generate Data from survival distribtuion having true hazard above
 u = runif(N)
 y = as.vector(unlist(sapply(-log(u),function(z) H.inv(z))))
 c = runif(N,0,5)
 t = pmin(y,c,rep(1,N))
 delta = ifelse(t==y,1,0)
 tilde.t = c(0:K/K)

 ############
 # PE Model #
 ############
 #JAGS Setup
 dta <- list("N", "t", "delta", "K", "tilde.t")
 ints <- function() {list(sigma.gamma=runif(1,0.01,100), gamma=rnorm(K))}
 pars <- c("gamma") 
 
 #Call JAGS
 PE.fit <- jags(data=dta,inits=ints,model.file="PEMod.txt",parameters=pars,
 n.chains=2,n.iter=12000,n.burnin=2000,n.thin=1)

 ############
 # PL Model #
 ############
 #JAGS Setup
 zeros = rep(0,N); C = 1000
 dta <- list("N", "K", "tilde.t", "T", "TT", "U", "delta", "zeros", "C")
 ints <- function() {list(sigma.alpha=runif(1,0.01,100), alpha=rnorm(K+1))}
 pars <- c("alpha")

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
 U_K <- matrix(1, K, K)
 U_K[upper.tri(U_K)] <- -1
 U <- U_K%*%inv.D[-1,]

 #Call JAGS
 PL.fit <- jags(data=dta,inits=ints,model.file="PLMod.txt",parameters=pars,
 n.chains=2,n.iter=12000,n.burnin=2000,n.thin=1)

 #Write mean parameter estimates to a results file
 write.table(rbind(round(PE.fit$BUGSoutput$mean$gamma,4)),paste("PEres",scen,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep=',',append=TRUE)
 write.table(rbind(round(PL.fit$BUGSoutput$mean$alpha,4)),paste("PLres",scen,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep=',',append=TRUE)
}

###############################
## Parallel Computation Code ##
###############################
library(snowfall)
sfInit(socketHosts=rep("localhost",8),cpus=8,type='SOCK',parallel=TRUE)
sfLibrary(MASS);sfLibrary(survival);sfLibrary(R2jags)
sfExport('sim')
wrapper = function(x){sim(N=200,K=20,scen=3)}
res = sfLapply(1:200,wrapper)
sfStop()

#####################
## Compile Results ##
#####################
#Scenario
scen=1
if(scen==1){
 ##Scenario 1 (Valley)
 h = function(x) 4*(.5*cos(pi*x+pi/2) + 1.5/(x+1) - log(1.125))
 H = function(x) 4*(.5*cos(pi*x)/pi + 1.5*log(x+1) - log(1.125)*x - 1/(2*pi))
}
if(scen==2){
 ##Scenario 2 (Peak)
 h = function(x) ifelse(x<.5,exp(8*x)/5,exp(-8*(x-1))/5)
 H = function(x) ifelse(x<.5,exp(8*x)/40-1/40,exp(4)/20-1/40-exp(-8*(x-1))/40)
}
if(scen==3){
 ##Scenario 3 (Double-Hill)
 h = function(x) 1.4*(sin(pi*3*x) + 1.5)
 H = function(x) 1.4*(-cos(pi*3*x)/(3*pi) + 1/(3*pi) + 1.5*x)
}

#Read in mean parameter estimates for all simulation runs
PEres = as.matrix(read.csv(paste("PEres",scen,".txt",sep=""),header=FALSE))
PLres = as.matrix(read.csv(paste("PLres",scen,".txt",sep=""),header=FALSE))

#Time Partition
K = dim(PEres)[2]; tilde.t = 0:K/K

#log-hazard estimates
 t = seq(0,1,length.out=10000)
 #PE models
 S = t(sapply(t,function(t) ifelse(t<c(tilde.t[-c(1,K+1)],1.1),1,0)*ifelse(t>=tilde.t[-(K+1)],1,0)))
 PE.lhaz = t(S%*%t(PEres))
 #PL models
 #Time Transformation Matrix
 OMEGA_alpha <- abs(outer(tilde.t[-c(1,K+1)],tilde.t[-c(1,K+1)],"-"))
 svd.OMEGA_alpha<-svd(OMEGA_alpha)
 sqrt.OMEGA_alpha<-t(svd.OMEGA_alpha$v%*%(t(svd.OMEGA_alpha$u)*sqrt(svd.OMEGA_alpha$d)))
 inv.D <- solve(cbind(c(1,rep(0,K)),c(0,1,rep(0,K-1)),rbind(rep(0,K-1),rep(0,K-1),sqrt.OMEGA_alpha)))
 T_K = cbind(1,t(sapply(t,function(z) abs(z - tilde.t[-c(K+1)]) - tilde.t[-c(K+1)])))
 T = T_K%*%inv.D
 PL.lhaz = t(T%*%t(PLres))

#Survival estimates
 #PE models
 SS = t(sapply(t,function(t) (pmin(t, tilde.t[-1]) - tilde.t[-(K+1)])*ifelse(t > tilde.t[-(K+1)],1,0)))
 PE.Surv = t(exp(-SS%*%exp(t(PEres))))
 #PL models
 sk = t(sapply(t,function(z) pmax(pmin(z,tilde.t[-1]),tilde.t[-(K+1)])))
 TT_K = TT = array(NA,c(length(t),K,K+1))
 for(k in 1:K) TT_K[,k,] = cbind(1,t(sapply(sk[,k],function(z) abs(z - tilde.t[-c(K+1)]) - tilde.t[-c(K+1)])))
 for(i in 1:length(t)) TT[i,,] = TT_K[i,,]%*%inv.D
 U_K <- matrix(1, K, K); U_K[upper.tri(U_K)] <- -1; U <- U_K%*%inv.D[-1,]
 PL.Surv = t(sapply(1:dim(PLres)[1],function(r) sapply(1:length(t),function(i) exp(-sum(exp(TT[i,,]%*%cbind(PLres[r,]))*(1-exp(-(TT[i,,2]-tilde.t[-(K+1)])*(U%*%cbind(PLres[r,]))))/(U%*%cbind(PLres[r,])))))))

#Average and Quantile Estimates
PE.lhaz.hat = apply(PE.lhaz,2,function(x) c(mean(x),quantile(x,c(0.025,0.975))))
PL.lhaz.hat = apply(PL.lhaz,2,function(x) c(mean(x),quantile(x,c(0.025,0.975))))
PE.Surv.hat = apply(PE.Surv,2,function(x) c(mean(x),quantile(x,c(0.025,0.975))))
PL.Surv.hat = apply(PL.Surv,2,function(x) c(mean(x),quantile(x,c(0.025,0.975))))

#RISE of average estimate and best individual fit
RISE = c(sqrt(mean((PE.lhaz.hat[1,]-log(h(t)))**2)),sqrt(mean((PL.lhaz.hat[1,]-log(h(t)))**2)),sqrt(mean((PE.Surv.hat[1,]-exp(-H(t)))**2)),sqrt(mean((PL.Surv.hat[1,]-exp(-H(t)))**2)))
best.fit = c(which.min(apply(PE.lhaz,1,function(x) sqrt(mean((x-log(h(t)))**2)))),which.min(apply(PL.lhaz,1,function(x) sqrt(mean((x-log(h(t)))**2)))))

###############################
# Generate panels in Figure 1 #
###############################
#Piecewise Exponential log-hazard results
pdf(paste("PElhaz",scen,".pdf",sep=""),family="serif")
par(mfrow=c(1,1),mar=c(2.7,3.2,.2,.2),las=1,mgp=c(1.5,.5,0),cex.lab=1.5,cex=2.5)
curve(log(h(x)),main=" ",xlab=expression(t),ylab=expression(paste("log{h(t)}")),ylim=c(-2.5,2.5),col='white')
polygon(c(c(t[1],rep(t[-1],each=2)),rev(c(rep(t[-length(t)],each=2),t[length(t)]))),
 c(c(rep(PE.lhaz.hat[3,-length(t)],each=2),PE.lhaz.hat[3,length(t)]),c(rev(rep(PE.lhaz.hat[2,-1],each=2)),PE.lhaz.hat[2,1])),col='grey60',border='grey60')
lines(t,PE.lhaz.hat[1,],lwd=3,type='s')
lines(t,PE.lhaz[best.fit[1],],lwd=3,type='s',col='grey40')
curve(log(h(x)),from=0,to=1,add=TRUE,lwd=3,lty=2)
if(scen%in%c(1,3)) text(.7, 2.4, paste("RISE =",formatC(round(RISE[1],2),format='f',digits=2)),cex=1.5)
if(scen==2) text(.5, -2.4, paste("RISE =",formatC(round(RISE[1],2),format='f',digits=2)),cex=1.5)
if(scen==1) legend("bottomleft",c("Truth","Average Estimate","Best Estimate"),lwd=3,col=c("black","black","grey40"),lty=c(2,1,1),bty='n',cex=1.2,seg.len=1.2)
dev.off()

#Piecewise Linear log-hazard results
pdf(paste("PLlhaz",scen,".pdf",sep=""),family="serif")
par(mfrow=c(1,1),mar=c(2.7,3.2,.2,.2),las=1,mgp=c(1.5,.5,0),cex.lab=1.5,cex=2.5)
curve(log(h(x)),main=" ",xlab=expression(t),ylab=expression(paste("log{h(t)}")),ylim=c(-2.5,2.5),col='white')
polygon(c(c(t[1],rep(t[-1],each=2)),rev(c(rep(t[-length(t)],each=2),t[length(t)]))),
 c(c(rep(PL.lhaz.hat[3,-length(t)],each=2),PL.lhaz.hat[3,length(t)]),c(rev(rep(PL.lhaz.hat[2,-1],each=2)),PL.lhaz.hat[2,1])),col='grey60',border='grey60')
lines(t,PL.lhaz.hat[1,],lwd=3)
lines(t,PL.lhaz[best.fit[2],],lwd=3,col='grey40')
curve(log(h(x)),from=0,to=1,add=TRUE,lwd=3,lty=2)
if(scen%in%c(1,3)) text(.7, 2.4, paste("RISE =",formatC(round(RISE[2],2),format='f',digits=2)),cex=1.5)
if(scen==2) text(.5, -2.4, paste("RISE =",formatC(round(RISE[2],2),format='f',digits=2)),cex=1.5)
dev.off()

#Piecewise Exponential survival results
pdf(paste("PESurv",scen,".pdf",sep=""),family="serif")
par(mfrow=c(1,1),mar=c(2.7,3.2,.2,.2),las=1,mgp=c(1.5,.5,0),cex.lab=1.5,cex=2.5)
curve(exp(-H(x)),main=" ",xlab=expression(t),ylab=expression(S(t)),ylim=c(0,1),col='white')
polygon(c(c(t[1],rep(t[-1],each=2)),rev(c(rep(t[-length(t)],each=2),t[length(t)]))),
 c(c(rep(PE.Surv.hat[3,-length(t)],each=2),PE.Surv.hat[3,length(t)]),c(rev(rep(PE.Surv.hat[2,-1],each=2)),PE.Surv.hat[2,1])),col='grey60',border='grey60')
lines(t,PE.Surv.hat[1,],lwd=3)
lines(t,PE.Surv[best.fit[1],],lwd=3,col='grey40')
curve(exp(-H(x)),from=0,to=1,add=TRUE,lwd=3,lty=2)
text(.59, 0.99, paste("100xRISE =",formatC(round(100*RISE[3],2),format='f',digits=2)),cex=1.5)
dev.off()

#Piecewise Linear survival results
pdf(paste("PLSurv",scen,".pdf",sep=""),family="serif")
par(mfrow=c(1,1),mar=c(2.7,3.2,.2,.2),las=1,mgp=c(1.5,.5,0),cex.lab=1.5,cex=2.5)
curve(exp(-H(x)),main=" ",xlab=expression(t),ylab=expression(S(t)),ylim=c(0,1),col='white')
polygon(c(c(t[1],rep(t[-1],each=2)),rev(c(rep(t[-length(t)],each=2),t[length(t)]))),
 c(c(rep(PL.Surv.hat[3,-length(t)],each=2),PL.Surv.hat[3,length(t)]),c(rev(rep(PL.Surv.hat[2,-1],each=2)),PL.Surv.hat[2,1])),col='grey60',border='grey60')
lines(t,PL.Surv.hat[1,],lwd=3)
lines(t,PL.Surv[best.fit[2],],lwd=3,col='grey40')
curve(exp(-H(x)),from=0,to=1,add=TRUE,lwd=3,lty=2)
text(.59, 0.99, paste("100xRISE =",formatC(round(100*RISE[4],2),format='f',digits=2)),cex=1.5)
dev.off()