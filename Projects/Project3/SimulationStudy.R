##-------------------------------------
## NMF: Setting up a simulation study 
##-------------------------------------
## Simulate mutational signatures
##----------------------------------------
source("NQPAlgorithms.R")

## Assume K=4 signatures andd N=100 mutational types
H <- matrix(runif(4*100),nrow=4,ncol=100)
## Signature k has an overrepresentation of 
## mutation types (1+(k-1)*25):(k*25), i.e.
## 1st signature mutation type 1:25,
## 2nd 26:50, 3rd 51:75 and 4th 76:100
for (k in 1:4){ 
  mut.type <- (1+(k-1)*25):(k*25) 
  H[k,mut.type] <- H[k,mut.type]+1
  ## Normalize mutation types
  H[k,] <- H[k,]/sum(H[k,])
}
## Plot the signatures
col.muta <- c(rep("blue",25),rep("red",25),rep("green",25),rep("pink",25))
par(mfrow=c(4,1))
plot(H[1,],type="h",col=col.muta,lwd=2,xlab="",ylab="",main="Signature 1",ylim=c(0,0.03))
plot(H[2,],type="h",col=col.muta,lwd=2,xlab="",ylab="",main="Signature 2",ylim=c(0,0.03))
plot(H[3,],type="h",col=col.muta,lwd=2,xlab="",ylab="",main="Signature 3",ylim=c(0,0.03))
plot(H[4,],type="h",col=col.muta,lwd=2,xlab="",ylab="",main="Signature 4",ylim=c(0,0.03))
##-------------------
## Define loadings
##-------------------
## Here we just choose the loadings 
w <- 100000*c(0.5,0.3,0.25,0.05)
## Mean of the N=100 mutational types
mn.v <- as.vector(w%*%H)
##-------------------------------
## Normal model: Simulate data 
##-------------------------------
par(mfrow=c(1,1))
plot(mn.v,type="h",col=col.muta,lwd=1,
     xlab="mutation type",ylab="mutation count",main="Mutation counts")
obs.v <- round(abs(rnorm(100,mean=mn.v,sd=sqrt(mean(mn.v)))))
points(1:100+0.4,obs.v+0.5,type="h",col="black",lwd=1)
## Residual plot
raw.resi <- obs.v-mn.v
srt <- sort(mn.v,index=TRUE)
plot(srt$x,raw.resi[srt$ix],pch=19,cex=0.5,xlab="expected",ylab="residual")
a.coef <- mean(abs(raw.resi[srt$ix]))
abline(a=a.coef,b=0,col="red")
abline(a=-a.coef,b=0,col="red")
##-------------------------------
## Poisson model: Simulate data
##-------------------------------
par(mfrow=c(1,1))
plot(mn.v,type="h",col=col.muta,lwd=1,ylim=c(0,max(mn.v)),
     xlab="mutation type",ylab="mutation count",main="Mutation counts")
obs.v <- rpois(100,mn.v)
points(1:100+0.4,obs.v+0.5,type="h",col="black",lwd=1)
## Residual plot
raw.resi <- obs.v-mn.v
srt <- sort(mn.v,index=TRUE)
plot(srt$x,raw.resi[srt$ix],pch=19,cex=0.5,xlab="expected",ylab="residual")
lm.coef <- lm(abs(raw.resi[srt$ix])~srt$x+0)$coefficient
abline(a=0,b=lm.coef,col="red")
abline(a=0,b=-lm.coef,col="red")
##--------------------------------------------------
## MM algorithm for loadings based on observed data
##--------------------------------------------------
A <- 2*(H%*%t(H))
## Check on true mean 
b <- as.vector( -2*mn.v%*%t(H) )
as.vector(MajorizeMinimizeNQP(A,b)$x)
## Apply on observed data
b <- as.vector( -2*obs.v%*%t(H) )
w.est <- as.vector(MajorizeMinimizeNQP(A,b)$x)
## Compare objection function for estimated and true values
sum( (obs.v-w.est%*%H)^2 )
sum( (obs.v-w%*%H)^2 )
plot(obs.v,w%*%H,pch=19,col="red",cex=0.8,
     xlab="observed",ylab="expected")
points(obs.v,w.est%*%H,pch=19,col="blue",cex=0.5)
##-------------------------------------------------
