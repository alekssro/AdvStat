##--------------------------------------------------------
## Illustration of EM algorithm for mixture of two normals
##---------------------------------------------------------
## Simulate Data
alpha <- 0.25
nObs <- 100
x <- rbinom(n=nObs,size=1,alpha)
mu <- 5
y <- rnorm(n=nObs,mean=mu*x,sd=1)
## Illustrate data
s <- seq(-3,8,len=100)
hist(y,prob=TRUE,ylim=c(0,0.5))
lines(s,alpha*dnorm(s,mean=mu,sd=1),col="red",lwd=4,lty=1)
lines(s,(1-alpha)*dnorm(s,mean=0,sd=1),col="orange",lwd=4,lty=1)
lines(s,alpha*dnorm(s,mean=mu,sd=1)+(1-alpha)*dnorm(s,mean=0,sd=1),
      col="blue",lwd=2,lty=1)
##------------------------------------------------
## EM algorithm
##------------------------------------------------
## Starting values
alphaEst <- 0.1
muEst <- 3
threshold <- 0.001
## Number of iterations
nIter <- 20

# ## Run EM algorithm
# for (iter in 1:nIter){
#   xPrb <- alphaEst*dnorm(y,mean=muEst,sd=1)/
#     (alphaEst*dnorm(y,mean=muEst,sd=1)+
#      (1-alphaEst)*dnorm(y,mean=0,sd=1))
#   alphaEst <- sum(xPrb)/nObs
#   muEst <- sum(xPrb*y)/sum(xPrb)
#   logLk <- sum(log(alphaEst*dnorm(y,mean=muEst,sd=1)+
#                    (1-alphaEst)*dnorm(y,mean=0,sd=1)))
#   cat("Iteration:",iter,
#       "alpha:",alphaEst,"mu:",muEst,
#       "logLk",logLk,"\n")
# }
# ##--------------------------------------------------------------
# ## Illustrate final fit
# ##--------------------------------------------------------------
# lines(s,alphaEst*dnorm(s,mean=muEst,sd=1),col="red",lwd=4,lty=2)
# lines(s,(1-alphaEst)*dnorm(s,mean=0,sd=1),col="orange",lwd=4,lty=2)
# lines(s,alphaEst*dnorm(s,mean=muEst,sd=1)+(1-alphaEst)*dnorm(s,mean=0,sd=1),
#       col="blue",lwd=2,lty=2)
# ## Add legend to plot
# legend("topright",
#        c("True mixture distribution","Estimated mixture distribution"),
#        bty="n",col=c("black","black"),lty=c(1,2),lwd=3)

datA <- scan(file = "MixtNormA(1).dat", sep = "")
datB <- scan(file = "MixtNormB(1).dat", sep = "")

## Illustrate data
s <- seq(-3,8,len=length(datA))
hist(datA,prob=TRUE,ylim=c(0,0.5))
hist(datB,prob=TRUE,ylim=c(0,0.5), breaks = 15)

source("EM_function.R")
myEMfunction(datB, alphaEst = 0.5, muEst = 5, maxIter = 500, stopThres = 1e-100)

