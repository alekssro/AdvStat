##------------------------------------------------------------
## Make sure all the basic algorithms are available
##-----------------------------------------------------------
source("NQPAlgorithms.R")
##------------------------------------------------------------
## Lange, Chi and Zhou (2013) test problem
## Small change: 
## Use exponential with rate 1 instead of standard normal
## in the design matrix to ensure positive entries.
##------------------------------------------------------------
## Simulate the data
##------------------------------------------------------------
set.seed(10)
N = 100 ; K = 50
X = matrix(rexp(N*K,rate=1),nrow=N,ncol=K)
tht = runif(K)
err = rnorm(N,mean=0,sd=1)
Ytrue = X%*%tht
Y = Ytrue + err
plot(Y,pch=19,cex=0.8,ylim=c(0,max(c(Y,Ytrue))))
points(Ytrue,pch=19,cex=0.7,col="red")
for (i in 1:N) {
  lines(c(i,i),c(Y[i],Ytrue[i]),col="black")
  }
##------------------------------------------------------------
## Unconstrained solution
##------------------------------------------------------------
thtUncHat <- solve(t(X)%*%X)%*%t(X)%*%Y
plot(1:K,tht,pch=19,ylim=range(thtUncHat),cex=0.8)
abline(h=0,col="blue")
points(1:K,thtUncHat,pch=19,col="blue",cex=0.8)
for (i in 1:K) {
  lines(c(i,i),c(tht[i],thtUncHat[i]),col="blue")
}
##---------------------------------------------------------
## Objective function for original, 
## constrained and unconstrained theta
##---------------------------------------------------------
sum( (Y-X%*%tht)^2 )
sum( (Y-X%*%thtUncHat)^2 )
A = 2*t(X)%*%X
b = as.vector(-2*t(X)%*%Y)
##--------------------------------------------------------
## Objective function 1/2*t(x)%*%A%*%x+t(b)%*%x
## Plot RSS=Objective function+t(Y)Y
##---------------------------------------------------------
n.rep <- 100 ; tol <- 0.0001 ; mxIter <- 1100
##---------------------------------------------------------
## MM algorithm
## Note: The algorithms are so fast that we need to replicate
## the experiments in order to obtain stable times.
##-------------------------------------------------------
MM <- MajorizeMinimizeNQP(A=A,b=b,init=tht,maxIter=mxIter,tol=tol)
MM$time <- rowMeans(
  replicate(n.rep,MajorizeMinimizeNQP(A=A,b=b,init=tht,maxIter=mxIter,tol=tol)$time))
MM.RSS <- MM$objfctList+sum(Y^2)
MM.RSS.min <- sum( (Y-X%*%MM$x)^2 )
MM.RelLoss <- (MM.RSS-MM.RSS.min)/MM.RSS.min
plot(MM.RSS,type="l",col="black",xlab="Iteration",ylab="RSS")
plot(MM$time,MM.RSS,type='l',col='black',lwd=1,xlab='Time',ylab='RSS')
abline(h=MM.RSS.min)
##-------------------------------
## Projected Coordinate Descent
##-------------------------------
PCD <- CoordinateDescentNQP(A=A,b=b,init=tht,maxIter=mxIter,tol=tol)
PCD$time <- rowMeans(
  replicate(n.rep,CoordinateDescentNQP(A=A,b=b,init=tht,maxIter=mxIter,tol=tol)$time))
PCD.RSS <- PCD$objfctList+sum(Y^2)
PCD.RSS.min <- min(PCD.RSS)
PCD.RelLoss <- (PCD.RSS-PCD.RSS.min)/PCD.RSS.min
plot(PCD.RSS,type="l",col="black",xlab="Iteration",ylab="RSS")
plot(PCD$time,PCD.RSS,type="l",col="black",lwd=1,xlab="Time",ylab="RSS")
abline(h=PCD.RSS.min)
##--------------------------------------
## Projected Gradient Descent
##-------------------------------------
## Exact line search
PGDe <- ProjectedGradientDescentNQP(A=A,b=b,init=tht,maxIter=mxIter,
                                    tol=tol,exactLineSearch=T)
PGDe$time <- rowMeans(
  replicate(n.rep,ProjectedGradientDescentNQP(A=A,b=b,init=tht,maxIter=mxIter,tol=tol,exactLineSearch=T)$time))
PGDe.RSS <- PGDe$objfctList+sum(Y^2)
PGDe.RSS.min <- min(PGDe.RSS)
PGDe.RelLoss <- (PGDe.RSS-PGDe.RSS.min)/PGDe.RSS.min
plot(PGDe.RSS,type="l",col="black",xlab="Iteration",ylab="RSS")
plot(PGDe$time,PGDe.RSS,type="l",col="black",xlab="Time",ylab="RSS")
## r=0.5 (small)
PGDs <- ProjectedGradientDescentNQP(A=A,b=b,init=tht,maxIter=mxIter,
                                    tol=tol,exactLineSearch=F,r=0.5)
PGDs$time <- rowMeans(
  replicate(n.rep,ProjectedGradientDescentNQP(A=A,b=b,init=tht,maxIter=mxIter,tol=tol,exactLineSearch=F,r=0.5)$time))
PGDs.RSS <- PGDs$objfctList+sum(Y^2)
PGDs.RSS.min <- min(PGDs.RSS)
PGDs.RelLoss <- (PGDs.RSS-PGDs.RSS.min)/PGDs.RSS.min
plot(PGDs.RSS,type="l",col="black",xlab="Iteration",ylab="RSS")
plot(PGDs$time,PGDs.RSS,type="l",col="black",xlab="Time",ylab="RSS")
## r=2 (large)
PGDl <- ProjectedGradientDescentNQP(A=A,b=b,init=tht,tol=tol,maxIter=mxIter,
                                    exactLineSearch=F,r=2)
PGDl$time <- rowMeans(
  replicate(n.rep,ProjectedGradientDescentNQP(A=A,b=b,init=tht,tol=tol,maxIter=mxIter,exactLineSearch=F,r=2)$time))
PGDl.RSS <- PGDl$objfctList+sum(Y^2)
PGDl.RSS.min <- min(PGDl.RSS)
PGDl.RelLoss <- (PGDl.RSS-PGDl.RSS.min)/PGDl.RSS.min
plot(PGDl.RelLoss,type="l",col="black",xlab="Iteration",ylab="RSS")
plot(PGDl$time,PGDl.RSS,type="l",col="black",xlab="Time",ylab="RSS")
##-----------------------
## Multiplicative EM
##-----------------------
MuEM <- MultEMNQP(A=A,b=b,init=tht,maxIter=mxIter,tol=tol)
MuEM$time <- rowMeans(
  replicate(n.rep,MultEMNQP(A=A,b=b,init=tht,maxIter=mxIter,tol=tol)$time))
MuEM.RSS <- MuEM$objfctList+sum(Y^2)
MuEM.RSS.min <- min(MuEM.RSS)
MuEM.RelLoss <- (MuEM.RSS-min(MuEM.RSS))/min(MuEM.RSS)
plot(MuEM.RSS,type="l",col="black",xlab="Iteration",ylab="RSS")
plot(MuEM$time,MuEM.RSS,type="l",col="black",xlab="Time",ylab="RSS")
##-----------------------
## Fevotte-Cemgil EM
##-----------------------
FCEM <- FevotteCemgilNQP(A=A,b=b,init=tht,maxIter=mxIter,tol=tol)
FCEM$time <- rowMeans(
  replicate(n.rep,FevotteCemgilNQP(A=A,b=b,init=tht,maxIter=mxIter,tol=tol)$time))
FCEM.RSS <- FCEM$objfctList+sum(Y^2)
FCEM.RSS.min <- min(FCEM.RSS)
FCEM.RelLoss <- (FCEM.RSS-FCEM.RSS.min)/FCEM.RSS.min
plot(FCEM.RSS,type="l",col="black",xlab="Iteration",ylab="RSS")
plot(FCEM$time,FCEM.RSS,type="l",col="black",xlab="Time",ylab="RSS")
##-----------------------
## Cone projection
##-----------------------
CNP <- ConeProjNonNegRegr(Y,X)
CNP$time <- mean(replicate(n.rep,ConeProjNonNegRegr(Y,X)$time))
##-------------------------------------------------------------------
## Figures for illustration
##-------------------------------------------------------------------
## Algorithms
algs <- c("MM","PCD","PGD r=0.5","PGD r=2.0","PGD Exact",
          "MuEM","FC-EM","CNP")
## Corresponding colors
cols <- c("red","orange","lightblue","blue","darkblue",
          "lightgreen","darkgreen","purple")
##---------------------------------------------------------------------
## Plot RSS versus iteration: Most simple plot
##---------------------------------------------------------------------
#pdf("TestProblemRSSvsIter.pdf",width=7.0,height=7.0)
plot(MM.RSS,log="x",
     xlim=c(1,1200),ylim=c(90,130),type="l",lwd=2,col=cols[1],
     xlab="Iteration",ylab="RSS",cex=1.5,cex.lab=1.4,cex.axis=1.3)
points(PCD.RSS,type="l",lwd=2,col=cols[2])
points(PGDs.RSS,type="l",lwd=5,col=cols[3])
points(PGDl.RSS,type="l",lwd=2,col=cols[4])
points(PGDe.RSS,type="l",lwd=2,col=cols[5])
points(MuEM.RSS,type="l",lwd=2,col=cols[6])
lines(FCEM.RSS,lty=2,lwd=2,col=cols[7])
points(1,CNP$objfct,col=cols[8],pch=4,lwd=2,cex=1.5)
legend('topright',algs,col=cols,
       lty=c(rep(1,6),2,NA),pch=c(rep(NA,7),4),pt.cex=1.5,
       lwd=c(2,2,5,rep(2,5)),bty='n')
abline(h=CNP$objfct,col=cols[8],lwd=1,lty=2)
#dev.off()
##---------------------------------------------------------------------
## Plot RSS versus iteration 
##---------------------------------------------------------------------
#pdf("TestProblemRelRSSvsIter.pdf",width=7.0,height=7.0)
plot((MM.RSS-CNP$objfct)/CNP$objfct,log="y",
     xlim=c(0,1200),ylim=c(1e-4,1),type="l",lwd=3,col=cols[1],
     xlab="Iteration",ylab="Relative RSS",cex=1.5,cex.lab=1.4,cex.axis=1.3)
points((PCD.RSS-CNP$objfct)/CNP$objfct,type="l",lwd=2,col=cols[2])
points((PGDs.RSS-CNP$objfct)/CNP$objfct,type="l",lwd=2,col=cols[3])
points((PGDl.RSS-CNP$objfct)/CNP$objfct,type="l",lwd=2,col=cols[4])
points((PGDe.RSS-CNP$objfct)/CNP$objfct,type="l",lwd=1,col=cols[5])
points((MuEM.RSS-CNP$objfct)/CNP$objfct,type="l",lwd=2,col=cols[6])
points((FCEM.RSS-CNP$objfct)/CNP$objfct,type="l",lwd=1,col=cols[7])
legend('topright',algs[1:7],col=cols[1:7],lty=1,lwd=2,bty='n')
#dev.off()
##-------------------------------------------------------------------
## Zoom in: Demonstrate that stopping criteria is non-trivial 
##-------------------------------------------------------------------
#pdf("TestProblemRSSvsIterZoom.pdf",width=7.0,height=7.0)
tole <- 0.001 
MMindx <- min( which( (diff(-MM$objfctList)<tole)==TRUE ) ) 
plot(MM.RSS[1:MMindx],xlim=c(0,1200),ylim=c(90,92),type="l",lwd=3,col=cols[1],
     xlab="Iteration",ylab="RSS",cex=1.5,cex.lab=1.4,cex.axis=1.3)
PCDindx <- min( which( (diff(-PCD$objfctList)<tole)==TRUE ) )
points(PCD.RSS[1:PCDindx],type="l",lwd=2,col=cols[2])
PGDsindx <- min( which( (diff(-PGDs$objfctList)<tole)==TRUE ) )
points(PGDs.RSS[1:PGDsindx],type="l",lwd=2,col=cols[3])
PGDlindx <- min( which( (diff(-PGDl$objfctList)<tole)==TRUE ) )
points(PGDl.RSS[1:PGDlindx],type="l",lwd=2,col=cols[4])
PGDeindx <- length(PGDe$objfctList)
points(PGDe.RSS[1:PGDeindx],type="l",lwd=1,col=cols[5])
MuEMindx <- length(MuEM$objfctList)
points(MuEM.RSS[1:MuEMindx],type="l",lwd=2,col=cols[6])
FCEMindx <- min( which( (diff(-FCEM$objfctList)<tole)==TRUE ) )
lines(FCEM.RSS[1:FCEMindx],lty=2,lwd=2,col=cols[7])
points(1,CNP$objfct,col=cols[8],pch=4,lwd=2,cex=1.5)
abline(h=CNP$objfct,lty=2,col=cols[8],lwd=2)
legend('topright',algs,col=cols,
       lty=c(rep(1,6),2,NA),pch=c(rep(NA,7),4),pt.cex=1.5,
       lwd=2,bty='n')
#dev.off()
##--------------------------------------------------------------------
## Plot loss versus time 
##--------------------------------------------------------------------
#pdf("TestProblemRSSvsTime.pdf",width=7.0,height=7.0)
plot(MM$time[1:MMindx],MM.RSS[1:MMindx],ylim=c(90,92),
     xlim=c(0,max(MuEM$time+0.01)),
     type="l",lwd=3,col=cols[1],
     xlab="Time",ylab="RSS",cex=1.5,cex.lab=1.4,cex.axis=1.3)
points(PGDs$time[1:PGDsindx],PGDs.RSS[1:PGDsindx],type="l",lwd=2,col=cols[3])
points(PGDl$time[1:PGDlindx],PGDl.RSS[1:PGDlindx],type="l",lwd=2,col=cols[4])
points(PGDe$time[1:PGDeindx],PGDe.RSS[1:PGDeindx],type="l",lwd=2,col=cols[5])
points(PCD$time[1:PCDindx],PCD.RSS[1:PCDindx],type="l",lwd=2,col=cols[2])
points(MuEM$time[1:MuEMindx],MuEM.RSS[1:MuEMindx],type="l",lwd=2,col=cols[6])
points(FCEM$time[1:FCEMindx],FCEM.RSS[1:FCEMindx],type="l",lwd=2,col=cols[7])
points(CNP$time,CNP$objfct,col=cols[8],pch=4,lwd=2,cex=1.5)
legend('topright',algs,col=cols,
       lty=c(rep(1,7),NA),pch=c(rep(NA,7),4),pt.cex=1.5,
       lwd=2,bty='n')
abline(h=CNP$objfct,lty=2,col=cols[8],lwd=2)
##--------------------------------------------------------------------
## Stopping criteria as a function of time
##--------------------------------------------------------------------
#pdf("TestProblemRSSDiffvsTime.pdf",width=7.0,height=7.0)
plot(MM$time[-1],diff(-MM$objfctList),log="y",xlim=c(0,0.16),ylim=c(1e-4,1),
     type="l",lwd=3,col=cols[1],xlab="Time",ylab="RSS difference",
     cex=1.5,cex.lab=1.4,cex.axis=1.3)
points(PCD$time[-1],diff(-PCD.RSS),type="l",lwd=2,col=cols[2])
points(PGDs$time[-1],diff(-PGDs.RSS),type="l",lwd=2,col=cols[3])
points(PGDl$time[-1],diff(-PGDl.RSS),type="l",lwd=2,col=cols[4])
points(MuEM$time[-1],diff(-MuEM.RSS),type="l",lwd=2,col=cols[6])
points(FCEM$time[-1],diff(-FCEM.RSS),type="l",lwd=2,col=cols[7])
abline(v=CNP$time,col=cols[8],lty=2,pch=4,lwd=2,cex=1.5)
legend('topright',algs[-5],col=cols[-5],
       lty=c(rep(1,7),2)[-5],lwd=2,bty='n')
lines(c(-0.1,0.12),c(0.1,0.1),lty=2,col="black",lwd=1)
abline(h=0.01,lty=2,col="black",lwd=1)
abline(h=0.001,lty=2,col="black",lwd=1)
#dev.off()
