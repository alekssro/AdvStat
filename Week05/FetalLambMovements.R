##-------------------------------------------
## Fetal Lamb Movements
##-------------------------------------------
## Data: Guttorp (1995) Exercise d9 page 124
##-------------------------------------------
ObsSeq <- c(0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,
            0,1,0,0,0,0,2,2,0,0,0,0,1,0,0,
            1,1,0,0,1,1,1,0,0,1,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,
            1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,1,0,0,0,0,0,7,3,2,3,2,4,
            0,0,0,0,1,0,0,0,0,0,0,0,1,0,2,
            0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
            2,1,0,0,1,0,0,0,1,0,1,1,0,0,0,
            1,0,0,1,0,0,0,1,2,0,0,0,1,0,1,
            1,0,1,0,0,2,0,1,2,1,1,2,1,0,1,
            1,0,0,1,1,0,0,0,1,1,1,0,4,0,0,
            2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
plot(ObsSeq,type="h")
table(ObsSeq)
##---------------------------------------------
## A single Poisson is not appropriate
## for example compare observed and 
## simulated variance
##---------------------------------------------
var(ObsSeq) ###0.6925
sim <- rpois(n=length(ObsSeq),mean(ObsSeq)) ###Simulating poisson distribution
var(sim) ###0.333
##---------------------------------------------
## Analyse the data using a mixture of two 
## Poission distributions. 
## Use EM for parameter estimation
##--------------------------------------------
## Starting values
rate1Est <- 0.1
rate2Est <- 3
alphaEst <- 0.5
## Number of iterations
nIter <- 100
## Run EM algorithm
for (iter in 1:nIter){
    xPrb <- alphaEst*dpois(ObsSeq,lam=rate1Est)/
            (alphaEst*dpois(ObsSeq,lam=rate1Est)+
            (1-alphaEst)*dpois(ObsSeq,lam=rate2Est))
    alphaEst <- sum(xPrb)/length(ObsSeq)
    rate1Est <- sum(xPrb*ObsSeq)/sum(xPrb)
    rate2Est <- sum((1-xPrb)*ObsSeq)/sum(1-xPrb)
    logLk <- sum(log(alphaEst*dpois(ObsSeq,lam=rate1Est)+
                 (1-alphaEst)*dpois(ObsSeq,lam=rate2Est)))
    cat("Iteration:",iter,
      "Rate1:",rate1Est,"Rate2:",rate2Est,"alpha:",alphaEst,
      "logLk:",logLk,"\n")
}
##-------------------------------------------------
## Illustrate final fit
##-------------------------------------------------
cnt <- tabulate(ObsSeq+1,8)
plot(0:7,cnt,pch=19)
points(0:7,alphaEst*length(ObsSeq)*dpois(0:7,lam=rate1Est),col="red",pch=19)
points(0:7,(1-alphaEst)*length(ObsSeq)*dpois(0:7,lam=rate2Est),col="blue",pch=19)
mixtEst <-  alphaEst*length(ObsSeq)*dpois(0:7,lam=rate1Est)+
                (1-alphaEst)*length(ObsSeq)*dpois(0:7,lam=rate2Est)
points(0:7,mixtEst,col="green",pch=19)
##-----------------------------------------------------
## Time series modelling
##------------------------------------------------------
## Plot along the sequence
##------------------------------------------------------
plot(ObsSeq,pch=19,cex=.5)
points(1:length(ObsSeq),1-xPrb,pch=19,col="red",cex=.5)
##-------------------------------------------------------
## Analyse the data using a two-state HMM
## with Poisson emissions
##-------------------------------------------------------
source("HMMPoisExpectations.R")
a <- 0.5 ; b <- 0.5 
TransProb <- matrix(c(a,1-a,1-b,b),byrow=TRUE,nrow=2)
lam1 <- 1 ; lam2 <- 3
Lambda <- c(lam1,lam2)
HMMRes <- HMMPoisExpectationsFct(InitProb=c(0.5,0.5),TransProb=TransProb,Lambda=Lambda,ObsSeq=ObsSeq)
lala <- ObsSeq * HMMRes$PostProb[,1]
## Run EM algorithm
nIter <- 100
for (i in 1:nIter){
    ## Estimation of rates
    lam1Est <- sum(ObsSeq*HMMRes$PostProb[,1])/sum(HMMRes$PostProb[,1])
    lam2Est <- sum(ObsSeq*HMMRes$PostProb[,2])/sum(HMMRes$PostProb[,2])
    Lambda <- c(lam1Est,lam2Est)
    ## Estimation of transitions
    aEst <- HMMRes$TransCnt[1,1] / sum(HMMRes$TransCnt[1, ])
    bEst <- HMMRes$TransCnt[2,2] / sum(HMMRes$TransCnt[2, ]) 
    TransProb <- matrix(c(aEst, 1-aEst, 1-bEst, bEst), byrow=TRUE, nrow=2)
    HMMRes <- HMMPoisExpectationsFct(InitProb=c(0.5,0.5),TransProb=TransProb,Lambda=Lambda,ObsSeq=ObsSeq)
    cat("Iteration", i, HMMRes$Lk, "\n")
}
HMMRes$TransCnt     # Expected counts
HMMRes$PostProb
points(1:length(ObsSeq), HMMRes$PostProb[,2], pch=19, col="blue", cex=.5)

##---------------------------------------------------------
## Viterbi algorithm
##---------------------------------------------------------
source("ViterbiPois.R")
HMMVit <- ViterbiPoisFct(InitProb=c(0.5,0.5),TransProb=TransProb,Lambda=Lambda,ObsSeq=ObsSeq)
HMMVit$BackTrack
points(1:length(ObsSeq),HMMVit$BackTrack-0.9,col="green",cex=.5,pch=19)
##-----------------------------------------------------------
## Expected under the three models
## M1: Poisson
## M2: Poisson mixture
## M3: Poisson HMM
##-----------------------------------------------------------
CntHMMest <- (1-bEst) / ( (1-aEst)+(1-bEst) )*dpois(0:7,lam=lam1Est)*length(ObsSeq)+
                (1-aEst) / ( (1-aEst)+(1-bEst) )*dpois(0:7,lam=lam2Est)*length(ObsSeq)
cnt
round(length(ObsSeq)*dpois(0:7,mean(ObsSeq)), digits=1)
round(CntHMMest, digits=1)
round(mixtEst, digits=1)     # Seems to be closer to the observations
