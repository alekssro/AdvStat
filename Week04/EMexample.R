source("HMMsim.R")
##------------------------------
## Simulate sequence from HMM
##------------------------------
## Three hidden states
## Initial probabilities
IntPrb <- c(1/3,1/3,1/3)
## Transition probabilities
TrnsPrb <- matrix(c(.6,.2,.2,
                    .3,.4,.3,
                    .1,.1,.8),
                  byrow=TRUE,nrow=3,ncol=3)
## Four possible emissions
EmsPrb <- matrix(c(1/2,1/4,1/8,1/8,
                   1/8,3/4,1/16,1/16,
                   1/10,1/10,2/5,2/5),
                 byrow=TRUE,nrow=3,ncol=4)
## Simulated sequence of length 100
Ln <- 300
HMMsim <- HMMsimFct(IntPrb,TrnsPrb,EmsPrb,Ln)
cat("Simulated sequence:","\n")
print(HMMsim)
##-------------------------------------
## Estimate transition matrix using EM
##-------------------------------------
source("HMMexpectations.R")
## Initial parameter values
TrnsPrb <- matrix(c(.5,.25,0.25,
                    .25,0.5,0.25,
                    0.25,0.25,0.5),
                  nrow=3,ncol=3,byrow=TRUE)
## Number of iterations
nIter <- 25
for (iter in 1:nIter){
  cat("--------------------------------------","\n")
  HMMexpct <- HMMexpectationsFct(IntPrb,TrnsPrb,EmsPrb,HMMsim$ObsSeq)
  cat("\n")
  cat("Iteration:",iter,"log-likelihood:",log(HMMexpct$Lk),"\n")
  TrnsPrb <- HMMexpct$TransCnt/rowSums(HMMexpct$TransCnt)
  cat("Updated transition probability matrix:","\n")
  print(TrnsPrb)
}


