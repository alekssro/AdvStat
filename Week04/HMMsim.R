# Name: HMMsimFct
# Author: Asger Hobolth
# Purpose:
# Simulation from a HMM of specified length
#--------------------------------------------------------
# Input:
# InitProb: Initial probabilities for hidden states
# TransProb: Transition probabilities;
#            Probabilities between hidden states
#            nHS times nHS matrix
# EmisProb: Emission probabilities
# Len: Length of hidden and emitted sequences
#
# Output:
# Sim: List of hidden and emitted sequence
# $HidSeq and $ObsSeq
#---------------------------------------------------------
HMMsimFct <- function(InitProb,TransProb,EmisProb,len){
  HidSeq <- rep(0,len)
  ObsSeq <- rep(0,len)
  nHS <- nrow(TransProb)
  nE <- ncol(EmisProb)
  HidSeq[1] <- sample(1:nHS,size=1,replace=TRUE,InitProb)
  ObsSeq[1] <- sample(1:nE,size=1,replace=TRUE,EmisProb[HidSeq[1],])
  for (i in 2:len){
    HidSeq[i] <- sample(1:nHS,size=1,replace=TRUE,TransProb[HidSeq[i-1],])
    ObsSeq[i] <- sample(1:nE,size=1,replace=TRUE,EmisProb[HidSeq[i],])
  }
  Sim <- list()
  Sim$HidSeq <- HidSeq
  Sim$ObsSeq <- ObsSeq
  return(Sim)
}
##------------------------------
## Example of HMM simulation
##------------------------------
## Three hidden states
## Initial probabilities
IntPrb <- c(1/3,1/3,1/3)
## Transition probabilities
TrnsPrb <- matrix(c(0.5,0.2,0.3,
                    0.3,0.4,0.3,
                    0.1,0.1,0.8),
                  byrow=TRUE,nrow=3,ncol=3)
## Four possible emissions
EmsPrb <- matrix(c(1/2,1/4,1/8,1/8,
                   1/8,3/4,1/16,1/16,
                   1/10,1/10,2/5,2/5),
                 byrow=TRUE,nrow=3,ncol=4)
## Simulated sequence of length 100
Ln <- 100
HMMsim <- HMMsimFct(IntPrb,TrnsPrb,EmsPrb,Ln)
print(HMMsim)
plot(HMMsim$HidSeq,type="l",
     xlab="Sequence index",ylab="Hidden state")
