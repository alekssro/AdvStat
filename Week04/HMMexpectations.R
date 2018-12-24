# Name: HMMexpectationsFct.R
# Author: Asger Hobolth
# Purpose:
# Calculates expected transition count matrix
# Calculates likelihood
# Calculates posterior probability
#--------------------------------------------------------
# Input:
# InitProb: Initial probabilities
# TransProb: Transition probabilities;
#            Probabilities between hidden states
#            nHS times nHS matrix
# EmisProb: Emission probabilities
# ObsSeq: Observed sequence
#
# Output:
# $TransCnt: Expected transition count matrix
# $PostProb: Posterior probability matrix
# $Lk: Log likelihood
#---------------------------------------------------------
HMMexpectationsFct <- function(InitProb,TransProb,EmisProb,ObsSeq){
  len <- length(ObsSeq)
  nHS <- nrow(TransProb)
  #-------------------
  # Forward algorithm 
  #-------------------
  # Define ForwardLik matrix
  ForwardLik <- matrix(0,nrow=len,ncol=nHS)
  # Start condition
  ForwardLik[1,] <- InitProb*EmisProb[,ObsSeq[1]]
  # Determine ForwardLik by recursion
  for (k in 2:len){
    for (j in 1:nHS){
      ForwardLik[k,j] <- sum(TransProb[,j]*
                             rep(EmisProb[j,ObsSeq[k]],nHS)*
                             ForwardLik[k-1,])
    }
  }
  ForwardLikVal <- sum(ForwardLik[len,])
  cat("Likelihood from Forward algorithm:",ForwardLikVal,"\n")
  #--------------------
  # Backward algorithm
  #--------------------
  # Define BackwardLik
  BackwardLik <- matrix(0,nrow=len,ncol=nHS)
  # Start condition
  BackwardLik[len,] <- rep(1,nHS)
  # Determine logBackwardLik by recursion
  for (k in (len-1):1){
    for (j in 1:nHS){
      BackwardLik[k,j] <- sum(TransProb[j,1:nHS]*
                              EmisProb[,ObsSeq[k+1]]*
                              BackwardLik[k+1,])
    }
  }
  BackwardLikVal <- sum(InitProb*
                        EmisProb[,ObsSeq[k]]*
                        BackwardLik[1,])
  cat("Likelihood from Backward algorithm:",BackwardLikVal,"\n")
  ##-----------------------
  ## Posterior probability
  ##-----------------------
  PostProb <- exp(log(BackwardLik)+log(ForwardLik)-log(BackwardLikVal))
  ##----------------------------
  ## Expected transition counts
  ##----------------------------
  TransCnt <- matrix(0,nrow=nHS,ncol=nHS)
  for (k in 1:nHS){
    for (l in 1:nHS){
      Probkl <- ForwardLik[1:(len-1),k]*BackwardLik[2:len,l]/
        BackwardLikVal*TransProb[k,l]*EmisProb[l,ObsSeq[2:len]]
      TransCnt[k,l] <- sum(Probkl)
    }
  }
  output <- list()
  output$TransCnt <- TransCnt
  output$PostProb <- PostProb
  output$Lk <- BackwardLikVal
  return(output)
}

