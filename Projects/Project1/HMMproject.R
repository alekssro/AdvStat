## NOTE:
## The HMM expectation function HMMexpectationsFct
## (can be downloaded from the homepage: HMMexpectations.R)
## need to be active for the EM algorithm below to work
# source("HMMexpectations.R")
##----------------------------------
## Read sequences from file
## Analyse sequences
##----------------------------------
## Read sequences from file
browser()
dat <- dget("PhyloFoot.dat")
S1 <- dat$S1
S2 <- dat$S2

S1.num <- as.numeric(strsplit(chartr("AGCT","1234",S1),"")[[1]])
S2.num <- as.numeric(strsplit(chartr("AGCT","1234",S2),"")[[1]])

obs.seq <- 4*(S1.num-1)+S2.num
## Estimate transition and emission matrix using EM
diff.indx <- rep(0,16)
for (i in 1:4){  ## (A,G,C,T)
  for (j in 1:4){  ## (A,G,C,T)
    if (i!=j) diff.indx[4*(i-1)+j] <- 1
  }
}
diff.seq <- diff.indx[obs.seq]
cat("Number of different sites:",sum(diff.seq),"\n")
##-------------------------------------------------
## Emission probability function
prb.di.1 <- 0.1 ## Functionally important
prb.di.2 <- 0.3 ## Neutrally evolving
EmsPrbFct <- function(prb.di.1,prb.di.2){
  EmsPrb <- matrix(0,nrow=2,ncol=16)
  for (i in 1:4){ ## (A,G,C,T)
    for (j in 1:4){ ## (A,G,C,T)
      if (i!=j) EmsPrb[1,4*(i-1)+j] <- prb.di.1/12
      if (i!=j) EmsPrb[2,4*(i-1)+j] <- prb.di.2/12
      if (i==j) EmsPrb[1,4*(i-1)+j] <- (1-prb.di.1)/4
      if (i==j) EmsPrb[2,4*(i-1)+j] <- (1-prb.di.2)/4
    }
  }
  return(EmsPrb)
}
EmsPrb <- EmsPrbFct(prb.di.1,prb.di.2)
##-------------------------------------------------------
## Specify initial parameter values
IntPrb <- c(1/2,1/2)
TrnsPrb <- matrix(c(0.8,0.2,
                    0.1,0.9),byrow=TRUE,nrow=2,ncol=2)
prb.di.1 <- 0.2
prb.di.2 <- 0.4
EmsPrb <- EmsPrbFct(prb.di.1,prb.di.2)
## Number of iterations
nIter <- 20
for (iter in 1:nIter){
  cat("--------------------------------------","\n")
  HMMexpct <- HMMexpectationsFct(IntPrb,TrnsPrb,EmsPrb,obs.seq)
  cat("\n")
  cat("Iteration:",iter,"log-likelihood:",log(HMMexpct$Lk),"\n")
  TrnsPrb <- HMMexpct$TransCnt/rowSums(HMMexpct$TransCnt)
  cat("Updated transition probability matrix:","\n")
  print(TrnsPrb)
  N.di.1 <- sum(HMMexpct$PostProb[,1]*diff.seq)
  N.sa.1 <- sum(HMMexpct$PostProb[,1]*(1-diff.seq))
  N.di.2 <- sum(HMMexpct$PostProb[,2]*diff.seq)
  N.sa.2 <- sum(HMMexpct$PostProb[,2]*(1-diff.seq))
  cat("Updated emission probabilities:","\n")
  cat(prb.di.1,prb.di.2,"\n")
  prb.di.1 <- N.di.1/(N.di.1+N.sa.1)
  prb.di.2 <- N.di.2/(N.di.2+N.sa.2)
  EmsPrb <- EmsPrbFct(prb.di.1,prb.di.2)
}

