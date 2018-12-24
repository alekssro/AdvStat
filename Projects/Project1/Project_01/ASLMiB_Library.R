

# library(stringr)

# EM algorithm, estimates mu and alpha given the observed values
# Stops when logLk difference is smaller then stopThres
myEM <- function(y, alphaEst, muEst, maxIter, stopThres = 0.0001){
    
    old_logLk <- 0
    ## Run EM algorithm
    for (iter in 1:maxIter){
        xPrb <- alphaEst*dnorm(y,mean=muEst,sd=1)/
            (alphaEst*dnorm(y,mean=muEst,sd=1)+
                 (1-alphaEst)*dnorm(y,mean=0,sd=1))
        alphaEst <- sum(xPrb)/length(y)
        muEst <- sum(xPrb*y)/sum(xPrb)
        logLk <- sum(log(alphaEst*dnorm(y,mean=muEst,sd=1)+
                             (1-alphaEst)*dnorm(y,mean=0,sd=1)))
        
        if ((abs(old_logLk - logLk) < stopThres) & iter > 5) {
            cat("Converged after", iter, "iteration:\n",
                "alpha:",alphaEst,"mu:",muEst,
                "logLk",logLk,"\n")
            break
        }
        
        old_logLk <- logLk
        
        cat("Iteration:",iter,
            "alpha:",alphaEst,"mu:",muEst,
            "logLk",logLk,"\n")
    }
    plot(xPrb, y)
}

# Simulates HMM and return hidden states and observed states
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


# Generates markov chain with specified arguments.
# Example: 
  #generateMC(initDist = c(1,0,0),
  #           trans_Matrix = matrix( c(0.1,0.5,0.4,
  #                                    0.3,0.2,0.5,
  #                                    0.2,0.4,0.4), 
  #                                  nrow=3, byrow =T),
  #           n = 10)
generateMC <- function(initDist, trans_Matrix, n){
  chain <- c()
  chain <- c(chain,sample(1:length(initDist),1,replace=T, prob=initDist))
  for(i in 1:n){
    chain <- c(chain,sample(c(1:length(initDist)), 1, replace = T, prob = trans_Matrix[chain[i],]))
  }
  return(chain)
}

# Uses viterbi algorithm to estimate probability of observing seq and also
# returns mosst likely sequences given the arguments
calcViterbi<- function(initProb, transProb, emissionProb, obsSeq, DEBUG = F){
  
  # initial params
  lenSeq <- length(obsSeq)
  numHiddenStates <- nrow(transProb)
  
  MaxArrow <- matrix(0,nrow=lenSeq,ncol=numHiddenStates)
  
  # empty sequence for backtracking
  seqOfStates <- rep(0, lenSeq)
  # initialize empty matrix of size seq len x hidden states
  score_matrix <- matrix(0, nrow = lenSeq, ncol = numHiddenStates)
  # first row of matrix
  score_matrix[1,] <- initProb*emissionProb[,obsSeq[1]]
  
  
  for(charInSeq in 2:lenSeq){
    for(hidState in 1:numHiddenStates){
      # look at previous row
      previousRow <- score_matrix[charInSeq-1, ]
      # look at emission prob given state to obs char
      emission4givenState <- emissionProb[hidState, obsSeq[charInSeq]]
      # trans probs of coming to specified hidState
      transitionProb <- transProb[, hidState ]
      
      product <- previousRow  * transitionProb * rep(emission4givenState, numHiddenStates)
      # save the max value to the position in matrix
      score_matrix[charInSeq, hidState] <- max(product)
      # Remember where did max value came from and store it in matrix
      MaxArrow[charInSeq,hidState] <- which.max(product)
      
      if(DEBUG){
        cat("##################\n")
        cat("TransProb j (1 to hiddenstates): ", transitionProb[,hidState])
        cat("\nRow before: \n", score_matrix[charInSeq-1,])
        cat("\nEMisson prob j, char at k: repeated nHS times:  ", rep(emissionProb[charInSeq,obsSeq[charInSeq]],numHiddenStates))
        cat("\n",product)
        cat("\nMax Value", max(product))
        print(score_matrix)
        cat("\n##################\n")
      }
    }
    
  }
  # Empty seq
  BackTrack <- rep(0,lenSeq)
  # We start from last row and take the max
  BackTrack[lenSeq] <- which.max(score_matrix[lenSeq,])
  # For every seq in obs seq we look at index where the max came from(char at i, row from MaxArrow)
  for (i in lenSeq:2) {
    BackTrack[i-1] <- MaxArrow[i,BackTrack[i]]
  }
  out <- list()
  out$backtrack <- BackTrack
  out$score_matrix <- score_matrix
  out$jointProb <- max( score_matrix[charInSeq, ])
  
  return(out)
}

# Generates markov chain with provided arguments and 
# looks for specific pattern and calculates time to observe
generateMC_and_estMeanTimeToObserve <- function(iterNum, 
                                                seqLen, 
                                                transProb = matrix(rep(c(0.3,0.2,0.2,0.3),4), nrow=4, byrow=T),
                                                initProb = c(1/4,1/4,1/4,1/4),
                                                pattern = '244'){
  
  time_to_obs <- c()
  n <- seqLen
  for(i in 1:iterNum)
  {
    # generate sequence
    chain <- generateMC(initProb,transProb, n)
    # count number off occurances of word CTT(244)
    count <- str_count(paste0(collapse = '',chain), pattern = pattern)
    #Get the indexes in order to cheat and cut off the seq after the last occurance
    index <- gregexpr('244', paste0(collapse = '',chain))[[1]][-1] + 2
    # this will yield faster time to get to the aimed mean, 
    # so we dont have to iterate 1000 time with seq of len 10000... but it is cheating
    
    # ime_to_obs[i] <- (n-index)/count
    
    # This is proper way to do it
    time_to_obs[i] <- n/count
    
  }
  out = c()
  out$time_of_observations = time_to_obs
  out$meanTime = mean(time_to_obs)
  return(out)
}

## Generate transition matrix in time t
# Change cols and row for more flexible usage
transMatrix <- function(t){
  trans_matrix <- matrix(0, ncol = 16, nrow = 16)
  for(i in 1:15){
    if(i %% 2 == 0){
      trans_matrix[i, c(i-1, i+1)] = c(1-P(t), P(t)) 
    }else
    {
      trans_matrix[i, c(i,i+1)] = c(1-P(t), P(t)) 
    }
  }
  trans_matrix[16,16] = 1
  return(trans_matrix)
}

## Backward forward algorithm 
forwardBackward <- function(InitProb,TransProb,EmisProb,ObsSeq, verbose = F){
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
    if(verbose)
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
    if(verbose)
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

