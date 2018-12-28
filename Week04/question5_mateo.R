##################################################
## Project: Week 4 exc
## Script purpose: 
## Date: 20 Sep 2018
## Author: Mateo Sokac
##################################################

source("ASLMiB_package.R")

# Stupid way of estiamting prob matrix estimate...
# 1) do a forward backword for obs seq with initalProbMatrix (random guess or something smarter)
# 2) we get Trans count which can be used in formule 12.13 page 415
# 3) generate transMatrix from formula above
# 4) Use the newly generated transMatrix for backward forward,... iterate
estimateProbMatrix<- function(ObsSeq, initTransProb, initProb, emisProb, n_iter){
  
  simProbMatrices = data.frame(prob1=double(), prob2=double())
  for(i in 1:n_iter){
    if(i == 1){ # only first time we start with initTransProb
      res <- forwardBackward(initProb,
                             ObsSeq = ObsSeq, 
                             TransProb = initTransProb,
                             EmisProb = emisPrb)
    } else { # we use newly generated transProb
      res <- forwardBackward(initProb,
                             ObsSeq = ObsSeq, 
                             TransProb = currentProbMatrix,
                             EmisProb = emisPrb)
      
    }
    # formula 12.13 page 415
    N_jk_1 <- res$TransCnt[1,1]/(res$TransCnt[1,1] + res$TransCnt[1,2])
    N_jk_2 <- res$TransCnt[2,1]/(res$TransCnt[2,1] + res$TransCnt[2,2])
    
    # just to keep track of probs
    simProbMatrices <- rbind(simProbMatrices, data.frame(prob1 = N_jk_1, prob2 = N_jk_2))
    
    currentProbMatrix <- matrix(c(N_jk_1, 1 - N_jk_1, 
                                  N_jk_2, 1 - N_jk_2),
                                nrow = 2, ncol = 2, byrow = T
    )
    print(currentProbMatrix)
    
  }
  
  
  estTransMatrix <- matrix(c(mean(simProbMatrices$prob1), 1-mean(simProbMatrices$prob1),
                             mean(simProbMatrices$prob2), 1-mean(simProbMatrices$prob2)), 
                           byrow = T,ncol = 2,nrow = 2)
  
  out <- c()
  out$estTransMatrix <- estTransMatrix
  out$currentProbMatrix <- currentProbMatrix
  return(out)
}


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



#############################################################

# a)
CpGdat <- readLines("CpGisland.dat")
ObsSeq <-as.numeric(strsplit(chartr("AGCT","1234",CpGdat),"")[[1]])

# b) 

initProb = c(1/2,1/2)
emisPrb = matrix(c(1/8,3/8,3/8,1/8,
                   1/4,1/4,1/4,1/4), byrow = T, nrow = 2, ncol = 4)

randomGuess <- matrix(c(1/2,1/2,1/2,1/2),
                      byrow = T, ncol = 2, nrow = 2)

estimatedProbMatrix <- estimateProbMatrix(ObsSeq = ObsSeq, 
                                          initTransProb = randomGuess,
                                          initProb = initProb, 
                                          emisProb = emisPrb,
                                          n_iter = 2000)

estimatedProbMatrix$currentProbMatrix
estimatedProbMatrix$estTransMatrix

## C)
viterbi <- calcViterbi(initProb = initProb, 
                       transProb = estimatedProbMatrix$currentProbMatrix, 
                       emissionProb = emisPrb, 
                       obsSeq = ObsSeq)

viterbi$backtrack

