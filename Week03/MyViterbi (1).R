


calcViterbi<- function(initProb, transProb, emissionProb, obsSeq, DEBUG = F){
  
  # initial params
  lenSeq <- length(obsSeq)
  numHiddenStates <- nrow(transProb)
  
  MaxArrow <- matrix(0,nrow=lenSeq,ncol=numHiddenStates)
  
  # empty sequence for backtracking
  seqOfStates <- rep(0, lenSeq)
  # initialize empty matrix of size seq len x hidden states
  score_matrix <- matrix(0, nrow = lenSeq, ncol = numHiddenStates)
  backtrack <- c(0:lenSeq)
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
    
    # remember the index where the max values came from for the computed character
    backtrack[charInSeq-1] <- which.max( score_matrix[charInSeq-1, ] )
  }
  
  
  BackTrack <- rep(0,lenSeq)
  BackTrack[lenSeq] <- which.max(score_matrix[lenSeq,])
  for (i in lenSeq:2) {
    BackTrack[i-1] <- MaxArrow[i,BackTrack[i]]
  }
  out <- list()
  out$backtrack <- BackTrack
  out$myBacktrack <- backtrack
  out$score_matrix <- score_matrix
  out$jointProb <- max( score_matrix[charInSeq, ])
  
  return(out)
}