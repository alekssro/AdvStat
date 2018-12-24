
calcViterbi<- function(initProb, transProb, emissionProb, obsSeq){
  
  # initial params
  lenSeq <- length(obsSeq)
  numHiddenStates <- nrow(transProb)
  # initialize empty matrix of size seq len x hidden states
  score_matrix <- matrix(0, nrow = lenSeq, ncol = numHiddenStates)
  # initialize variables for backtracking
  BackTrack <- rep(0,lenSeq)
  MaxArrow <- matrix(0,nrow=lenSeq,ncol=numHiddenStates)
  

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
      
    }
    
  }
  
  
  BackTrack[lenSeq] <- which.max(score_matrix[lenSeq,])
  for (i in lenSeq:2) {
    BackTrack[i-1] <- MaxArrow[i,BackTrack[i]]
  }
  out <- list()
  out$backtrack <- BackTrack
  out$score_matrix <- score_matrix
  out$jointProb <- max( score_matrix[charInSeq, ])
  
  return(out)
}