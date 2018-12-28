################################################################
## EXERCISES FROM HANDOUT 2 (Exercises in Week 36)
################################################################
#
#
####################### Exc 1 ##################################
#   1.Write a computer program in R where you reproduce the numbers
#           for N5;2 in Table 1 in Fu and Koutras (1994).
#

# to get possible cases we do n/k rounded down so N5,2 = 2, so we have to cover cases 0,1,2
# N5,2(x = 0) = ?
# N5,2(x = 1) = ?
# N5,2(x = 2) = ?
# N5,2(x = 0) = 'initial state prob' * transMatrix^5 * 'state we are looking for'
initalState <- c(1,0,0,0,0,0)
gettingSOnBegining <- c(1,1,0,0,0,0)
gettingSInMid <- c(0,0,1,1,0,0)
gettingSOnEnd <- c(0,0,0,0,1,1)

P <- function(t){
  return( 1/(t + 1))    # see "Distribution of N_n,k
}

# generating trans matrix in time
# following the pattern from paper, moving p and q by one to the right and always calc 1/t+1
transMatrix <- function(t){
  trans_matrix <- matrix(0, ncol = 6, nrow = 6)     # initialize with 0
  trans_matrix[1, 1:2] = trans_matrix[2, c(1,3)] = trans_matrix[3, c(3,4)] =
      trans_matrix[4, c(3,5)] = trans_matrix[5, c(5,6)] = c(1-P(t), P(t)) 
  trans_matrix[6,6] = 1
  return(trans_matrix)
}
# starting with time 1
trans_Matrix <- transMatrix(1)
# Generate trans mat in time 5
for(t in 2:5) {
  trans_Matrix = trans_Matrix %*% transMatrix(t)
}

# Package expm to take the power of matrix...

#using formula up there we can compute prob of states in time 5
# N5,2(x = 0) = 
initalState %*% trans_Matrix %*% gettingSOnBegining

# N5,2(x = 1) = 
initalState %*% trans_Matrix %*% gettingSInMid

# N5,2(x = 2) = 
initalState %*% trans_Matrix %*% gettingSOnEnd


####################### Exc 2 ##################################
#   2. Write a computer program in R where you reproduce the numbers 
#           for N15;2 in Table 2 in Fu and Koutras (1994).

# Same as first, but n/k rounded up is 7 so...
# N15,2(x = 0) = ?
# N15,2(x = 1) = ?
# N15,2(x = 2) = ?
# N15,2(x = 3) = ?
# N15,2(x = 4) = ?
# N15,2(x = 5) = ?
# N15,2(x = 6) = ?
# N15,2(x = 7) = ?

# generating pattern for trans matrix
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

# making vectors of needed states
initalState <- rep(0,16)
initalState[1] <- 1

## all possible ending states(non overlapping SS)
everySecond = T
cnt = 0
for(j in 1:16){
    if(everySecond){
        
        finalState <- rep(0,16)
        finalState[j] <- 1
        finalState[j+1] <- 1
        everySecond = F
            
        # starting with time 1
        trans_Matrix <- transMatrix(1)
        
        # Generate trans mat in time 5
        for(t in 2:15) {
            trans_Matrix = trans_Matrix %*% transMatrix(t)
        }
        
        #using formula up there we can compute prob of states in time 5
        # N5,2(x = 0) = 
        prob <- initalState %*% trans_Matrix %*% finalState
        cat("N15,2(X = ",cnt,") = ",prob, "\n", sep = "")
        cnt = cnt +1
    }
    else{
        everySecond = T
    }
}



####################### Exc 3 ##################################
#   3. Write a computer program in R that simulates a homogeneous finite state Markov chain.
#       Input should be an initial distribution, a transition matrix, and the length of the desired chain.
#       Output should be a simulation from the Markov chain with the desired length.

generateMC <- function(initDist, trans_Matrix, n){
  chain <- c()
  chain <- c(chain,sample(1:length(initDist),1,replace=T, prob=initDist))
  for(i in 1:n){
    chain <- c(chain,sample(c(1:length(initDist)), 1, replace = T, prob = trans_Matrix[chain[i],]))
  }
  return(chain)
}

# Example:
generateMC(initDist = c(1,0,0),
           trans_Matrix = matrix( c(0.1,0.5,0.4,
                                    0.3,0.2,0.5,
                                    0.2,0.4,0.4), 
                                nrow=3, byrow =T),
           n = 10)


####################### Exc 4 ##################################
#   4. We are interested in the word 'CTT' in a number of DNA sequences. 
#       The null model is that letters A,C,G,T are independent and occur 
#       with frequencies (0.3,0.2,0.2,0.3). Write a computer program in R for 
#       calculating the probability for 0, 1, 2,  >=3 number of occurences of the word 'CTT'
#           in a DNA sequence of arbitrary length.
#       Hint: Use similar ideas as in Fu and Koutras (1994) to define an appropriate Markov chain.
#       Verify your calculations by a simulation study. 
#           (E.g. use the computer program from the previous exercise.)
                                                 
# Expected mean waiting time to see the word first time
# Theoretically is 1/PcPtPt 
Pc = 0.2
Pa = 0.3
Pg = 0.2
Pt = 0.3

1/(Pc*Pt*Pt) # 55.55

# we have to get this number but generating the very long chain thousand of 
# times counting the occurance of CTT and calculating mean the occurances. 
# We tried with 10000 chain length and 200 iterations
# Also chain generating function was used from exc 3
library(stringr)
time_to_obs <- c()
n <- 5000
for(i in 1:300){
    
    # generate sequence
    chain <- generateMC(c(1/4,1/4,1/4,1/4),matrix(rep(c(0.3,0.2,0.2,0.3),4), nrow=4, byrow=T), n)
    # count number off occurances of word CTT(244)
    count <- str_count(paste0(collapse = '',chain), pattern = '244')
    #Get the indexes in order to cheat and cut off the seq after the last occurance
    index <- gregexpr('244', paste0(collapse = '',chain))[[1]][-1] + 2
    # this will yield faster time to get to the aimed mean, 
    # so we dont have to iterate 1000 time with seq of len 10000... but it is cheating
    
    # ime_to_obs[i] <- (n-index)/count
    
    # This is proper way to do it
    time_to_obs[i] <- n/count
  
}

hist(time_to_obs)
abline(v=mean(time_to_obs), col="black")
abline(v=1/(Pc*Pt*Pt), col="red")


