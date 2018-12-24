
# Functions
##############################################################
make_matrix <- function(N,t){
    
    P <- P <- matrix(0, nrow = N+1, ncol = N+1)
    
    pt <- 1/(t+1)
    # print(pt)
    qt <- 1-pt
    
    for(j in seq(1, dim(P)[1], by=k)){
        
        print(j)
        P[j,j] <- qt
        P[j,j+1] <- pt
        P[j+1,j] <- qt
        if (j!=1) {
            P[j-1,j] = pt
        }
    }
    
    j = dim(P)[1]
    P[j,j] <- 1
    P[j,j-1] <- 0
    
    for (t in 1:N) {
        
        
        
    }
    
    return(P)
    
}

# End of Functions
#######################################################################

# Write a computer program in R where you reproduce the numbers for N5,2 
# in Table 1 in Fu and Koutras (1994).

# We need to find the probabilities and the mean 

t <- c(1:5)
pt <- 1/(t+1)
qt <- 1-pt

InitProb <- c(1, 0, 0, 0, 0, 0)
U_prime <- c(1, 1, 0, 0, 0, 0)


P <- make_matrix(5,2)


# P <- matrix(0, nrow = 6, ncol = 6)
# P[1,1:2] <- c(qt[1], pt[1])
# P[2,1:3] <- c(qt[2], 0, pt[2])
# P[3,3:4] <- c(qt[3], pt[3])
# P[4,3:5] <- c(qt[4], 0, pt[4])
# P[5,5:6] <- c(qt[5], pt[5])
# P[6,6] <- 1

apply(P, 1, sum)  # check if is a real transition matrix 

result <- InitProb %*% P
probs <- c()
for (i in 1:4) {
  probs <- append(probs, result[4])
  result <- result %*% P
}
result %*% U_prime
plot(probs)
