ObsSeq <- c(0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,
            0,1,0,0,0,0,2,2,0,0,0,0,1,0,0,
            1,1,0,0,1,1,1,0,0,1,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,
            1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,1,0,0,0,0,0,7,3,2,3,2,4,
            0,0,0,0,1,0,0,0,0,0,0,0,1,0,2,
            0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
            2,1,0,0,1,0,0,0,1,0,1,1,0,0,0,
            1,0,0,1,0,0,0,1,2,0,0,0,1,0,1,
            1,0,1,0,0,2,0,1,2,1,1,2,1,0,1,
            1,0,0,1,1,0,0,0,1,1,1,0,4,0,0,
            2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)


# The following three models seem relevant and natural:

# i) Single Poisson (get lambda)

obs_var <- sd(ObsSeq)
one_poiss <- rpois(length(ObsSeq), mean(ObsSeq))
sd(one_poiss); obs_var
# observed sd doesn't match to a sd of a single poisson distribution with lambda=mean(observations)

# ii) A mixture of two Poisson distributions (corresponding to low or high activity)

# States: 0 -> not moving ; 1 -> moving
# We have a different lambda for each state (2 different poisson distributions) - lambda_1 , lambda_2
# Probability of not moving is defined by alpha

# Initial parameter values:
alphaEst <- 0.5
lambda_1 <- 1
lambda_2 <- 1.5

# EM algorithm, estimates alpha, lambda_1 and lambda_2 given the observed values
# Stops when logLk difference is smaller than stopThres
myEM <- function(y, alphaEst, lambda1Est, lambda2Est, maxIter, stopThres = 0.0001){
    
    old_logLk <- 0
    ## Run EM algorithm
    for (iter in 1:maxIter){
        xPrb <- alphaEst*dpois(y,lambda = lambda1Est)/
            ((alphaEst*dpois(y,lambda = lambda1Est))+
                 (1-alphaEst)*dpois(y, lambda = lambda2Est))
        alphaEst <- sum(xPrb)/length(y)
        lambda1Est <- sum(xPrb*y)/sum(xPrb)
        lambda2Est <- sum((1-xPrb)*y)/sum(1-xPrb)
        # FOLLOW FROM HEREEE
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



# iii) A HMM with two hidden states (corresponding to low or high activity) and Poisson emissions.

# Fit and compare these three models,
# and discuss how they are different. Classify (i.e. decode) the observations into low or high activity
# in model (ii) and (iii). Are the classifiers different? Why? Which models do you prefer? Which
# classifications do you prefer?
