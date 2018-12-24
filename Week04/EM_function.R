myEMfunction <- function(y, alphaEst, muEst, maxIter, stopThres = 0.0001){
    
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