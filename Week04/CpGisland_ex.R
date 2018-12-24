source("../AsgerLibrary.R")

## Read sequence from file
CpGdat <- readLines("CpGisland.dat")
## Split and translate raw sequence from AGCT to 1234
## and make the vector numeric
ObsSeq <-
    as.numeric(strsplit(chartr("AGCT","1234",CpGdat),"")[[1]])
pi <- c(0.5, 0.5)
p <- myEM(ObsSeq, 1000, 0.5, 1)

##-------------------------------------
## Estimate transition matrix using EM
##-------------------------------------
source("HMMexpectations.R")
## Initial parameter values
TrnsPrb <- matrix(c(.5,.25,0.25,
                    .25,0.5,0.25,
                    0.25,0.25,0.5),
                  nrow=3,ncol=3,byrow=TRUE)
## Number of iterations
nIter <- 25
for (iter in 1:nIter){
    cat("--------------------------------------","\n")
    HMMexpct <- HMMexpectationsFct(IntPrb,TrnsPrb,EmsPrb,HMMsim$ObsSeq)
    cat("\n")
    cat("Iteration:",iter,"log-likelihood:",log(HMMexpct$Lk),"\n")
    TrnsPrb <- HMMexpct$TransCnt/rowSums(HMMexpct$TransCnt)
    cat("Updated transition probability matrix:","\n")
    print(TrnsPrb)
}