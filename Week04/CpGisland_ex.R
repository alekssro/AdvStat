# source("HMMexpectations.R")
# source("HMMsim.R")
source("ASLMiB_package.R")

EmisPr <- matrix(c(1/8, 3/8, 3/8, 1/8,
                   1/4, 1/4, 1/4, 1/4), byrow = T, nrow = 2)

# a)
## Read sequence from file
CpGdat <- readLines("CpGisland.dat")
## Split and translate raw sequence from AGCT to 1234
## and make the vector numeric
ObsSeq <-
    as.numeric(strsplit(chartr("AGCT","1234",CpGdat),"")[[1]])
pi <- c(0.5, 0.5)

# b)
##-------------------------------------
## Estimate transition matrix using EM
##-------------------------------------
# source("HMMexpectations.R")
## Initial parameter values
TrnsPrb <- matrix(c(.5, .5,
                    .5, .5),
                  nrow=2,ncol=2,byrow=TRUE)
## Number of iterations
nIter <- 25
for (iter in 1:nIter){
    cat("--------------------------------------","\n")
    HMMexpct <- HMMexpectationsFct(pi, TrnsPrb, EmisPr, ObsSeq)
    cat("\n")
    cat("Iteration:",iter,"log-likelihood:",log(HMMexpct$Lk),"\n")
    TrnsPrb <- HMMexpct$TransCnt/rowSums(HMMexpct$TransCnt)
    cat("Updated transition probability matrix:","\n")
    print(TrnsPrb)      # Estimated transition probability matrix
}

# c)
# Study CpG islands
par(mfrow=c(2,1))

# Viterbi decoding
viterbiRes <- calcViterbi(pi, TrnsPrb, EmisPr, ObsSeq)
plot(viterbiRes$backtrack, type="l",  main = "Viterbi Decoding",
     xlab="Sequence index",ylab="Hidden state")

# Posterior decoding
posterior <- apply(HMMexpct$PostProb, 1, which.max)
plot(posterior, type="l", main = "Posterior Decoding",
     xlab="Sequence index", ylab="Hidden state")