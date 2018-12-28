
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

# DISHONEST CASINO

initPr <- c(0.5, 0.5)
transPr <- matrix(c(3/4, 1/4, 1/2, 1/2), nrow = 2, byrow = T)
EmisPr <- matrix(c(1/2, 1/2, 2/3, 1/3), nrow = 2, byrow = T)

simRes <- HMMsimFct(InitProb = initPr, TransProb = transPr, EmisProb = EmisPr, 100)
ifelse(simRes$ObsSeq == 1, 'H', 'T')
