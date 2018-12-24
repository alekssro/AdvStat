source("viterbi.R")
source("MyViterbi.R")
##------------------------------
## Example of Viterbi algorithm
##------------------------------
## Three hidden states
## Initial probabilities
IntPrb <- c(1/3,1/3,1/3)
## Transition probabilities
TrnsPrb <- matrix(c(0.5,0.2,0.3,
                    0.3,0.4,0.3,
                    0.1,0.1,0.8),
                  byrow=TRUE,nrow=3,ncol=3)
## Four possible emissions
EmsPrb <- matrix(c(1/2,1/4,1/8,1/8,
                   1/8,3/4,1/16,1/16,
                   1/10,1/10,2/5,2/5),
                 byrow=TRUE,nrow=3,ncol=4)
## Observed sequence
ObsSq <- c(1,1,1,1,2,2,3,2,2,1,4,4,4,3,3,2,2,1,2,4,3,2,1,3)
## Apply Viterbi algorithm
ViterbiRes <- ViterbiFct(IntPrb,TrnsPrb,EmsPrb,ObsSq)
MyViterbiRes <- calcViterbi(IntPrb,TrnsPrb,EmsPrb,ObsSq)

ViterbiRes$MaxValue == MyViterbiRes$score_matrix
ViterbiRes$BackTrack == MyViterbiRes$backtrack

## Visualize result

Ln <- length(ObsSq)
par(mfrow=c(2,1))
plot(1:Ln,ObsSq,
     col="blue",type="l",
     main="Observed sequence",
     xlab="sequence index",ylab="Observed state", ylim = range(0:4))
# plot(1:Ln,ViterbiRes$BackTrack,
#      col="red",type="l",
#      main="Viterbi sequence",
#      xlab="sequence index",ylab="Decoded state",ylim = range(0:4))
plot(1:Ln,MyViterbiRes$backtrack,
     col="black",type="l",
     main="My Viterbi sequence",
     xlab="sequence index",ylab="Decoded state",ylim = range(0:4))
