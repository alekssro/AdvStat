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

setwd("C:/Users/Ky/Downloads/")

alphaEst <- 0.1
lb1<-3
lb2<-1
y<-ObsSeq
## Number of iterations
nIter <- 100
## Run EM algorithm
for (iter in 1:nIter){
  xPrb <- alphaEst*dpois(y,lambda = lb1)/
    (alphaEst*dpois(y,lambda = lb1)+
       (1-alphaEst)*dpois(y,lambda = lb2))
  alphaEst <- sum(xPrb)/length(y)
  lb1 <- sum(xPrb*y)/sum(xPrb)
  lb2 <-sum((1-xPrb)*y)/(length(y)-sum(xPrb))
  logLk <- sum(log(alphaEst*dpois(y,mean=muEst,sd=1)+(1-alphaEst)*dnorm(y,mean=0,sd=1)))
  cat("Iteration:",iter,
      "alpha:",alphaEst,"Lamba1:",lb1,"Lamba2", lb2,"\n")
  
}


###Part 3: HMM

## Specify initial parameter values
source("HMMexpectation.R")
IntPrb <- c(1/2,1/2)
TrnsPrb <- matrix(c(0.8,0.2,
                    0.1,0.9),byrow=TRUE,nrow=2,ncol=2)

lb<-c(1,3)

poisEM<-matrix(0,nrow = 2, ncol = 8)

for(i in 0:7){
  for(j in 1:2){
    poisEM[j,(i+1)]<-dpois(i, lambda =lb[j] )
  }
}
## Number of iterations
nIter <- 100
for (iter in 1:nIter){
  cat("--------------------------------------","\n")
  HMMexpct <- HMMexpectationsFct(IntPrb,TrnsPrb,poisEM,y)
  cat("\n")
  cat("Iteration:",iter,"log-likelihood:",log(HMMexpct$Lk),"\n")
  TrnsPrb <- HMMexpct$TransCnt/rowSums(HMMexpct$TransCnt)
  cat("Updated transition probability matrix:","\n")
  print(TrnsPrb)
  lb[1] <- sum(HMMexpct$PostProb[,1]*y/sum(HMMexpct$PostProb[,1]))
  lb[2] <- sum(HMMexpct$PostProb[,2]*y/sum(HMMexpct$PostProb[,2]))
  cat("Updated emission probabilities:","\n")
  cat(lb[1],lb[2],"\n")

  for(i in 0:7){
    for(j in 1:2){
      poisEM[j,(i+1)]<-dpois(i, lambda =lb[j] )
    }
  }

}

lines(s,alphaEst*dnorm(s,mean=muEst,sd=1),col="red",lwd=4,lty=2)
lines(s,(1-alphaEst)*dnorm(s,mean=0,sd=1),col="orange",lwd=4,lty=2)
lines(s,alphaEst*dnorm(s,mean=muEst,sd=1)+(1-alphaEst)*dorm(s,mean=0,sd=1),
      col="blue",lwd=2,lty=2)
## Add legend to plot
legend("topright",
       c("True mixture distribution","Estimated mixture distribution"),
       bty="n",col=c("black","black"),lty=c(1,2),lwd=3)


source("myLibrary.R")



emat <- matrix(c(sapply(unique(ObsSeq),dpois,lambda=lb[1]),sapply(unique(ObsSeq),dpois,lambda=lb[2])),byrow = T,nrow=2)

vit1 <- calcViterbi(c(0.5,0.5),TrnsPrb,emat,ObsSeq3)


mtransprobs <- function(t,j,i,backtrack,transition,emission,obsseq){
  
  out <- (backtrack[t,j]*transition[i,j]*emission[j,y[t]])/backtrack[t-1,i]
  return(out)
  
}

mtprobsA <-  sapply(c(2:225), mtransprobs, j = 1,i=1,backtrack = HMMexpct$backward,transition = TrnsPrb,emission = poisEM,obsseq = y)
mtprobsB <-  sapply(c(2:225), mtransprobs, j = 2,i=2,backtrack = HMMexpct$backward,transition = TrnsPrb,emission = poisEM,obsseq = y)

simHMM <- function(len,transmat){
  
  HSeq <- c(sample(c(1,2),size = 1,prob = c(0.5,0.5)))
  
  for(i in (2:len)){
    HSeq <- append(HSeq,sample(c(1,2),size = 1,prob = transmat[HSeq[i-1],]))
  }
  return(HSeq)
}

myDistlop <- function(len,len2){
  disti <- c()
  for(i in (1:len)){
    simA <- simHMM(len2,TrnsPrb)
    fracHigh <- table(simA)[2]/(table(simA)[1]+table(simA)[2])
    print(table(simA))
    disti <-  append(disti,fracHigh)
  }
  return(disti)
}





