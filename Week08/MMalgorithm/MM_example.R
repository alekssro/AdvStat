source("NQPAlgorithms.R")
library(tidyverse)

# Example data
# X:=W, Y:=v, tht:=H
W = matrix(c(10,1,5,2), nrow=2, ncol=2) # weights
v = c(1,8) # data matrix
h_init = c(2,2) # signature (what we want to find)
A = 2*t(W)%*%W
b = as.vector(-2*t(W)%*%v)
true_val = c(-0.4,5)

#### Unconstrained
thtUncHat <- solve(t(W)%*%W)%*%t(W)%*%v # solution, h free
plot(1:2,h_init,pch=19,ylim=range(thtUncHat),cex=0.8, xlab="", ylab="", main="Least Square Unconstrained Solution")
abline(h=0,col="red")
points(1:2,thtUncHat,pch=19,col="blue",cex=0.8)
for (i in 1:2) {
  lines(c(i,i),c(h_init[i],thtUncHat[i]),col="blue")
}

plot(c(h_init[1], thtUncHat[1]), c(h_init[2], thtUncHat[2]), pch=19,col="blue",cex=0.8, xlab="", ylab="", main="Least Square Unconstrained Solution")
abline(v=0,col="red")
lines(c(h_init[1], thtUncHat[1]), c(h_init[2], thtUncHat[2]))

# Run MM algo
n.rep <- 100 ; tol <- 0.0001 ; mxIter <- 1100
MM <- MajorizeMinimizeNQP(A=A,b=b,init=h_init,maxIter=mxIter,tol=tol)
xList <- MM$xList
K <- length(xList[1,])
x <- c()
y <- c()
xList[1,][17]
xList[2,][17]

for(i in 1:K){
  png(paste0('img/NMF/test',i,'.png'), width=480, height=480)
  x <- c(x, xList[1,][i])
  y <- c(y, xList[2,][i])
  plot(c(x, true_val[1]), c(y, true_val[2]), pch=19,col="blue",cex=0.8, main="Majorize-Minimize", xlab="", ylab="", ylim=c(0, 5))
  lines(x, y, col="blue")
  abline(v=0,col="red") 
  dev.off()
}

####
MM.RSS <- MM$objfctList+sum(v^2)
MM.RSS.min <- sum( (v-W%*%MM$x)^2 )
MM.RelLoss <- (MM.RSS-MM.RSS.min)/MM.RSS.min
plot(MM.RSS,type="l",col="black",xlab="Iteration",ylab="RSS")
plot(MM$time,MM.RSS,type='l',col='black',lwd=1,xlab='Time',ylab='RSS')
abline(h=MM.RSS.min, col="red")


### Test problem outside NLS ###

kidsdata <- read.csv("http://www.murraylax.org/datasets/gradeschool.csv")
grades <- kidsdata$Grades
median(grades)

MajorizeMinimizeMedian <- function(y,x,maxIter=10000,tol=0.00001){
  ## Objective function
  f = function(x) as.numeric(sum(abs(y-x)))
  ## Initial values of xList, objFctList and time
  xList = x
  objfctList = c(f(x))
  time = c(0)
  start = Sys.time()
  ##-------------------------------------------
  ## Iterative procedure for solving the NQP
  ##-------------------------------------------
  ## Minimization 
  for (i in 1:maxIter){
    #x = sum(y/(abs(y-x)))/sum(abs(y-x)) # update function
    w = 1/abs(y-x)
    x = sum(w*y)/sum(w)
    objfct = f(x)
    if (abs(objfct - objfctList[i]) < tol) break
    xList = cbind(xList, x)
    objfctList = c(objfctList, objfct)
    time = c(time, Sys.time() - start)
  } 
  return (list(x=x, objfct=objfct, xList = xList, objfctList = objfctList, time=time))
}

mm <- MajorizeMinimizeMedian(grades, 0)
mm
mm$x
step <- 1
# Illustrate step by step; run from here to line 114 multiple times for multiple steps
# Don't mind warnings (some results are NA)
rang <- seq(0, max(grades), length.out = 50)
init_theta <- mm$xList[1,]
thetas <- rang
ys <- c() ; gs <- c()
for (i in 1:length(thetas)) {
    ys[i] <- sum(abs(grades - thetas[i]))
    gs[i] <- 1/2*sum((grades-thetas[i])^2/abs(grades-init_theta[step]) + abs(grades-init_theta[step]))
}

df <- data.frame(thet = thetas, f = ys, g = gs) 

p <- df %>% ggplot() +
    geom_bar(data = kidsdata, aes(x = grades)) +
    geom_line(aes(x = thet, y = f, color = "f")) +
    geom_line(aes(x = thet, y = g, color = "g")) +
    geom_vline(aes(xintercept = init_theta[step])) +
    ylim(c(0, 1500)) + ylab("y") + ggtitle(paste0('Step ', step)) +
    annotate("text", x = init_theta[step] + 0.2, y = 1400, label = expression('h'^'t'))

ggsave(paste0('img/test',step,'.png'), p)
p <- p + geom_vline(aes(xintercept = init_theta[step+1])) +
    annotate("text", x = init_theta[step+1] + 0.2, y = 1300, label = expression('h'^'t+1'))
ggsave(paste0('img/test',step,'min.png'), p)
step = step + 1


xList <- mm$xList
xList
true_val <- median(y)

mm_kids <- MajorizeMinimizeMedian(grades, 10)
xList <- mm_kids$xList
true_val <- median(grades)

y_axis <- c()
x_axis <- c()
for(i in 1:length(xList)){
  jpeg(paste0('img/kids_ex',i,'.jpg'))
  x_axis <- c(x_axis, i)
  y_axis <- c(y_axis, xList[i])
  plot(c(x_axis, 50), c(y_axis, true_val), pch=19,col="blue",cex=0.8, main="Finding a Sample Median", xlab="", ylab="", ylim=c(0, 50), xlim=c(0,50))
  lines(x_axis, y_axis, col="blue")
  abline(v=0,col="red") 
  dev.off()
}

###
hist(grades)
median(grades)
mm_kids <- MajorizeMinimizeMedian(grades, 0)
mm_kids$x
mm_kids
