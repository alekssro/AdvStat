######################################################
## All the basic algorithms for the NQP problem:
## MM 
## PCD
## PGD
## MuEM
## FCEM
## CNP
######################################################
##	MAJORIZE MINIMIZE ALGORITHM FOR THE NQP PROBLEM ##
######################################################
##
## Date: June 2018
## Authors: Hobolth, Guo, Kousholt, Jensen
##
## The algorithm minimizes the quadratic form 
## f(x) = 1/2x'Ax + b'x subject to x >= 0 
## for a symmetric positive definite matrix A and 
## a negative vector b. 
## 
## The solution to the 
## non-negative quadratic programming problem (NQP) is
## based on the MM algorithm (Lee and Seung, 1999)
## 
## Name: MajorizeMinimizeNQP
## 
## Input: 
## A: positive definite K by K matrix with positive entries.
## b: negative vector of length K. 
## init: initial value of x, vector of length K. Optional. 
## maxIter: maximum number of iterations.
## tol: tolerance. Algorithm breaks when the objective function
## decreases tol.
##
## Output:
## x: minimizer of objective function f.
## objfct: final value of objective function f(x).
## xList: values of x in each iteration.
## objfctList: values of f(x) in each iteration.
## time: accumulated time after each iteration. 
##
###############################################################
MajorizeMinimizeNQP <- function(A,b,init,maxIter=1000,tol=0.001){
  ## Objective function
  f = function(x) as.numeric(1/2*crossprod(x, A)%*%x+crossprod(b, x))
  ## Initial values of xList, objFctList and time
  K = ncol(A)
  if (missing(init)) {
    x = matrix(runif(K), nrow=K, byrow=T)
  } else {
    x = init
  }
  xList = x
  objfctList = c(f(x))
  time = c(0)
  start = Sys.time()
  ##-------------------------------------------
  ## Iterative procedure for solving the NQP
  ##-------------------------------------------
  ## Minimization 
  for (i in 1:maxIter){
    x = -b*x/(A %*% x) 
    objfct = f(x)
    if (abs(objfct - objfctList[i]) < tol) break
    xList = cbind(xList, x)
    objfctList = c(objfctList, objfct)
    time = c(time, Sys.time() - start)
  } 
  return (list(x=x, objfct=objfct, xList = xList, objfctList = objfctList, time=time))
}
#####################################################
## COORDINATE DESCENT ALGORITHM FOR THE NQP PROBLEM 
#####################################################
##
## Date: June 2018
## Authors: Hobolth, Guo, Kousholt, Jensen
## 
## The algorithm minimizes the quadratic form 
## f(x) = 1/2x'Ax + b'x subject to x >= 0 
## for a symmetric positive definite matrix A 
## and a negative vector b.
## 
## The solution to the 
## non-negative quadratic programming problem (NQP)
## is based on coordinate descent (CD) algorithm
## 
## Name: CoordinateDescentNQP
##
## Input: 
## A: positive definite K by K matrix with positive entries.
## b: negative vector of length K. 
## init: initial value of x, vector of length K. Optional. 
## maxIter: maximum number of iteration
## tol: tolerance. Algorithm breaks when objective function 
##      decreases less than tol.      		
##
## Output:
## x: minimizer of objective function f.
## objfct: final value of objective function f(x).
## xList: values of x in each iteration.
## objfctList: values of f(x) in each iteration.
## time: accumulated time after each iteration. 
##
###########################################################################
CoordinateDescentNQP = function(A, b, init, maxIter=1000, tol=0.001) {
  ## Objective function
  f = function(x) as.numeric(1/2*crossprod(x, A)%*%x+crossprod(b, x))
  ## Initial value
  K = ncol(A)
  if (missing(init)) {
    x = matrix(runif(K),nrow=K,byrow=T)
  } else {
    x = init
  }
  xList = x
  objfctList = c(f(x))
  time = c(0)
  start = Sys.time()
  ##---------------------------------------------
  ## Iterative procedure for solving the NQP
  ##---------------------------------------------
  ## Minimization  
  for (i in 1:maxIter) {
    for (j in 1:K) {
      grdj = b[j] + crossprod(A[,j], x) 
      x[j] = max(x[j] - grdj / A[j,j],0)
    } ## end for j
    objfct = f(x)
    if (abs(objfct - objfctList[i]) < tol) break 
    xList = cbind(xList, x)
    objfctList = c(objfctList, objfct)
    time = c(time, Sys.time() - start)    
  } ## end for i
  return (list(x=x,objfct=objfct,xList=xList,objfctList=objfctList,time=time))
}
#################################################################
##	PROJECTED GRADIENT DESCENT ALGORITHM FOR THE NQP PROBLEM 
#################################################################
##
## Date: June 2018
## Authors: Hobolth, Guo, Kousholt, Jensen
##
## The algorithm minimizes the quadratic form 
## f(x) = 1/2x'Ax + b'x subject to x >= 0 
## for a symmetric positive definite matrix A 
## and a negative vector b. 
## 
## The solution to the 
## non-negative quadratic programming problem (NQP)
## is based on the projected gradient descent (PGD) algorithm
## Two choices are possible:
## a. Exact line search
## b. Stepsize r/L where L is max eigenvalue of A 
##    and r is a value between 0 and 2.
##    (see Lange, Chi and Zhou, 2014, page 50)
##
## Name: ProjectedGradientDescentNQP
##
## Input: 
## A: positive definite K by K matrix with positive entries.
## b: negtaive vector of length K. 
## init: initial value of x, vector of length K. Optional argument. 
## maxIter: maximum number of iterations.
## tol: tolerance. Algorithm breaks when the objective function 
##      decreases less than tol.      		
## exactLineSearch: Logical. If TRUE, the (unconstrained) optimal 
##                  stepsize is calculated and used in each iteration.
## r: If exactLineSearch = FALSE, then stepsize is r/L 
##   (where L is the max eigenvalue of A)
##
## Output:
## x: minimizer of objective function f.
## objfct: final value of objective function f(x).
## xList: values of x in each iteration.
## objfctList: values of f(x) in each iteration.
## time: accumulated time after each iteration. 
##
#########################################################################
ProjectedGradientDescentNQP = 
  function(A,b,init,maxIter=1000,tol=0.001,exactLineSearch=F,r=1){
    ## Objective function
    f = function(x) as.numeric(1/2*crossprod(x, A) %*% x + crossprod(b, x))
    ## Initial value
    K = ncol(A)
    if (missing(init)) {
      x = matrix(runif(K), nrow=K, byrow=T)
    } else {
      x = init
    }
    xList = x
    objfctList = c(f(x))
    time = c(0)
    start = Sys.time()
    ## Gradient function
    Grad = function(x){
      return (crossprod(A, x) + b)
    }
    ## Exact line search for s
    if (exactLineSearch){	
      Eta = function(grad) norm(grad,type='F')^2 / as.numeric((crossprod(grad, A) %*% grad))
    }
    ## Lange, Chi and Zhou (2014): s=r/L 
    if (exactLineSearch==FALSE){
      L = max(eigen(A)$values)
      s = r / L
    }
    ##-----------------------------------------------------
    ## Iterative procedure for solving the NQP
    ##-----------------------------------------------------
    ## Minimization
    for (i in 1:maxIter) {
      g = Grad(x)
      if (exactLineSearch) {
        s = Eta(g)
      }
      x = x - s * g
      x[x < 0] = 0
      objfct = f(x)
      if (abs(objfct - objfctList[i]) < tol) break
      xList = cbind(xList, x)
      objfctList = c(objfctList, objfct)
      time = c(time, Sys.time() - start)
    } #end for i
    return(list(x=x,objfct=objfct,xList=xList,objfctList=objfctList,time=time))
  }
#####################################################
## Multiplicative EM ALGORITHM FOR THE NQP PROBLEM
#####################################################
##
## Date: June 2018
## Authors: Hobolth, Guo, Kousholt, Jensen
##
## The algorithm minimizes the quadratic form 
## f(x) = 1/2x'Ax + b'x subject to x >= 0 
## for a symmetric positive definite matrix A 
## and a negative vector b.
##
## The solution to the 
## non-negative quadratic programming problem (NQP) is
## based on the EM algorithm with multiplicative updates.
##
## Name: MultEMNQP
##
## Input: 
## A: positive definite K by K matrix
## b: negtaive vector of length K 
## init: initial value of x, vector of length K. Optional argument. 
## maxIter: maximum number of iteration
## tol: tolerance. Algorithm breaks when the objective function 
## decreases less than tol.      		
##
## Output:
## x: minimizer of objective function f.
## objfct: value of objective function in x, f(x).
## xList: values of x in each iteration
## objfctList: values of f(x) in each iteration.
## time: accumulated time after each iteration. 
##
#########################################################################
MultEMNQP = function(A, b, init, maxIter=1000, tol=0.001){
  ## Objective function
  f = function(x) as.numeric(1/2*crossprod(x, A)%*%x+crossprod(b,x))
  ## Initial value
  K=nrow(A)
  if (missing(init)) {
    x = matrix(runif(K), nrow=K, byrow=T)
  } else {
    x = init
  }
  xList = x
  objfctList = c(f(x))
  time = c(0)
  start=Sys.time()
  ##---------------------------------------------------------------------
  ## Iterative procedure of NQP
  ##---------------------------------------------------------------------
  for (i in 1:maxIter){
    sigmaSquare=as.numeric(crossprod(diag(A),x))
    grad=crossprod(t(A),x) + b      		
    tauSquare=ifelse( max(grad)>sigmaSquare , max(grad)-sigmaSquare , 0)		 
    x=x*(1-grad/(tauSquare+sigmaSquare))
    objfct = f(x)
    if (abs(objfct - objfctList[i]) < tol) break 
    xList = cbind(xList, x)
    objfctList = c(objfctList, objfct)
    time = c(time, Sys.time() - start)    
  } 
  return(list(x=x,objfct=objfct,xList=xList,objfctList=objfctList,time=time))
}
########################################################
## FEVOTTE AND CEMGIL EM ALGORITHM FOR THE NQP PROBLEM
########################################################
## 
## Date: June 2018
## Authors: Hobolth, Guo, Kousholt, Jensen
##
## The algorithm minimizes a quadratic form 
## f(x) = 1/2x'Ax + b'x subject to x >= 0 for 
## a symmetric positive definite matrix A 
## and a negative vector b.
## 
## The solution to the 
## non-negative quadratic programming problem (NQP) is
## based on the specific EM algorithm described in 
## Fevotte and Cemgil (2009).
## 
## Name: FevotteCemgilNQP
## 
## Input: 
## A: positive definite K by K matrix
## b: negative vector of length K 
## init: initial value of x, vector of length K. Optional argument. 
## maxIter: maximum number of iteration
## tol: tolerance. Algorithm breaks when the objective function 
## decreases less than tol.      		
##
## Output:
## x: minimizer of objective function f.
## objfct: value of objective function in x, f(x).
## xList: values of x in each iteration
## objfctList: values of f(x) in each iteration.
## time: accumulated time in each iteration step. 
##
###################################################################################
FevotteCemgilNQP = function(A, b, init, maxIter=1000, tol=0.001){
  ## Objective function
  f = function(x) as.numeric(1/2*crossprod(x, A) %*% x + crossprod(b, x))
  ## Initial value
  K=nrow(A)
  if (missing(init)) {
    x = matrix(runif(K), nrow=K, byrow=T)
  } else {
    x = init
  }
  xList = x
  objfctList = c(f(x))
  time = c(0)
  start = Sys.time()
  ##--------------------------------------------
  ## Iterative procedure for solving the NQP
  ##--------------------------------------------
  ## Minimization
  for (i in  (1:maxIter)){
    grad <- A %*%x+b
    x=x-grad/diag(A)/K
    x=pmax(x,0)
    objfct = f(x)
    if (abs(objfct - objfctList[i]) < tol) break
    xList = cbind(xList, x)     
    objfctList = c(objfctList, objfct)
    time = c(time, Sys.time() - start)     
  } 
  return(list(x=x,objfct=objfct,xList=xList,objfctList=objfctList,time=time))
}
#############################################
##  CONE PROJECTION FOR THE NLS PROBLEM    ##
#############################################
##
## Date: June 2018
## Authors: Hobolth, Guo, Kousholt, Jensen
##
## The algorithm finds the non-negative vector hhat
## that minimizes the residual sum of squares
## RSS(h)=(v-W h)'(v-W h)
##
## The solution to the problem is based on a modification
## of the cone projection algorithm in Meyer (2013).
##
## Name: ConeProjNonNegRegr
## 
## Input:
## v: vector of length n with positive entries
## W: n times k matrix with positive entries
## 
## Output:
## x: minimizer of the RSS 
##    (entries constrained to be non-negative)
## objfct: final value of the RSS
## time: time for the projection
##
##############################################################
ConeProjNonNegRegr <- function(v,W){
  tol=1e-8 ; n=length(v) ; k=length(W)/n
  start=Sys.time()
  ##-------------------
  ## Initial step
  ##-------------------
  a=t(W)%*%v
  if(max(a)<tol){
    h=rep(0,k) ; chck=0    ## Finish and return h=0
  }else{
    J0=which.max(a) ; chck=1
  }
  while(chck==1){
    ##------------------
    ## Step 1
    ##------------------
    W0 = W[,J0,drop=FALSE]
    a = solve( t(W0)%*%W0 )%*%t(W0)%*%v
    cc = t(W)%*%(v-W0%*%a)
    if(max(cc)<tol){
      h=rep(0,k) ; h[J0]=a ; chck=0  ## Finish and return h=a
    }else{
      j1=which.max(cc) ; J1=sort(c(J0,j1))
      hw=rep(0,k) ; hw[J0]=a ; h1=hw[J1]
      chck=2
      while(chck==2){
        ##------------------
        ## Step 2
        ##------------------
        W1=W[,J1,drop=FALSE]
        d1=solve( t(W1)%*%W1 )%*%t(W1)%*%v
        if(min(d1)<(-tol)){
          xhat=min( h1[d1<0]/(h1[d1<0]-d1[d1<0]) )
          h1hat=xhat*d1+(1-xhat)*h1
          J1=J1[h1hat>tol]
          h1=h1hat[h1hat>tol]
        }else{
          J0=J1[d1>tol] ; chck=1
        }
      }
    }
  }
  ##---------------------------------
  time=Sys.time()-start
  RSS=sum( (v-W%*%h)^2 )
  return(list(x=h,objfct=RSS,time=time))
} 