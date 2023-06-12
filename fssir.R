#############################################################
##  Adjusting high-dimensional sparse SIR for clustered X: ##
##  By Wei Luo and Yan Guo (2023.5.29)                     ##
##  Submitted to EJS (Electronic Journal of Statistics)    ##
##                                                         ##
##  The main function "ssir" is provided by Tan (2018),    ##
##  for algorithm details, see their paper                 ##
##  "A convex formulation for high-dimensional             ##
##   sparse sliced inverse regression "                    ##
#############################################################


################################################
##   Main function for Convex Sparse SIR      ##
##     M - the kernel matrix of SDR methods   ##
##         function covxy - for SSIR          ##
##         function covxy_M - for SSIR_M      ##
##         function covxy_RM - for SSIR_RM    ##
##     Covx - marginal covariance of x        ##
##     Lambda - sparsity tuning parameter     ##
##     K - structural dimension               ##
################################################
ssir <- function(M,covx,lambda,K,nu=1,epsilon=1e-3,maxiter=1000,trace=FALSE,init=FALSE,initPi=NULL,initH=NULL,initGamma=NULL){
  p <- nrow(covx)
  eigencovx <- eigen(covx)
  sqcovx <- eigencovx$vectors%*%sqrt(diag(pmax(eigencovx$values,0)))%*%t(eigencovx$vectors)	
  covxy <- M
  tau <- 4*nu*eigencovx$values[1]^2	
  criteria <- 1e10
  i <- 1
  
  # Initialize parameters
  H <- Pi <- oldPi <-  diag(1,p,p)
  Gamma <- matrix(0,p,p)
  
  if(init==TRUE){
    H <- initH
    Pi <- initPi
    Gamma <- initGamma
  }
  
  # While loop for the iterations
  while(criteria > epsilon && i <= maxiter){
    Pi <- updatePi(covx,sqcovx,covxy,H,Gamma,nu,lambda,Pi,tau)
    
    H <- updateH(sqcovx,Gamma,nu,Pi,K)
    Gamma <- Gamma + sqcovx%*%Pi%*%sqcovx-H	
    criteria <- sqrt(sum((Pi-oldPi)^2))
    oldPi <- Pi
    i <- i+1
    if(trace==TRUE)
    {
      print(i)
      print(criteria)
    }
  }
  return(list(Pi=Pi,H=H,Gamma=Gamma,iteration=i,convergence=criteria))
}



######################################################
##                 Update Pi                        ##
######################################################
updatePi <- function(covx,sqcovx,covxy,H,Gamma,nu,lambda,Pi,tau){
  
  A <- Pi + 1/tau*covxy-nu/tau*covx%*%Pi%*%covx+nu/tau*sqcovx%*%(H-Gamma)%*%sqcovx
  B <- lambda/tau
  return(Soft(A,B))
}


######################################################
##                Update H                          ## 
######################################################
updateH <- function(sqcovx,Gamma,nu,Pi,K){
  
  temp <- Gamma + sqcovx%*%Pi%*%sqcovx
  temp <- (temp+t(temp))/2
  svdtemp <- eigen(temp)
  d <- svdtemp$values
  p <- length(d)
  
  if(sum(pmin(1,pmax(d,0)))<=K){
    dfinal <- pmin(1,pmax(d,0))
    return(svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors))
  }
  
  fr <- function(x){
    sum(pmin(1,pmax(d-x,0)))
  }
  # Vincent Vu Fantope Projection
  knots <- unique(c((d-1),d))
  knots <- sort(knots,decreasing=TRUE)
  temp <- which(sapply(knots,fr)<=K)
  lentemp <- tail(temp,1)
  a=knots[lentemp]
  b=knots[lentemp+1]
  fa <- sum(pmin(pmax(d-a,0),1))
  fb <- sum(pmin(pmax(d-b,0),1))
  theta <- a+ (b-a)*(K-fa)/(fb-fa)
  dfinal <- pmin(1,pmax(d-theta,0))
  res <- svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors)
  return(res)
}	



######################################################
##          Soft-thresholding Operator              ##
######################################################
Soft <- function(a,b){
  if(b<0) stop("Can soft-threshold by a nonnegative quantity only.")
  return(sign(a)*pmax(0,abs(a)-b))
}


##########################################################
###    slicing function: discretizing y into h slices  ###
###    Input: y - n-dim vector                         ###
###           h - the number of slices                 ###
###   Output: ys - n*h matrix                          ###
##########################################################
slicing<-function(y,H)
{
  if (length(levels(as.factor(y)))>H)
  {
    ytilde<-rep(0,H+1)
    ytilde[1]<-min(y)
    for (h in 1:(H-1))
    {
      ytilde[h+1]<-quantile(y,h/H)
    }
  }
  if (length(levels(as.factor(y)))<=H)
  {
    H <- length(levels(as.factor(y)))
    ytilde<-rep(0,H+1)
    ytilde[1]=min(y)
    for (h in 1:(H-1))
    {
      ytilde[h+1]<-min(y[y>ytilde[h]])
    }
  }
  ytilde[H+1]=max(y)+1
  prop<-rep(1,H)
  for (i in 1:H)
  {
    prop[i] = sum((y >= ytilde[i])&(y < ytilde[i+1]))/n
  }
  ys <- matrix(0,n,H)
  for (i in 1:n)
  {
    for (j in 1:H)
    {
      if ((y[i] >= ytilde[j])&(y[i] < ytilde[j+1]))
      {
        ys[i,j] <- 1
      }
    }
  }
  res<-list()
  res$H<-H
  res$ytilde<-ytilde
  res$prop<-prop
  res$ys<-ys
  return(res)
}

#############################################################################
##    A function to calculate the conditional covariance for SSIR          ##
##     If y is categorical, nslice is # of unique values in y              ##
##     If y is continuous, then the conditional covariance is constructed  ##
##        with Sigma_x - T where T is constructed based on Equation (7)    ##
##     X is first mean centered                                            ##
#############################################################################
covxy <- function(y,X,nslice,standardize=TRUE){
  n <- length(y)
  if(standardize==TRUE){ 	
    X <- scale(X,T,F)
  }
  
  # Test if this is integer or continuous
  if(sum(y-round(y))==0){
    nslice <- length(unique(y))
    quany <- sort(unique(y))
  }
  
  # For continuous response y
  else if(sum(y-round(y))!=0){
    quany <- quantile(y,seq(1/nslice,1,length.out=nslice))
  }
  
  indexy <- vector("list",nslice)
  
  # Assign indices into indexy
  indexy[[1]] <- which(y<=quany[1])
  for(k in 2:nslice){
    indexy[[k]] <- which(y >quany[k-1] & y<=quany[k])
  }
  
  nindexy <- lapply(indexy,length)
  
  # f matrix 
  f <- matrix(0,n,nslice)
  for(k1 in 1:(nslice-1)){
    for(k2 in 1:nslice){
      if(k1==k2){
        f[indexy[[k1]],k2] <- 1 - nindexy[[k2]]/n
      }
      if(k1!=k2){
        f[indexy[[k1]],k2] <- -nindexy[[k2]]/n			
      }
    }
  }
  for(k in 1:nslice){
    f[indexy[[nslice]],k] <- -nindexy[[k]]/n
  }
  
  bigF <- f%*%solve(t(f)%*%f)%*%t(f)
  Sigmafit <- t(X)%*%bigF%*%X/(n)
  return(Sigmafit)
}




####################################################################
##    A function to calculate the kernel matrix for SSIR_M        ##
####################################################################
covxy_M <- function(x,y_s,h,q,pix,mux)
{
  p<-ncol(x)
  n<-nrow(x)
  newy <- y_s
  M1 <- matrix(0,p,q*h)
  for (k in 1:q)
  {
    M1[,(h*(k-1)+1):(k*h)]<-((t(x)-mux[,k])%*%(newy*pix[,k]))/(n/h)
  }  
  return(M1)
}


####################################################################
##    A function to calculate the kernel matrix for SSIR_RM       ##
####################################################################
covxy_RM <- function(x,y_s,h,q,beta,pix,mux,sigx)
{
  p<-ncol(x)
  n<-nrow(x)
  rx<-x%*%beta
  d<-ncol(rx)
  newy <- y_s
  
  mubx <- list(as.numeric(mux[,1]%*%beta),as.numeric(mux[,2]%*%beta))
  sigbx <- t(beta)%*%sigx%*%beta
  yy <- as.vector(pix[,1])
  xx <- as.vector(rx)
  bw1<-npregbw(yy~xx,bwscaling=FALSE,ckertype ="gaussian")$bw
  fit<-npreg(yy~xx,bws=bw1,residuals =TRUE,ckertype="gaussian")
  pre <- predict(fit)
  pibx <- cbind(pre,1-pre)
 
  u <- matrix(0,p,d*h)
  for (i in 1:n)
  {
    ### ed
    ed <- matrix(0,d,h)
    for (j in 1:q)
      for (k in 1:q)
      {
        ed <- ed + matpower(sigbx,-1)%*%
          (t(rx)-mubx[[k]])%*%(newy*pibx[,j]*pibx[,k]/(n/h))*pibx[i,j]
      }
    xi <- x[i,]
    for (j in 1:q)
      xi <- xi-(pibx[i,j]*mux[,j])
    u <- u + xi %*% t( as.vector( ed%*%diag(newy[i,])) )
  }
  return(u/(n/h))
}





