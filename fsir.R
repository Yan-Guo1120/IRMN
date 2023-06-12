##########################################################
##  Adjusting SIR for clustered X:                      ##
##     SIR_M  (M in the suffix for the mixture model)   ##
##     SIR_RM (R in the suffix for refined)             ##
##  By Wei Luo and Yan Guo (2023.5.29)                  ##
##  Submitted to EJS (Electronic Journal of Statistics) ##
##########################################################


#######################################################
##   Some basic functions：                         ###
##    1. power of matrix (matpower)                 ###
##    2. the projection matrix (proj)               ###
##    3. the distance between two spaces (dv)       ###
##    4. discretize Y into slices for SDR (ys)      ###
##    5. mixture parameter for X (mnEM)             ###
##    6. mixture parameter for beta*X (parabx)      ###
##   Main functions:                                ###
##    1. SIR_M                                      ###
##    2. SIR_RM                                     ###
#######################################################


#######################################
###    power of a matrix            ###
###    Input: a - original matrix   ###
###           alpha - the power     ###
###   Output: new matrix            ###
#######################################
matpower<-function(a,alpha){
      small <- .00000001
a<-as.matrix(a)
p1<-nrow(a)
if (p1==1) {return(matrix(c(a)^alpha,1,1))} else {
  eva<-eigen(a)$values
  eve<-eigen(a)$vectors
  eve<-eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
  index<-(1:p1)[eva>small]
  evai<-eva
  evai[index]<-(eva[index])^(alpha)
  ai<-eve%*%diag(evai)%*%t(eve)
  return(ai)
  }
}

#########################################################################
###    the projection matrix                                          ###
###    Input: u - p*d matrix                                          ###
###           Lam - p*p positive definite matrix                      ###
###   Output: the projection matrix of u under the Lam inner product  ###
#########################################################################
proj<-function(u,Lam)
 {
   if (is.vector(u)==TRUE) {return(u%*%t(u)%*%Lam/c(t(u)%*%Lam%*%u))} else {
   if (ncol(u)==1) {return(u%*%t(u)%*%Lam/c(t(u)%*%Lam%*%u))} else {
                    return(u%*%matpower(t(u)%*%Lam%*%u,-1)%*%t(u)%*%Lam)}
      }
 }

##########################################################
###    the distance between two spaces                 ###
###    Input: u - p*d basis matrix                     ###
###           v - p*d basis matrix                     ###
###   Output: the distance between Span(u) and Span(v) ###
##########################################################
dv<-function(u,v)
 {
if (is.vector(u)==TRUE) {p=length(u)} else {p=nrow(u)}
   return(sqrt(sum((proj(u,diag(1,p))-proj(v,diag(1,p)))^2)))
}

##########################################################
###    slicing function: discretizing y into h slices  ###
###    Input: x - n*p matrix                           ###
###           y - n-dim vector                         ###
###           h - the number of slices                 ###
###   Output: n*h matrix                               ###
##########################################################
ys <- function(y,x,h)
{
  n<-nrow(x)
  yx <- data.frame(y,x)
  beta<-dr(y~x,data=yx,nslices=h)           
  newy<-matrix(0,n,h)
  for (i in 1:n)
    for (j in 1:h)
    {
      if (beta$slice.info$`slice.indicator`[i]==j) 
        newy[i,j]<-1
      else
        newy[i,j]<-0
    }
  return(newy)
}

################################################################################
####  mnEM function                                                          ###
##    (EM algorithm to estimate parameters for mixture normal distribution)  ###
##    Input: x - nxp matrix                                                  ###
##           q - number of clusters                                          ###
##    Output: mux - a list of the mean vectors                               ###
##            sigx - a list of the covariance matrices                       ###
##            pix - n*q matrix of posterior probabilities                    ###
##    Remark: we use Kmeans to give an initial value for EM algorithm        ###
################################################################################
mnEM <- function(x,q)
{
  s<-kmeans(x,centers=q)
  lam0<-s$size/n
  mu0<-list()
  for (k in 1:q)
  {mu0[[k]]=s$centers[k,]}
  sx<-try(mvnormalmixEM(x,k=q,lambda=lam0,mu=mu0,maxit=100,arbmean=TRUE),silent=TRUE)
  if('try-error'%in%class(sx))next
  mux<-sx$mu
  sigx<-sx$sigma
  pix<-sx$posterior
  res <- list()
  res$mux<-mux
  res$sigx<-sigx
  res$pix<-pix
  return(res)
}

################################################################################
##    parabx function                                                        ###
##    (parameters for mixture normal distribution - beta*X)                  ###
##    Input: d - structural dimension                                        ###
##           q - number of clusters                                          ###
##           nmlist - parameters for mixture normal distribution - X         ###
##           beta - basis of central subspace                                ###
##    Output: mubx - a list of the mean vectors                              ###
##            sigbx - a list of the covariance matrices                      ###
##            pibx - n*q matrix of posterior probabilities                   ###
##    Remark: we use nonlinear estimation for pibx:                          ###
##            regress pi(beta*X) on pi(X)                                    ###
################################################################################
parabx <- function(d,q,nmlist,beta)
{
  mux <- nmlist$mux
  sigx <- nmlist$sigx
  pix <- nmlist$pix
  
  mubx <- list()
  sigbx <- list()
  for (i in 1:q)
  {
    mubx[[i]] <- as.numeric(mux[[i]]%*%beta)
    sigbx[[i]] <- t(beta)%*%sigx[[i]]%*%beta
  }
  
  pibx_np <- numeric(0)
  pibx_sum <- 0
  for (i in 1:(q-1))
  {
    yy <- as.vector(pix[,i])
    xx <- as.matrix(x%*%beta)
    if (d==1)
    {
      bw1<-npregbw(yy~xx[,1],bwscaling=FALSE,ckertype ="gaussian")$bw
      fit<-npreg(yy~xx[,1],bws=bw1,residuals =TRUE,ckertype="gaussian")
    }else
      if (d==2)
      {
        bw1<-npregbw(yy~xx[,1]+xx[,2],bwscaling=FALSE,ckertype ="gaussian")$bw
        fit<-npreg(yy~xx[,1]+xx[,2],bws=bw1,residuals =TRUE,ckertype="gaussian")
      }
    pre <- predict(fit)
    pibx_np <- cbind(pibx_np, pre)
    pibx_sum <- pibx_sum+pre
  }
  pibx <- cbind(pibx_np,1-pibx_sum)
  res <- list()
  res$mubx<-mubx
  res$sigbx<-sigbx
  res$pibx<-pibx
  return(res)
}


#######################################################################
##    SIR_M function                                                ###
##    (SIR for predictors with mixture normal distribution)         ###
##    Input: x - predictor, nxp matrix                              ###
##           y - response, n-dim vector                             ###
##           q - number of clusters                                 ###
##           d - structural dimension                               ###
##           h - number of slices                                   ###
##   Output: beta - estimated basis function for CS                 ###
##           para - a list of parameters (mux, sigx, pix)           ###
##                       for mixture normal distribution            ###
#######################################################################
SIR_M<-function(x,y,q,d,h)
{
  p<-ncol(x)
  n<-nrow(x)
  newy <- ys(y,x,h)
  normalmix <- mnEM(x,q)
  mux <- normalmix$mux
  sigx <- normalmix$sigx
  pix <- normalmix$pix
  mols<-matrix(0,p,q*h)
  for (k in 1:q)
  {
    mols[,(h*(k-1)+1):(k*h)]<-matpower(sigx[[k]],-1)%*%((t(x)-mux[[k]])%*%(newy*pix[,k]))/n
  }  
  b<-svd(mols)$u[,1:d]
  res<-list()
  res$para<-normalmix
  res$beta<-b
  return(res)
}


#################################################################################
##    SIR_RM function                                                         ###
##    (Revised SIR for the predictors with mixture normal distribution)       ###
##    Input: x - n*p matrix                                                   ###
##           y - n-dim vector                                                 ###
##           q - number of clusters                                           ###
##           d - structural dimension                                         ###
##           h - number of slices                                             ###
##           beta - the initial value (we take SIR_M as the initial value)    ###
##           nmlist - the parameters of mixture normal distribution           ###
##   Output: beta - estimated basis function for the central subspace         ###
#################################################################################
SIR_RM<-function(x,y,q,d,h,beta,nmlist)
{
  p<-ncol(x)
  n<-nrow(x)
  rx<-x%*%beta
  newy <- ys(y,x,h)
  
  mux <- nmlist$mux
  sigx <- nmlist$sigx
  pix <- nmlist$pix
  par <- parabx(d,q,nmlist,beta)
  mubx <- par$mubx
  sigbx <- par$sigbx
  pibx <- par$pibx
  
  dis<-1   #### The distance between the estimators obtained by the two iterations, 
           #### Set the initial value of dis as 1, and the subsequent value is given by the following iteration process. 
  e<-1e-4  #### A threshold. If dis<e, stop the iteration.
  m<-1     #### The number of iterations.
  betak<-beta
  while(dis>e)
  {
    ### u
    u <- matrix(0,p,d*h)
    for (i in 1:n)
    {
      ### ed
      ed <- matrix(0,d,h)
      for (j in 1:q)
        for (k in 1:q)
        {
          ed <- ed + matpower(sigbx[[k]],-1)%*%
            (t(rx)-mubx[[k]])%*%(newy*pibx[,j]*pibx[,k]/(n/h))*pibx[i,j]
        }
      xi <- x[i,]
      for (j in 1:q)
        xi <- xi-(pibx[i,j]*mux[[j]])
      u <- u + xi %*% t( as.vector( ed%*%diag(newy[i,])) )
    }
    u <- as.vector(u/(n/h))
    
    ### v
    v <- matrix(0,p*d*h,p*d)
    for (i in 1:n)
    {
      ### ed
      ed <- matrix(0,d,h)
      for (j in 1:q)
        for (k in 1:q)
        {
          ed <- ed + matpower(sigbx[[k]],-1)%*%
            (t(rx)-mubx[[k]])%*%(newy*pibx[,j]*pibx[,k]/(n/h))*pibx[i,j]
        }
      for (j in 1:q)
      {
        bk<-t(betak)%*%sigx[[j]]%*%betak
        bt<-t(t(as.vector(ed%*%diag(newy[i,]))))%*%(x[i,]-mux[[j]])%*%betak%*%t(matpower(bk,-1))
        a<-pibx[i,j]*sigx[[j]]
        v<-v+kronecker(bt,a)
      }
    }
    v <- v/(n/h)
    
    betak1<- matpower(t(v)%*%v,-1)%*%t(v)%*%u
    betak1<-matrix(betak1,p,d)
    betak1 <- betak1/sqrt(sum(betak1^2))
    dis<-dv(betak,betak1)
    betak<-betak1
    m <- m+1
    if(m>5)
      break
  }
  res<-list()
  res$beta<-betak
  return(res)
}










