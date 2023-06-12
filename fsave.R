##########################################################
##  Adjusting SAVE for clustered X:                     ##
##     SAVE_M  (M in the suffix for the mixture model)  ##
##     SAVE_RM (R in the suffix for refined)            ##
##  By Wei Luo and Yan Guo (2023.5.29)                  ##
##  Submitted to EJS (Electronic Journal of Statistics) ##
##########################################################


#######################################################
##   Some basic functions：                         ###
##    1. power of matrix (matpower)                 ###
##    2. the projection matrix (proj)               ###
##    3. the distance between two spaces (dv)       ###
##    4. discretize Y into slices for SDR (slicing) ###
##    5. mixture parameter for X (mnEM)             ###
##    6. mixture parameter for beta*X (parabx)      ###
##   Main functions:                                ###
##    1. SAVE_M                                     ###
##    2. SAVE_RM                                    ###
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
###    Input: y - n-dim vector                         ###
###           h - the number of slices                 ###
###   Output: ytilde - cut off points for slicing      ###
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
  res<-list()
  res$H<-H
  res$ytilde<-ytilde
  res$prop<-prop
  return(res)
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
##    SAVE_M function                                               ###
##    (SAVE for predictors with mixture normal distribution)        ###
##    Input: x - predictor, nxp matrix                              ###
##           y - response, n-dim vector                             ###
##           q - number of clusters                                 ###
##           d - structural dimension                               ###
##           h - number of slices                                   ###
##   Output: beta - estimated basis function for CS                 ###
##           para - a list of parameters (mux, sigx, pix)           ###
##                       for mixture normal distribution            ###
#######################################################################
SAVE_M<-function(x,y,q,d,h)
{
  p<-ncol(x)
  n<-nrow(x)
  ys<-slicing(y,h)
  ytilde<-ys$ytilde
  normalmix <- mnEM(x,q)
  mux <- normalmix$mux
  sigx <- normalmix$sigx
  pix <- normalmix$pix
  
  isigx<-sigx
  for (k in 1:q) {isigx[[k]]<-matpower(sigx[[k]],-1)} 
  
  msave<-numeric(0)
  for (j in 1:h)
  {    
    ind<- ((y >= ytilde[j])&(y < ytilde[j+1]))
    nj<-sum(ind)
    xj<-x[ind,]    
    pixj<-pix[ind,]  
    pij<-apply(pixj,2,mean)
    for (k in 1:q)
    {
      xjk<-t(t(xj)-mux[[k]])
      e2xbx<-t(xjk)%*%diag(pixj[,k])%*%xjk/nj
      e1xbx<-t(xjk)%*%pixj[,k]/nj
      vxbx<- e2xbx - e1xbx%*%t(e1xbx)
      msave<-cbind(msave, (pij[k]*diag(p)-isigx[[k]]%*%vxbx) )   
    }
  }  
  b<-svd(msave)$u[,1:d]
  b<-as.matrix(b)
  res<-list()
  res$para<-normalmix
  res$beta<-b
  return(res)
}



#################################################################################
##    SAVE_RM function                                                        ###
##    (Revised SAVE for the predictors with mixture normal distribution)      ###
##    Input: x - n*p matrix                                                   ###
##           y - n-dim vector                                                 ###
##           q - number of clusters                                           ###
##           d - structural dimension                                         ###
##           h - number of slices                                             ###
##           beta - the initial value (we take SIR_M as the initial value)    ###
##           nmlist - the parameters of mixture normal distribution           ###
##   Output: beta - estimated basis function for the central subspace         ###
#################################################################################
SAVE_RM<-function(x,y,q,d,h,beta,nmlist)
{
  n<-nrow(x)
  p<-ncol(x)
  d<-ncol(beta)
  sy<-slicing(y,h)
  ytilde<-sy$ytilde
  
  mux <- nmlist$mux
  sigx <- nmlist$sigx
  pix <- nmlist$pix
  par <- parabx(d,q,nmlist,beta)
  mubx <- par$mubx
  sigbx <- par$sigbx
  prx <- par$pibx
  
  dis<-1    #### The distance between the estimators obtained by the two iterations, 
            #### Set the initial value of dis as 1, and the subsequent value is given by the following iteration process. 
  e<-1e-2   #### A threshold. If dis<e, stop the iteration.
  m<-1      #### The number of iterations.
  betak<-beta
  while(dis>e)
  {
    Asi <- matrix(0,p,p)
    Dsi <- matrix(0,p,p)
    Vsi <- matrix(0,p*p,p*d)
    for (s in 1:h)
    {
      ind <- (y>=ytilde[s])&(y<ytilde[s+1])
      xj <- x[ind,]
      nj<-sum(ind)
      prxj<-prx[ind,]
      for (i in 1:q)
      {
        ### ai=E(pi^2*pj*Ys), bi=E(pi^2*mu2*Ys), ci= E(X^2*pi^2*Ys) ###
        ai <- matrix(0,p,p)
        bi <- matrix(0,p,p)
        ci <- matrix(0,p,p)
        for (j in 1:q)
        { ai <- ai+mean(prx[ind,i]^2*prx[ind,j])*sigx[[j]] }
        for (j in 1:nj)
        {
          ci <- ci + prxj[j,i]^2*xj[j,]%*%t(xj[j,])/nj
          di<-matrix(0,p,p)
          for (k in 1:q)
          {
            pro<-c(t(proj(betak,sigx[[k]]))%*%(x[j,]-mux[[k]])+mux[[k]])
            di<- di + prx[j,k]*pro%*%t(pro)
          }
          bi<- bi + prx[j,i]^2*di/nj 
        }
  
        
        ###  A, B, C, D and V  ###
        Asi <- Asi -ai - bi + ci
        for (j in 1:q)
        { 
          Dsi <- Dsi + mean(prx[ind,i]^2*prx[ind,j])*sigx[[j]]%*%proj(betak,sigx[[j]])
          Bsi <- mean(prx[ind,i]^2*prx[ind,j])*sigx[[j]]
          Csi <- matpower(t(betak)%*%sigx[[j]]%*%betak,-1)%*%t(betak)%*%sigx[[j]]
          Vsi <- Vsi - kronecker(Asi%*%t(Csi),Bsi)
        }
        Usi <- Asi%*%t(Asi) + Asi%*%t(Dsi) + Dsi%*%t(Dsi)
      }
    }
    #####  u   #######
    u <- as.vector(Usi)
    v <- Vsi
   
    #beta(k+1)
    betak1<-matpower(t(v)%*%v,-1)%*%t(v)%*%u
    betak1<-matrix(betak1,p,d)
    betak1 <- betak1/sqrt(sum(betak1^2))
    dis <- dv(betak,betak1)
    betak<-betak1
    m <- m+1
    if(m>1)
      break
  }
  res<-list()
  res$beta<-betak
  return(res)
}























