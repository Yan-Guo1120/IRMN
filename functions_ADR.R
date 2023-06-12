######################################################
#  Some functions for AIME      
#  Aggregate inverse mean estimation for SDR
#  Wang and Yin (2019)
#  submitted to Technometrics in Oct 2019
#  revision submitted in Mar 2020
######################################################

######################################################
####  SIR function    
##    Sliced inverse regression for dimension reduction
##       Li, K. C. 1991, JASA
##  tuning parameter: H - number of slices
######################################################
sir_func = function(X, Y, H){
  n=nrow(X)
  p=ncol(X)
  
  xbar<-colMeans(X)
  
  od<-order(Y)
  xnew<-X[od, ]
  
  msir<-matrix(0, p, p)
  div<-round(seq(from=1,to=n,length=H+1)[c(-1,-(H+1))])
  
  xnew1 <- colMeans(xnew[1:(div[1]-1), ])-xbar
  p1 <- (div[1]-1)/n
  msir <- msir + p1*xnew1%*%t(xnew1)
  
  for (sh in 1:(H-2))  {    
    xnewh<-xnew[div[sh]:(div[sh+1]-1), ]
    ph <- (div[sh+1]-div[sh])/n
    xbarh<-colMeans(xnewh)-xbar
    msir<-msir+ph*xbarh%*%t(xbarh)
  }
  
  xnewH<-xnew[div[H-1]:n, ]
  xbarH<-colMeans(xnewH)-xbar
  pH<-(n-div[H-1]+1)/n
  msir <- msir + pH*xbarH%*%t(xbarH)
  return(msir)
}

#########################################################
## CUME function
## Dimension reduction in regressions through 
##             cumulative slicing estimation
## Zhu, Zhu and Feng, JASA, 2010
#########################################################
knncume <- function(x, y, d) {
  od<-order(y)
  xx<-x[od, ]
  yy<-y[od]
  
  xbar <- colMeans(xx)
  
  kk<-length(y)
  ndiv<-floor(kk/1)
  div<-round(seq(from=1,to=kk,length=ndiv)[c(-1,-ndiv)])
  
  bb_cume<-matrix(0, p, p)
  for (ii in 1:length(div)){
    xbar_ii <- colMeans(xx[1:div[ii], ]) - xbar
    wt <- div[ii]/n
    bb_cume <- bb_cume + (wt)^2*xbar_ii %*% t(xbar_ii)
    # weight function to account for the imbalance of 2 slices
    # several options provided in Zhu et al (2010)
    # wt <- min(div[ii]/(n-div[ii]), (n-div[ii])/div[ii])      
    # bb_cume <- bb_cume + (wt)^2*xbar_ii %*% t(xbar_ii)
  }
  
  svd_cume <- svd(bb_cume/kk)
  bhat_cume <- svd_cume$u[ ,1:d]
  
  dsum <- sum(svd_cume$d)
  dmax <- svd_cume$d[1:d]/dsum
  
  deig <- svd_cume$d
  dratio <- deig[1:d]/deig[2:(d+1)]
  
  if (d==1) {
    bhat_d_cume <- c(bhat_cume, dmax)}
  else {
    bhat_d_cume <- rbind(bhat_cume, dmax, dratio)}
  return(bhat_d_cume)
}


#####################################################################
## Fused Sliced Inverse Regression function                       ###
##  Fused estimators of the central subspace in SDR
##  Cook and Zhang, 2014, JASA
##  tuning parameter: nh -- number of difference H's
#####################################################################
fsir <- function(x, y, nh, d){
  p <- dim(x)[2]
  
  sigmainv <- mppower(cov(x), -1, 1e-6)
  mhat <- matrix(0, p, p)
  
  for (h in 5:(nh+5)){
    msir_h <- sir_func(x, y, h)
    mhat <- mhat + msir_h
  }
  bhat_fsir <- svd(sigmainv%*%mhat/(nh+1))$u[ ,1:d]
  return(bhat_fsir)
}



###################################################
#  SUBROUTINE: Moore-Penrose type power           #
#  Taking power ignoring 0 eigenvalues;           #
#    ignoring criterion=ignore                    #
###################################################
mppower = function(matrix,power,ignore){
  eig = eigen(matrix)
  eval = eig$values
  evec = eig$vectors
  m = length(eval[abs(eval)>ignore])
  tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%
    t(evec[,1:m])
  return(tmp)
}

##########################################################
#  SUBROUTINE: knnsir                                    #
#                                                        #
#    Aggregate Sufficient Dimension Reduction
#    Wang, Yin, Li and Tang (2019)
#    Statistica Sinica
# 
#  tuning parameters: k -- size of kNN
#                     delta -- threshold for # of kNN
#                             (set as median)
#                     l -- the order in Sigma Envelope
###########################################################
knnsir = function(X, Y, d, k) {
  n=nrow(X)
  p=ncol(X)
  
  bb<-matrix(0, p, n)
  wt<-array(0, n)
  ncov<-array(0, dim=c(n, p, p))
  
  for (j in 1:n) {
    
    x<-rep(X[j, ], n)
    dim(x)<-c(p, n)
    xd<-t(X)-x
    dist<-sqrt(diag(t(xd)%*%xd))
    
    od<-order(dist)
    xx<-X[od, ]
    yy<-Y[od]
    dist1<-dist[od]
    
    xnew<-xx[2:(k+1), ]
    ynew<-yy[2:(k+1)]
    sigma1<-cov(xnew)
    ncov[j, , ]<-sigma1
    
    xynew<-cbind(ynew, xnew)
    oy<-order(ynew)
    xnew1<-xnew[oy, ]
    ynew1<-ynew[oy]
    xynew1<-cbind(ynew1, xnew1)
    
    xx1<-xynew1[1:floor(k/2), -1]
    xx2<-xynew1[(floor(k/2)+1):k, -1]
    xbar<-colMeans(xynew1[ ,-1])
    mu1<-colMeans(xx1)
    
    ksi<-mu1-xbar
    bb[ ,j]<-ksi
    wt[j]<-t(ksi)%*%ksi
  }
  
  bb_new<-matrix(0, p, p)
  e1_50 <- quantile(wt, p=0.5)
  
  count<-0
  
  for (j in 1:n) {
    if (wt[j]>e1_50)  {
      count<-count+1
      b1<-bb[ ,j]
      sigma1<-ncov[j, , ]
      hwt<-1
      
      rq<-cbind(b1, sigma1%*%b1, sigma1%*%sigma1%*%b1, sigma1%*%sigma1
                %*%sigma1%*%b1, sigma1%*%sigma1%*%sigma1%*%sigma1%*%b1, sigma1%*%sigma1
                %*%sigma1%*%sigma1%*%sigma1%*%b1)
      kernel<-t(rq)%*%sigma1%*%rq
      pb1<-rq%*%ginv(kernel)%*%t(rq)%*%b1
      pb1<-pb1/sqrt(t(pb1)%*%pb1)[1,1]
      bb_new<-bb_new+hwt[1]*pb1%*%t(pb1)
    }
  }
  
  bksir<-svd(bb_new/count)$u[ ,1:d]
  return(bksir)
}


############################################################
#  SUBROUTINE: adaptive knnsir                             #
#                                                          #
#    Aggregate Sufficient Dimension Reduction
#    Wang, Yin, Li and Tang (2019)
#    Statistica Sinica
#
#  tuning parameters: k -- size of kNN
#                     delta -- threshold for # of kNN
#                             (set as median)
#                     l -- the order in Sigma Envelope
#                     kappa -- softening parameter
#############################################################
rksir<-function(X, Y, k, b0)  {
  
  n<-dim(X)[1]
  p<-dim(X)[2]
  if (is.vector(b0)==TRUE) {d=1}
  else {d<-dim(b0)[2]}
  
  stop<-0
  iter<-0
  while (stop==0) {
    
    iter<-iter+1
    xdt<-X%*%b0
    
    bb<-matrix(0, p, n)
    wt<-array(0, n)
    ncov<-array(0, dim=c(n, p, p))
    
    for (j in 1:n) {
      
      if (d==1) { 
        x<-rep(xdt[j], n)
        xd<-xdt-x
        dist<-abs(xd)
      } else {
        x<-rep(xdt[j, ], n)
        dim(x)<-c(d,n)
        xd<-t(xdt)-x
        dist<-sqrt(diag(t(xd)%*%xd))}

      ## adaptive distance  ##
      x0<-rep(X[j, ], n)
      dim(x0)<-c(p, n)
      xd0<-t(X)-x0
      dist0<-sqrt(diag(t(xd0)%*%xd0))
      dist<-dist+dist0/(iter)^(1/1.5)
      
      od<-order(dist)
      xx<-X[od, ]
      yy<-Y[od]
      dist1<-dist[od]
      
      xnew<-xx[2:(k+1), ]
      ynew<-yy[2:(k+1)]
      sigma1<-cov(xnew)
      ncov[j, , ]<-sigma1
      
      xynew<-cbind(ynew, xnew)
      oy<-order(ynew)
      xnew1<-xnew[oy, ]
      ynew1<-ynew[oy]
      xynew1<-cbind(ynew1, xnew1)
      
      xx1<-xynew1[1:floor(k/2), -1]
      xx2<-xynew1[(floor(k/2)+1):k, -1]
      xbar<-colMeans(xynew1[ ,-1])
      mu1<-colMeans(xx1)
      
      ksi<-mu1-xbar
      bb[ ,j]<-ksi
      wt[j]<-t(ksi)%*%ksi
    }
    
    bb_new<-matrix(0, p, p)
    e1_25 <- quantile(wt, p=0.25)
    e1_50 <- quantile(wt, p=0.50)
    e1_75 <- quantile(wt, p=0.75)
    count<-0
    
    for (j in 1:n) {
      if (wt[j]>e1_50)  {
        count<-count+1
        b1<-bb[ ,j]
        sigma1<-ncov[j, , ]
        hwt<-1 

        rq<-cbind(b1, sigma1%*%b1, sigma1%*%sigma1%*%b1, sigma1%*%sigma1
                  %*%sigma1%*%b1, sigma1%*%sigma1%*%sigma1%*%sigma1%*%b1, sigma1%*%sigma1
                  %*%sigma1%*%sigma1%*%sigma1%*%b1)
        kernel<-t(rq)%*%sigma1%*%rq
        pb1<-rq%*%ginv(kernel)%*%t(rq)%*%b1
        pb1<-pb1/sqrt(t(pb1)%*%pb1)[1,1]
        bb_new<-bb_new+hwt[1]*pb1%*%t(pb1)
      }
    }
    
    bksir<-svd(bb_new/count)$u[ ,1:d]
    if ( max(eigen(b0%*%t(b0)-bksir%*%t(bksir))$values) < 0.01 | iter>10 ) stop=1
    b0<-bksir
  }
  return(bksir)
}



#####################################################################################
#    The main function for Aggregate invese mean estimation   
#                                   Wang and Yin (2019)
#    k - size of nearest neighborhood
#    b0 - initial estimate for B, a p*d matrix
#
#   tuning parameters: k -- size of kNN
#                     delta -- threshold for # of kNN
#                             (set as median)
#                     l -- the order in Sigma Envelope
#                     kappa -- softening parameter
#                     w_1, e_2 -- weight function in cumulative mean estimation
#   output -- estimated direction \hat \B for CS
#####################################################################################

rkcume_d<-function(X, Y, k, b0)  {
  
  n<-dim(X)[1]
  p<-dim(X)[2]
  if (is.vector(b0)==TRUE) {d=1}
  else {d<-dim(b0)[2]}
  
  stop<-0
  iter<-0
  while (stop==0) {
    
    iter<-iter+1
    xdt<-X%*%b0
    
    bb<-matrix(0, p, d*n)
    wt<-matrix(0, n, d)
    wtratio<-matrix(0, n, d)
    ncov<-array(0, dim=c(n, p, p))
    
    for (j in 1:n) {
      
      if (d==1) { 
        x<-rep(xdt[j], n)
        xd<-xdt-x
        dist<-abs(xd)
      } else {
        x<-rep(xdt[j, ], n)
        dim(x)<-c(d,n)
        xd<-t(xdt)-x
        dist<-sqrt(diag(t(xd)%*%xd))}

      ## adaptive distance  ##
      x0<-rep(X[j, ], n)
      dim(x0)<-c(p, n)
      xd0<-t(X)-x0
      dist0<-sqrt(diag(t(xd0)%*%xd0))

      dist<-dist+dist0/(iter)^(1/1.5)
      
      od<-order(dist)
      xx<-X[od, ]
      yy<-Y[od]
      dist1<-dist[od]
      
      xnew<-xx[2:(k+1), ]
      ynew<-yy[2:(k+1)]
      sigma1<-cov(xnew)
      ncov[j, , ]<-sigma1
      
      xynew<-cbind(ynew, xnew)
      oy<-order(ynew)
      xnew1<-xnew[oy, ]
      ynew1<-ynew[oy]
      xynew1<-cbind(ynew1, xnew1)
      
      cume_out<-knncume(xnew1, ynew1, d)
      ksi<-cume_out[1:p,]
      dmax<-cume_out[p+1,]
      dratio<-cume_out[p+2,]
      
      bb[ ,(j*d-d+1):(j*d)]<-ksi
      wt[j, ]<-dmax
      wtratio[j, ]<-dratio
    }
    
    bb_new<-matrix(0, p, p)
    e1_50 <- quantile(wt[,1], p=0.50)

    count<-0
    
    for (j in 1:n) {
      if (wt[j,1]>e1_50)  {
        count<-count+1
        b1<-bb[ ,(j*d-d+1)]
        sigma1<-ncov[j, , ]
        hwt<-1
        
        rq<-cbind(b1, sigma1%*%b1, sigma1%*%sigma1%*%b1, sigma1%*%sigma1
                  %*%sigma1%*%b1, sigma1%*%sigma1%*%sigma1%*%sigma1%*%b1, sigma1%*%sigma1
                  %*%sigma1%*%sigma1%*%sigma1%*%b1)
        kernel<-t(rq)%*%sigma1%*%rq
        pb1<-rq%*%ginv(kernel)%*%t(rq)%*%b1
        pb1<-pb1/sqrt(t(pb1)%*%pb1)[1,1]
        
        bb_new<-bb_new+hwt[1]*pb1%*%t(pb1)
      }
      
      if (wtratio[j,2]>wtratio[j,1] & wt[j,2]>wt[j,1]/3)  {
        count<-count+1
        b1<-bb[ ,(j*d-d+2)]
        sigma1<-ncov[j, , ]
        hwt<-1
        
        rq<-cbind(b1, sigma1%*%b1, sigma1%*%sigma1%*%b1, sigma1%*%sigma1
                  %*%sigma1%*%b1, sigma1%*%sigma1%*%sigma1%*%sigma1%*%b1, sigma1%*%sigma1
                  %*%sigma1%*%sigma1%*%sigma1%*%b1)
        kernel<-t(rq)%*%sigma1%*%rq
        pb1<-rq%*%ginv(kernel)%*%t(rq)%*%b1
        pb1<-pb1/sqrt(t(pb1)%*%pb1)[1,1]
        
        bb_new<-bb_new+hwt[1]*pb1%*%t(pb1)
      }
    }
    
    bkcume<-svd(bb_new/count)$u[ ,1:d]
    if ( max(eigen(b0%*%t(b0)-bkcume%*%t(bkcume))$values) < 0.01 | iter>10 ) stop=1
    b0<-bkcume
  }
  return(bkcume)
}


