#############################################################
##       This code is provided by Tony Cai (2019)          ##
##  Parameter estimation and clustering for a two-class    ##
##       Gaussian mixture via the EM algorithm.            ##
##  For algorithm details, see their paper                 ##
##  "CHIME: CLUSTERING OF HIGH-DIMENSIONAL GAUSSIAN        ##
##   MIXTURES WITH EM ALGORITHM AND ITS OPTIMALITY"        ##
#############################################################



######################################################################################
###   Main function is CHIME, which deals with sparse clustering of                 ##
###   high-dimensional Gaussian mixture models using EM algorithm.                  ##             ###
###   It is based on the EM algorithm, and iteratively estimates the mixing         ##
###   ratio omega, component means mu_1, mu_2 and beta=inv(Sigma)*(mu_1 - mu_2)     ##
###                                                                                 ##
###   Outputs:                                                                      ##                                   ##
###      omega: ratio of two clusters                                               ##
###      mu: mean vectors for different clusters                                    ##
###      Vn: common covariance matrix of both clusters                              ##
###      gamma: posterior probabilities that data comes from the first cluster      ##
###   Inputs:                                                                       ##
###      z: N by p data matrix                                                      ##
###      omega0: initialization for \omega                                          ##
###      mu0: p x 2, initialization of [\mu_1, \mu_2]                               ##
###      beta0: p x 1, initialization of beta,                                      ##
###      rho: a vector used as constant multiplier for the penalty parameter        ##
######################################################################################
CHIME <- function(z, omega0, mu0, beta0, rho)
{
  lambda <- 0.1
  maxIter <- 50
  tol <- 1e-6
  N <- nrow(z)
  p <- ncol(z)
  nrho <- length(rho)
  omega <- matrix(0,nrho,1)
  mu <- list(matrix(0,p,2))
  beta <- matrix(0,p,nrho)
  for (loop_rho in 1:nrho)
  {
    lam_c <- lambda + rho[loop_rho]*sqrt(log(p)/N)
    
    old_omega <- omega0
    old_mu <- mu0
    old_beta <- beta0 
    
    iter <- 1
    diff <- 100
    while ((diff >= tol) & (iter < maxIter))
    {
      ## E-step: calculate gamma
      gamma<-old_omega/((1-old_omega)*exp((z-matrix(1,N,1)%*%t(apply(old_mu,1,mean)))%*%old_beta) + old_omega )
      
      ## M-step: update omega, mu
      new_omega <- mean(gamma)
      tmp1 <- t(apply(diag(c(rep(1,N)-gamma))%*%z,2,mean))/(1-new_omega) 
      tmp2 <- t(apply(diag(c(gamma))%*%z,2,mean))/new_omega
      new_mu <- cbind(t(tmp1),t(tmp2))
      
      ## Update the empirical covariance matrix Vn
      x <- as.vector(sqrt(rep(1,N)-gamma))*t(t(z)-new_mu[,1])
      y <- as.vector(sqrt(gamma))*t(t(z)-new_mu[,2])
      Vn <- 1/N*t(x)%*%x + 1/N*t(y)%*%y
      while (cond(Vn) > 1e+6)
        Vn <- Vn + sqrt(log(p)/N)*diag(p)
      
      ## M-step: update beta
      delta <- tmp1 - tmp2
      beta_init <- solve(Vn)%*%t(delta)
      new_beta <- clime(beta_init, Vn, t(delta), lam_c)
      
      lam_c <- 0.7 * lam_c + rho[loop_rho]*sqrt(log(p)/N)
      
      ## Calculate the difference between the new value and the old value
      diff <- sqrt(sum((new_beta - old_beta)^2)) + sqrt(max(eigen(t(new_mu - old_mu)%*%(new_mu - old_mu))$values)) + abs(new_omega - old_omega)
      
      old_omega <- new_omega
      old_mu <- new_mu
      old_beta <- new_beta       
      iter <- iter + 1
    }
    ## Save the estimate
    omega[loop_rho] <- new_omega
    mu[[loop_rho]] <- new_mu
    beta[,loop_rho] <- new_beta
  }
  res <- list(omega, mu, Vn, gamma)
  return(res)
}



######################################################################
##                clime: L1 dantzig selector                        ##
## Solves min_x  ||x||_1  subject to  ||Ax-b||_\infty <= epsilon    ##
##     Recast as linear program                                     ##
##     min_{x,u}  sum(u)  s.t.  x - u <= 0                          ##
##                         -x - u <= 0                              ##
##            (Ax-b) - epsilon <= 0                                 ##
##            -(Ax-b) - epsilon <= 0                                ##
##     and use primal-dual interior point method.                   ##
##                                                                  ##
## Inputs:                                                          ##
##     x0: px1 vector, initial point.                               ##
##     A : a pxN matrix.                                            ##
##     b : px1 vector of observations.                              ##
##     epsilon: scalar or Nx1 vector of correlation constraints     ##
######################################################################
clime <- function(x0, A,  b, epsilon)
{
  pdtol <- 1e-3    ## Tolerance for primal-dual algorithm (algorithm terminates if the duality gap is less than pdtol).  
  
  pdmaxiter <- 50  ## Maximum number of primal-dual iterations.  
  
  N <- length(x0)
  alpha <- 0.01
  beta <- 0.5
  mu <- 10
  gradf0 <- matrix(c(rep(0,N),rep(1,N)),2*N,1)
  
  if (max(abs((A%*%x0 - b))) >  epsilon )
  {
    cat('Starting point infeasible: using x0 = At*inv(AAt)*y.')
    initrho <- 10*epsilon
    nrowA <- nrow(A)
    initcount <- 0
    while (max(abs((A%*%x0 - b))) >  epsilon )
    {
      x0 <- solve(A+diag(nrowA)*initrho, b)
      hcond <- 1/cond(A+diag(nrowA)*initrho)
      if (hcond < 1e-14)
      {
        cat('A*At is ill-conditioned: cannot find starting point, return initial value')
        xp <- x0
        return (xp)
      }
      initcount <- initcount + 1
      initrho <- initrho/1.2
      if (initcount > 50) 
        break
    }
    
    if ((hcond < 1e-14) | initcount > 50 )
    {
      cat('A*At is ill-conditioned: cannot find starting point, return initial value')
      xp <- x0
      return (xp)
    }
  }
  
  x <- x0
  u <- (0.95)*abs(x0) + (0.10)*max(abs(x0))
  
  ## set up for the first iteration
  Atr <- A%*%x - b
  fu1 <- x - u
  fu2 <- -x - u
  fe1 <- Atr - epsilon
  fe2 <- -Atr - epsilon
  lamu1 <- -1/fu1
  lamu2 <- -1/fu2
  lame1 <- -1/fe1
  lame2 <- -1/fe2
  
  AtAv <- t(A)%*%(lame1-lame2)
  
  ## sdg = surrogate duality gap
  sdg <- -t(rbind(fu1, fu2, fe1, fe2))%*%rbind(lamu1, lamu2, lame1, lame2)
  tau <- as.numeric(mu*(4*N)/sdg)
  
  ## residuals
  rdual <- gradf0 + rbind(lamu1-lamu2 + AtAv, -lamu1-lamu2)
  rcent <- -rbind(lamu1*fu1, lamu2*fu2, lame1*fe1, lame2*fe2) - (1/tau)
  resnorm <- sqrt( sum( (rbind(rdual,rcent))^2 ) )
  
  ## iterations
  pditer <- 0
  while( (sdg >= pdtol) & (pditer < pdmaxiter) )
  {
    ## solve for step direction
    w2 <- - 1 - as.numeric(1/tau)*(1/fu1 + 1/fu2)
    
    sig11 <- -lamu1/fu1 - lamu2/fu2
    sig12 <- lamu1/fu1 - lamu2/fu2
    siga <- -(lame1/fe1 + lame2/fe2)
    sigx <- sig11 - sig12^2/sig11
    
    w1 <- -as.numeric(1/tau)*( t(A)%*%(1/fe2-1/fe1) + 1/fu2 - 1/fu1 )
    w1p <- w1 - (sig12/sig11)*w2
    k <- nrow(siga)
    Hp <- t(A) %*% diag(as.vector(siga),k) %*% A + diag(as.vector(sigx),N)
    
    dx <- solve(Hp, w1p)
    hcond <- 1/cond(Hp)
    if (hcond < 1e-14)
    {
      cat('Matrix ill-conditioned.  Returning previous iterate.  (See Section 4 of notes for more information.)')
      xp <- x
      return (xp)
    }
    AtAdx <- A%*%dx
    du = w2/sig11 - (sig12/sig11)*dx
    
    dlamu1 <- -(lamu1/fu1)*(dx-du) - lamu1 - as.numeric(1/tau)/fu1
    dlamu2 <- -(lamu2/fu2)*(-dx-du) - lamu2 - as.numeric(1/tau)/fu2
    dlame1 <- -(lame1/fe1)*(AtAdx) - lame1 - as.numeric(1/tau)/fe1
    dlame2 <- -(lame2/fe2)*(-AtAdx) - lame2 - as.numeric(1/tau)/fe2
    
    AtAdv = t(A)%*%(dlame1-dlame2)
    
    
    ## find minimal step size that keeps ineq functions < 0, dual vars > 0
    iu1 <- dlamu1 < 0
    iu2 <- dlamu2 < 0
    ie1 <- dlame1 < 0
    ie2 <- dlame2 < 0
    ifu1 <- (dx-du) > 0
    ifu2 <- (-dx-du) > 0
    ife1 <- AtAdx > 0
    ife2 <- -AtAdx > 0
    smax = min(1,min(rbind( as.matrix(-lamu1[iu1]/dlamu1[iu1]), as.matrix(-lamu2[iu2]/dlamu2[iu2]),
                            as.matrix(-lame1[ie1]/dlame1[ie1]), as.matrix(-lame2[ie2]/dlame2[ie2]),
                            as.matrix(-fu1[ifu1]/(dx[ifu1]-du[ifu1])), as.matrix(-fu2[ifu2]/(-dx[ifu2]-du[ifu2])),
                            as.matrix(-fe1[ife1]/AtAdx[ife1]), as.matrix(-fe2[ife2]/(-AtAdx[ife2])) )))
    s <- 0.99*smax
    
    ## backtracking line search
    suffdec <- 0
    backiter <- 0
    while (suffdec==0)
    {
      xp <- x + s*dx 
      up <- u + s*du
      Atrp <- Atr + s*AtAdx 
      AtAvp <- AtAv + s*AtAdv
      fu1p <- fu1 + s*(dx-du)
      fu2p <- fu2 + s*(-dx-du)
      fe1p <- fe1 + s*AtAdx 
      fe2p <- fe2 + s*(-AtAdx)
      lamu1p <- lamu1 + s*dlamu1
      lamu2p <- lamu2 + s*dlamu2
      lame1p <- lame1 + s*dlame1 
      lame2p <- lame2 + s*dlame2
      rdp <- gradf0 + rbind(lamu1p-lamu2p + AtAvp, -lamu1p-lamu2p)
      rcp <- -rbind(lamu1p*fu1p, lamu2p*fu2p, lame1p*fe1p,lame2p*fe2p) - (1/tau)
      suffdec <- ( max(svd(rbind(rdp, rcp))$d) <= (1-alpha*s)*resnorm )
      s <- beta*s
      backiter <- backiter+1
      if (backiter > 32)
      {
        # print('Stuck backtracking, returning last iterate.')
        xp <- x
        return (xp)
      }
    }
    
    ## setup for next iteration
    x <- xp  
    u <- up
    Atr <- Atrp
    AtAv <- AtAvp
    fu1 <- fu1p
    fu2 <- fu2p 
    fe1 <- fe1p 
    fe2 <- fe2p
    lamu1 <- lamu1p 
    lamu2 <- lamu2p 
    lame1 <- lame1p 
    lame2 <- lame2p
    
    sdg <- -t(rbind(fu1, fu2, fe1, fe2))%*%rbind(lamu1, lamu2, lame1, lame2)
    tau <- as.numeric(mu*(4*N)/sdg)
    
    rdual <- rdp
    rcent <- -rbind(lamu1*fu1, lamu2*fu2, lame1*fe1, lame2*fe2) - as.numeric(1/tau)
    resnorm <- sqrt(sum(rbind(rdual, rcent)^2))
    
    pditer <- pditer+1
  }
  return (xp)
}


cond <- function(A)
{
  res <- max(svd(A)$d)/min(svd(A)$d)
  return(res)
}












































