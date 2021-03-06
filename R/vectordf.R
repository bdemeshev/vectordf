#' vectordf
#'
#' @name vectordf
#' @docType package
#' @author Boris Demeshev 
#' @import dplyr mvtnorm flexsurv
NULL

#' Generate Discrete Markov Chain
#'
#' Generate Discrete Markov Chain
#'
#' Generate Discrete Markov Chain given transition matrix
#'
#' @param trans square transition matrix
#' @param n number of states to generate
#' @param start either a starting state (default is state 1) or a vector of probabilities for starting state
#' @export
#' @return integer vector of simulated states
#' @examples
#' P <- 0.7*diag(3) + 0.1
#' mc <- mc_sim(P, n=100)
#' mc
mc_sim <- function(trans, n=1, start=1) {
  n_states <- ncol(trans) # determine the number of states
  
  sim <- rep(1, n) # reserve space for state sequence
  if (length(start)>1) {
    sim[1] <- sample(x = n_states, size = 1, prob = start) # starting state probabilities
  } else {
    sim[1] <- start # number of starting state
  }
  
  if (n>1) # do simulations
    for (i in 2:n) 
      sim[i] <- sample(x = n_states, size = 1, prob = trans[sim[i-1],])
  
  return(sim)
}


#' Generate Matrix Normal distribution
#'
#' Generate Matrix Normal distribution
#'
#' Generate Matrix Normal distribution using Wikipedia notations. Default is just scalar N(0,1).
#'
#' @param n number of draws from matrix normal distribution
#' @param M matrix of expected values (r x s)
#' @param U among-row scale covariance matrix (r x r)
#' @param V among-column scale covariance matrix (s x s)
#' @export
#' @return array of generated matrix normal (n x r x s)
#' @examples
#' X <- rmatnorm(n = 5, M = matrix(0, nrow=3, ncol=2))
#' X
rmatnorm <- function(n=1, M = matrix(0), U = diag(nrow(M)), V = diag(ncol(M))){
  # M, U, V should be matrices!!! (if scalar, then 1x1)
  mn <- rmvnorm(n = n, mean = as.vector(M), sigma = kronecker(V,U))
  return(array(mn, dim=c(n, nrow(M), ncol(M))))
}




#' Generate gamma distribution with vector degrees of freedom
#'
#' Generate gamma distribution with vector degrees of freedom
#'
#' Generate gamma distribution with vector degrees of freedom. The distribution
#' is described in Shvedov 
#' Algorithm by Alexey Balaev
#'
#' @param n number of draws from matrix gamma distribution with vector degrees of freedom
#' @param a vector of degrees of freedom (r x 1)
#' @param A scale matrix (r x r)
#' @export
#' @return array of generated matrix gamma (n x r x r)
#' @examples
#' X <- rgammadf(n = 5, a = c(1, 2, 3)) 
#' X
rgammadf <- function(n=1, a = rep(1,2), A = diag(length(a))){
  r <- length(a)

  ans <- array(0, dim = c(n, r, r))
  
  # some precalculations to avoid them in cycle
  P <- t(chol(A)) # lower triangular matrix
  Pinv <- solve(P)
  tPinv <- t(Pinv)

  # Balaev, phd thesis, page 116  
  
  L0 <- matrix(0, nrow=r, ncol=r)
  d <- 2*a[r:1] - r + 1:r # gen-gamma pars: a = 1, p = 2, d=(d_1, d_2, ..., d_r)
  
  
  for (i in 1:n) {
    L0[lower.tri(L0)] <- rnorm(r*(r-1)/2, mean=0, sd=1/sqrt(2))
    
    # diag(L0) <- flexsurv::rgengamma(n = r, mu = (log(d)-log(2))/2, 
    #                                sigma = 1/sqrt(d*2),
    #                                Q = sqrt(2/d))
    
    # a little bit faster :)
    diag(L0) <- sqrt(rgamma(n = r, shape = d/2, rate = 1))
    
    W0 <- t(L0) %*% L0 # gamma-Bellman (a, I)
    ans[i,,] <- tPinv %*% W0 %*% Pinv # gamma-Bellman (a, A)
  }
  return(ans)
}


#' Generate t distribution with vector degrees of freedom
#'
#' Generate t distribution with vector degrees of freedom
#'
#' Generate t distribution with vector degrees of freedom. The distribution
#' is described in Shvedov 
#' generalisation of Algorithm by Alexey Balaev
#'
#' @param n number of draws from matrix t distribution with vector degrees of freedom
#' @param a vector of degrees of freedom for underlying gamma Bellman distribution (r x 1)
#' @param A scale matrix for underlying gamma Bellman distribution (r x r)
#' @param B among-column scale covariance matrix for underlying standard normal (s x s)
#' @param M matrix of expected values (r x s)
#' @export
#' @return array of generated matrix t (n x r x s)
#' @examples
#' X <- rtdf(n = 5, a = c(1, 2, 3))
#' X
#' tt <- rtdf(n = 10^5, a=1.5)*sqrt(1.5) # one dimensional t with df=3
#' tt0 <- rt(n = 10^5, df=3)
#' quantile(tt0, probs=c(0.001,0.1,0.3,0.7,0.9,0.999))
#' quantile(tt, probs=c(0.001,0.1,0.3,0.7,0.9,0.999))
rtdf <- function(n = 1, M = matrix(0, nrow=length(a), ncol=ncol(B)), B = diag(1), 
                 a = rep(1, 2), A = diag(length(a)) ){
  r <- nrow(M)
  s <- ncol(M)
  ans <- array(0, dim = c(n, r, s))
  P <- t(chol(A)) # lower triangular matrix
  
  L0 <- matrix(0, nrow=r, ncol=r)
  d <- 2*a[r:1] - r + 1:r # gen-gamma pars: a = 1, p = 2, d=(d_1, d_2, ..., d_r)
  
  I_kron_B <- kronecker(diag(r), B)
  sd_sq2 <- 1/sqrt(2) # constant equal to 1/sqrt(2), for faster cycle
  
  # insides of rmvnorm:
  R <- kronecker(diag(r), chol(B))
  
  for (i in 1:n) {
    # Balaev, phd thesis, page 116    
    L0[lower.tri(L0)] <- rnorm(r*(r-1)/2, mean=0, sd=sd_sq2)
    
    
    # diag(L0) <- flexsurv::rgengamma(n = r, mu = (log(d)-log(2))/2, 
    #                                sigma = 1/sqrt(d*2),
    #                                Q = sqrt(2/d))
    
    # a little bit faster :)
    diag(L0) <- sqrt(rgamma(n = r, shape = d/2, rate = 1))
    
    # Shvedov, WP2/2010/01, page 8
    # vecZ <- mvtnorm::rmvnorm(n=1, mean = rep(0, r*s), sigma = I_kron_B)
    vecZ <- rnorm(r*s) %*% R # insides of rmvnorm
    
    Z <- matrix(vecZ, nrow=r)
    
    
    ans[i,,] <- P %*% solve(L0) %*% Z + M
  }
  return(ans)
}


#' Matrix Normal density function
#'
#' Matrix Normal density function
#'
#' Matrix Normal density function 
#'
#' @param X matrix-point, argument for density function
#' @param M matrix of expected values (r x s)
#' @param U among-row scale covariance matrix (r x r)
#' @param V among-column scale covariance matrix (s x s)
#' @export
#' @return scalar, density at the point X
#' @examples
#' d <- dmatnorm(X = matrix(1, nrow=3, ncol=2), M = matrix(0, nrow=3, ncol=2))
#' d
dmatnorm <- function(X, M = matrix(0), U = diag(nrow(M)), V = diag(ncol(M))){
  # M, U, V should be matrices!!! (if scalar, then 1x1)
  r <- nrow(M)
  s <- ncol(M)
  
  X0 <- X - M
  
  nominator <- exp(-0.5*sum(diag(solve(V) %*% t(X0) %*% solve(U) %*% X0 )))
  denominator <- (2*pi)^(r*s/2)*det(V)^(r/2)*det(U)^(s/2)
    
  return(nominator/denominator)
}


#' Density of matrix t distribution with vector degrees of freedom
#'
#' Density of matrix t distribution with vector degrees of freedom
#'
#' Density of matrix t distribution with vector degrees of freedom. The distribution
#' is described in Shvedov 
#' generalisation of Algorithm by Alexey Balaev
#'
#' @param X matrix-point, argument for density function
#' @param a vector of degrees of freedom for underlying gamma Bellman distribution (r x 1)
#' @param A scale matrix for underlying gamma Bellman distribution (r x r)
#' @param B among-column scale covariance matrix for underlying standard normal (s x s)
#' @param M matrix of expected values (r x s)
#' @export
#' @return scalar, value of density function
#' @examples
#' d <- dtdf(X = matrix(1, nrow=3, ncol=1), a = c(1, 2, 3))
#' d
dtdf <- function(X, M = matrix(0, nrow=length(a), ncol=ncol(B)), B = diag(1), 
                 a = rep(1, 2), A = diag(length(a)) ){
  
  P <- t(chol(A)) # lower triangular matrix
  Z <- solve(P) %*% (X - M) # standartized matrix t distribution
  
  
  r <- nrow(M)
  s <- ncol(M)
  
  message("Not yet implemented")
  
  return(FALSE)  
}



#' Density of matrix gamma distribution with vector degrees of freedom
#'
#' Density of matrix gamma distribution with vector degrees of freedom
#'
#' Density of matrix gamma distribution with vector degrees of freedom. The distribution
#' is described in Shvedov 
#' generalisation of Algorithm by Alexey Balaev
#'
#' @param X matrix-point, argument for density function (r x r)
#' @param a vector of degrees of freedom for underlying gamma Bellman distribution (r x 1)
#' @param A scale matrix for underlying gamma Bellman distribution (r x r)
#' @export
#' @return scalar, value of density function
#' @examples
#' d <- dgammadf(X = matrix(1, nrow=2, ncol=2), a=c(1,2))
#' d
dgammadf <- function(X,  
                 a = rep(1, 2), A = diag(length(a)) ){
  
  P <- t(chol(A)) # lower triangular matrix
  
  message("Not yet implemented")
  
  return(FALSE)  
}
