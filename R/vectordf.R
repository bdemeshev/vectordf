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
#' @param N number of states to generate
#' @param start either a starting state (default is state 1) or a vector of probabilities for starting state
#' @export
#' @return integer vector of simulated states
#' @examples
#' P <- 0.6*diag(3) + matrix(0.2, nrow=3, ncol=3)
#' mc <- mc_sim(P, N=100)
#' mc
mc_sim <- function(trans, N=1, start=1) {
  n_states <- ncol(trans) # determine the number of states
  
  sim <- rep(1, N) # reserve space for state sequence
  if (length(start)>1) {
    sim[1] <- sample(x = n_states, size = 1, prob = start) # starting state probabilities
  } else {
    sim[1] <- start # number of starting state
  }
  
  if (N>1) # do simulations
    for (i in 2:N) 
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
#' X <- rmatnorm(n = 5, Mu = matrix(0, nrow=3, ncol=2))
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
#' X <- rgammadf(n = 5, Mu = matrix(0, nrow=3, ncol=2))
#' X
rgammadf <- function(n=1, a = rep(1,2), A = diag(length(a))){
  r <- length(a)
  
  # Balaev, phd thesis, page 116  
  L0 <- matrix(0, nrow=r, ncol=r)
  L0[lower.tri(L0)] <- rnorm(r*(r-1)/2, mean=0, sd=1/sqrt(2))
  
  d <- 2*a[r:1] - r + 1:r # a = 1, p = 2
  diag(L0) <- flexsurv::rgengamma(n = r, mu = (log(d)-log(2))/2, 
                                  sigma = 1/sqrt(d*2),
                                  Q = sqrt(2/d))
    
  P <- t(chol(A))
  Pinv <- solve(P)
  
  W0 <- t(L0) %*% L0 # gamma-Bellman (a, I)
  W <- t(Pinv) %*% W0 %*% Pinv # gamma-Bellman (a, A)
  return(W)
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
#' X <- rtdf(n = 5, M = matrix(0, nrow=3, ncol=2))
#' X
rtdf <- function(n = 1, M = matrix(0, nrow=length(a), ncol=ncol(B)), B = diag(1), 
                 a = rep(1, 2), A = diag(length(a)) ){
  r <- nrow(M)
  s <- ncol(M)
  
  # Balaev, phd thesis, page 116
  L0 <- matrix(0, nrow=r, ncol=r)
  L0[lower.tri(L0)] <- rnorm(r*(r-1)/2, mean=0, sd=1/sqrt(2))

  d <- 2*a[r:1] - r + 1:r # a = 1, p = 2
  diag(L0) <- flexsurv::rgengamma(n = r, mu = (log(d)-log(2))/2, 
                                  sigma = 1/sqrt(d*2),
                                  Q = sqrt(2/d))
  
  # Shvedov, WP2/2010/01, page 8
  vecZ <- rmvnorm(n=1, mean = rep(0, r*s), sigma = kronecker(diag(r), B))
  Z <- matrix(vecZ, nrow=r)
  
  P <- t(chol(A))
  U <- P %*% solve(L0) %*% Z + M
  
  return(U)
}


#' Matrix Normal density function
#'
#' Matrix Normal density function
#'
#' Matrix Normal density function 
#'
#' @param X matrix-point, argument for density function
#' @param Mu matrix of expected values (r x s)
#' @param U among-row scale covariance matrix (r x r)
#' @param V among-column scale covariance matrix (s x s)
#' @export
#' @return scalar, density at the point X
#' @examples
#' d <- dmatnorm(X = matrix(1, nrow=3, ncol=2), Mu = matrix(0, nrow=3, ncol=2))
#' d
dmatnorm <- function(X, Mu = matrix(0), U = diag(nrow(Mu)), V = diag(ncol(Mu))){
  # Mu, U, V should be matrices!!! (if scalar, then 1x1)
  r <- nrow(Mu)
  s <- ncol(Mu)
  
  X0 <- X - Mu
  
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
#' @return array of generated matrix t (n x r x s)
#' @examples
#' d <- dtdf(X = matrix(1, nrow=3, ncol=2), Mu = matrix(0, nrow=3, ncol=2))
#' d
dtdf <- function(X, M = matrix(0), B = diag(ncol(M)), 
                 a = rep(1, nrow(M)), A = diag(nrow(M)) ){
  
  P <- t(chol(A))
  Z <- solve(P) %*% (X - M) # standartized matrix t distribution
  
  
  r <- nrow(M)
  s <- ncol(M)
  
  message("Not yet implemented")
  
  return(FALSE)  
}