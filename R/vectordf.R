#' vectordf
#'
#' @name vectordf
#' @docType package
#' @author Boris Demeshev 
#' @import dplyr mvtnorm
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
#' @param Mu matrix of expected values (r x s)
#' @param U among-row scale covariance matrix (r x r)
#' @param V among-column scale covariance matrix (s x s)
#' @export
#' @return array of generated matrix normal (n x r x s)
#' @examples
#' X <- rmatnorm(n = 5, Mu = matrix(0, nrow=3, ncol=2))
#' X
rmatnorm <- function(n=1, Mu = matrix(0), U = diag(nrow(Mu)), V = diag(ncol(Mu))){
  # Mu, U, V should be matrices!!! (if scalar, then 1x1)
  mn <- rmvnorm(n = n, mean = as.vector(Mu), sigma = kronecker(V,U))
  return(array(mn, dim=c(n, nrow(Mu), ncol(Mu))))
}
