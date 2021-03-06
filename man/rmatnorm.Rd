% Generated by roxygen2 (4.0.2): do not edit by hand
\name{rmatnorm}
\alias{rmatnorm}
\title{Generate Matrix Normal distribution}
\usage{
rmatnorm(n = 1, M = matrix(0), U = diag(nrow(M)), V = diag(ncol(M)))
}
\arguments{
\item{n}{number of draws from matrix normal distribution}

\item{M}{matrix of expected values (r x s)}

\item{U}{among-row scale covariance matrix (r x r)}

\item{V}{among-column scale covariance matrix (s x s)}
}
\value{
array of generated matrix normal (n x r x s)
}
\description{
Generate Matrix Normal distribution
}
\details{
Generate Matrix Normal distribution using Wikipedia notations. Default is just scalar N(0,1).
}
\examples{
X <- rmatnorm(n = 5, M = matrix(0, nrow=3, ncol=2))
X
}

