% Generated by roxygen2 (4.0.2): do not edit by hand
\name{rtdf}
\alias{rtdf}
\title{Generate t distribution with vector degrees of freedom}
\usage{
rtdf(n = 1, M = matrix(0, nrow = length(a), ncol = ncol(B)), B = diag(1),
  a = rep(1, 2), A = diag(length(a)))
}
\arguments{
\item{n}{number of draws from matrix t distribution with vector degrees of freedom}

\item{M}{matrix of expected values (r x s)}

\item{B}{among-column scale covariance matrix for underlying standard normal (s x s)}

\item{a}{vector of degrees of freedom for underlying gamma Bellman distribution (r x 1)}

\item{A}{scale matrix for underlying gamma Bellman distribution (r x r)}
}
\value{
array of generated matrix t (n x r x s)
}
\description{
Generate t distribution with vector degrees of freedom
}
\details{
Generate t distribution with vector degrees of freedom. The distribution
is described in Shvedov
generalisation of Algorithm by Alexey Balaev
}
\examples{
X <- rtdf(n = 5, a = c(1, 2, 3))
X
tt <- rtdf(n = 10^5, a=1.5)*sqrt(1.5) # one dimensional t with df=3
tt0 <- rt(n = 10^5, df=3)
quantile(tt0, probs=c(0.001,0.1,0.3,0.7,0.9,0.999))
quantile(tt, probs=c(0.001,0.1,0.3,0.7,0.9,0.999))
}

