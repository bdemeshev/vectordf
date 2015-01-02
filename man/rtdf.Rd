% Generated by roxygen2 (4.0.2): do not edit by hand
\name{rtdf}
\alias{rtdf}
\title{Generate t distribution with vector degrees of freedom}
\usage{
rtdf(n = 1, M = matrix(0), B = diag(ncol(M)), a = rep(1, nrow(M)),
  A = diag(nrow(M)))
}
\arguments{
\item{n}{number of draws from matrix t distribution with vector degrees of freedom}

\item{M}{matrix of expected values (r x s)}

\item{a}{vector of degrees of freedom (r x 1)}

\item{A}{scale matrix (r x r)}
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
X <- rtdf(n = 5, M = matrix(0, nrow=3, ncol=2))
X
}
