% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM2.R
\name{EM2}
\alias{EM2}
\title{Two-dimensional mixed Gaussian density estimation}
\usage{
EM2(data, mu, Sig, tol)
}
\arguments{
\item{data}{The data as a matrix or data frame with two Two columns.}

\item{mu}{The initial mean value of the mixed Gaussian model with two lines.}

\item{Sig}{The initial standard deviation value of the mixed Gaussian model with two lines.}

\item{tol}{The tolerable deviation of mean, When it finally converges.}
}
\value{
a list of mean, standard error,weight of every Two-dimensional Gaussian distribution and times of iterations.
}
\description{
Estimate the two-dimensional joint density with gaussian mixture model by EM algorithm.
}
\examples{
\dontrun{
library(mvtnorm)
data(faithful)
mu <- matrix(c(2, 50, 4, 80),2, 2, byrow = T) 
Sig <- matrix(rep(c(1, 0, 0, 1), 2), ncol = 2, byrow = T)
EM2(faithful,mu, Sig,tol=1e-3)
}
}
