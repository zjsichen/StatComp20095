% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dic.vs.new.R
\name{dic.vs.new}
\alias{dic.vs.new}
\title{Dichotomy versus Newtonian algorithms}
\usage{
dic.vs.new(f, fgrad, xrange, x0, eps, it_max)
}
\arguments{
\item{f}{The equation that needs to be solved.}

\item{fgrad}{The first derivative of the equation to be solved.}

\item{xrange}{Dichotomy requires the interval of root.}

\item{x0}{The initial value of Newton's method root}

\item{eps}{The deviation of the root that can be tolerated when it finally converges.}

\item{it_max}{Maximum iteration number}
}
\value{
A contingency table containing the root and times of iterations obtained by dichotomy and Newton methods.
}
\description{
Compare the iterative speed of solution of dichotomy method and Newton method.
}
\examples{
\dontrun{
f=function(x){
  2*x^3-4*x^2+3*x-6
}
fgrad=function(x){
  6*x^2-8*x+3
}
dic.vs.new(f,fgrad,xrange=c(0,3),x0=1,eps=1e-5,it_max=1000)
}
}
