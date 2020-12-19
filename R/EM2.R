#' @title Two-dimensional mixed Gaussian density estimation
#' @name EM2
#' @description Estimate the two-dimensional joint density with gaussian mixture model by EM algorithm.
#' @param data The data as a matrix or data frame with two Two columns.
#' @param mu The initial mean value of the mixed Gaussian model with two lines.
#' @param Sig The initial standard deviation value of the mixed Gaussian model with two lines.
#' @param tol The tolerable deviation of mean, When it finally converges.
#' @return a list of mean, standard error,weight of every Two-dimensional Gaussian distribution and times of iterations.
#' @examples
#' \dontrun{
#' library(mvtnorm)
#' data(faithful)
#' mu <- matrix(c(2, 50, 4, 80),2, 2, byrow = T) 
#' Sig <- matrix(rep(c(1, 0, 0, 1), 2), ncol = 2, byrow = T)
#' EM2(faithful,mu, Sig,tol=1e-3)
#' }
#' @importFrom mvtnorm dmvnorm
#' @importFrom  mixtools mvnormalmixEM
#' @importFrom  Rcpp evalCpp
#' @import knitr
#' @import scales
#' @useDynLib StatComp20095
#' @export
EM2 <- function(data,mu, Sig,tol){
  C=ncol(mu)
  data=as.matrix(data)
  n=nrow(data);p=ncol(data)
  Z <- matrix(c(rep(1, n), rep(0, (n - 1) * C)), n, C)
  mu0=mu
  pi <- rep(1 / C, C)
  dif=1;iter=0
  while (dif >tol) { 
    for (k in 1:n) { 
      for (l in 1:C) { 
        Z[k,l] <- Z[k,l] <- pi[l] * dmvnorm(data[k,], mu[l,], Sig[(p * (l - 1) + 1):(p * l),])
      }
      Z[k,] <- Z[k,] / sum(Z[k,]) 
    } 
    # update Z (E-step) 
    pi <- colMeans(Z) 
    # update pi (M-step-1) 
    for (i in 1:C) { 
      mu[i,] <- t(Z[,i]) %*% data / sum(Z[,i]) 
      # update mu (M-step-2) 
      sumsig <- Z[1,i] * (data[1,] - mu[i,]) %*% t(data[1,] - mu[i,]) 
      for (k in 2:n) { 
        sumsig <- sumsig + Z[k,i] * (data[k,] - mu[i,]) %*% t(data[k,] - mu[i,]) 
      }
      Sig[(p * (i - 1) + 1):(p * i),] <- sumsig / sum(Z[,i]) 
      # update Sigma (M-step-3) 
    }
    dif=max(mu-mu0)
    mu0=mu
    iter=iter+1
  }
  #result of EM
  list("mu"=mu,'Sigma'=Sig,'lambda'=pi,'iters'=iter)
}
