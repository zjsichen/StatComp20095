#' @title skewness coeff
#' @description Computes the sample skewness coeff. 
#' @param x A contingency table matrix.
#' @return skewness coeff of x
#' @examples
#' \dontrun{
#' x=runif(30)
#' sk(x)
#' }
#' @import stats
#' @export
sk <- function(x) {
  xbar <- mean(x) 
  m3 <- mean((x - xbar)^3) 
  m2 <- mean((x - xbar)^2) 
  return( m3 / m2^1.5 )
}


#' @title multivariate skewness test
#' @description Computes the multivariate skewness coeff. 
#' @param x A contingency table matrix.
#' @return skewness coeff of x
#' @examples
#' \dontrun{
#' d=2
#' sigma=matrix(rep(0,d*d),nrow = d)
#' for(i in 1:d){
#'   for (j in 1:d){
#'       sigma[i,j]=0.5^(abs(i-j))
#'   }
#' }
#' mu=matrix(rep(0,d),nrow = 1)
#' x <- mvrnorm(n[i],mu,sigma) 
#' B_1d(x)
#' }
#' @importFrom  stats cov
#' @export
B_1d<- function(x) {
  #computes the sample skewness coeff. 
  xbar <- apply(x,2,mean)
  n=length(x[,1])
  sigma=cov(x)
  b_1d=matrix(rep(0,n*n),nrow = n)
  b_1d=(as.matrix(x-xbar)%*%solve(sigma)%*%t(as.matrix(x-xbar)))^3
  return(mean(b_1d))
}



#' @title Count Five test
#' @description Count Five test is used to test whether the two sets of data have the same variance
#' @param x A vector.
#' @param y Another vector with the same length of x.
#' @return The value of a test statistic.
#' @examples
#' \dontrun{
#' x <- rnorm(n=50, 0, 1) 
#' y <- rnorm(n=50, 0, 2) 
#' count5test(x,y)
#' }
#' @export
count5test <- function(x, y) {
  X <- x - mean(x) 
  Y <- y - mean(y) 
  outx <- sum(X > max(Y)) + sum(X < min(Y)) 
  outy <- sum(Y > max(X)) + sum(Y < min(X)) # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}


#' @title Count Five test and F test
#' @description Both calculate the F statistic and Count Five test statistic.
#' @param x A vector.
#' @param y Another vector with the same length of x.
#' @return The p value of the F statistic and Count Five test statistic.
#' @examples
#' \dontrun{
#' x <- rnorm(n=50, 0, 1) 
#' y <- rnorm(n=50, 0, 2) 
#' count5.F(x,y)
#' }
#' @importFrom stats var.test
#' @export
count5.F=function(x,y){
  p.reject=numeric(2)
  p.reject[1]=count5test(x,y)
  p.reject[2]=as.integer(var.test(x,y)$p.value<0.055)#F test
  return(p.reject)
}

f_Lap <- function(x, mu,lambda) {
  f=(1/(2*lambda))*exp(-abs(x-mu)/lambda)
}
#' @title Implement a random walk Metropolis sampler with R
#' @description Implement a random walk Metropolis sampler for generating the standard Laplace distribution.
#' @param mu The mean of Laplace distribution.
#' @param lambda The lambda parameter of Laplace distribution.
#' @param sigma The standard deviation of the norm target distribution.
#' @param x0 The initial value of a random walk sequence.
#' @param N The length of a random walk sequence.
#' @return a random walk sequence and the number of times it was accepted.
#' @examples
#' \dontrun{
#' set.seed(1234)
#' mu=0;lambda=1;sigma=1
#' x0=10
#' N=15000
#' rw<- rw.Metropolis(mu,lambda, sigma, x0, N)
#' }
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @export
rw.Metropolis <- function(mu, lambda, sigma,x0, N) {
  x <- numeric(N) 
  x[1] <- x0 
  u <- runif(N) 
  k<-0 
  for (i in 2:N) { 
    y <- rnorm(1, x[i-1], sigma) 
    if (u[i] <= (f_Lap(y, mu,lambda) / f_Lap(x[i-1], mu,lambda))) x[i] <- y else { 
      x[i] <- x[i-1] 
      k<-k+1
    } 
  }
  return(list(x=x, k=k))
}


Tn <- function(z, ix, sizes,k) { 
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2 
  if(is.vector(z)) z <- data.frame(z,0); 
  z <- z[ix, ]; 
  NN <- nn2(data=z, k=k+1) 
  block1 <- NN$nn.idx[1:n1,-1] 
  block2 <- NN$nn.idx[(n1+1):n,-1] 
  i1 <- sum(block1 < n1+.5); 
  i2 <- sum(block2 > n1+.5) 
  (i1 + i2) / (k * n)
}
#' @title evaluating the performance of the NN
#' @description evaluating the performance of the NN used to evaluate whether the variance is homophase
#' @param z Sample data to be tested
#' @param sizes Total sample size
#' @param k Nearest neighbor quantity selection
#' @param R The number of the bootstrap
#' @return KNN statistic and its P value
#' @examples
#' \dontrun{
#' m <- 1e3; k<-3; p<-2; set.seed(12345) 
#' n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#' x <- matrix(rnorm(n1*p),ncol=p); y <- cbind(rnorm(n2),rnorm(n2,0,1.8)); z <- rbind(x,y) 
#' eqdist.nn(z,N,k,R)
#' }
#' @import RANN
#' @import boot
#' @export
eqdist.nn <- function(z,sizes,k,R){ 
  boot.obj <- boot(data=z,statistic=Tn,R=R, sim = "permutation", sizes = sizes,k=k) 
  ts <- c(boot.obj$t0,boot.obj$t) 
  p.value <- mean(ts>=ts[1]) 
  list(statistic=ts[1],p.value=p.value)
}  