#include <Rcpp.h>
using namespace Rcpp;

//' @title Implement a random walk Metropolis sampler.
//' @description Implement a random walk Metropolis sampler for generating the standard Laplace distribution.
//' @param mu The mean of Laplace distribution.
//' @param lambda The lambda parameter of Laplace distribution.
//' @param sigma The standard deviation of the norm target distribution.
//' @param x0 The initial value of a random walk sequence.
//' @param N The length of a random walk sequence.
//' @return a random walk sequence and the number of times it was accepted.
//' @examples
//' \dontrun{
//' set.seed(1234)
//' mu=0;lambda=1;sigma=1
//' x0=10
//' N=15000
//' rw<- rwMetropolisC(mu,lambda, sigma, x0, N)
//' }
//' @useDynLib StatComp20095
//' @export
// [[Rcpp::export]]
List rwMetropolisC(double mu,double lambda, double sigma,double x0,int N) {
  NumericVector x(N);
  x[0]=x0;
  NumericVector u=runif(N,0,1);
  int k=0;
  for(int j =1;j<N;j++){
    double y=rnorm(1,x[j-1],sigma)[0];
    double fy=(1/(2*lambda))*exp(-abs(y-mu)/lambda);
    double fx=(1/(2*lambda))*exp(-abs(x[j-1]-mu)/lambda);
    double rate=fy/fx;
    if (u[j] <=rate) {
      x[j]=y;
    } else {
      x[j]=x[j-1];
      k++;
    }
  }
  return List::create(
    _["x"] = x,
    _["k"] = k
  );
}
