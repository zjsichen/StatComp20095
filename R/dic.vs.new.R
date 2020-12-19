#' @title Dichotomy versus Newtonian algorithms
#' @name dic.vs.new
#' @description Compare the iterative speed of solution of dichotomy method and Newton method.
#' @param f The equation that needs to be solved.
#' @param fgrad The first derivative of the equation to be solved.
#' @param xrange Dichotomy requires the interval of root.
#' @param x0 The initial value of Newton's method root
#' @param eps The deviation of the root that can be tolerated when it finally converges.
#' @param it_max Maximum iteration number
#' @return A contingency table containing the root and times of iterations obtained by dichotomy and Newton methods.
#' @examples
#' \dontrun{
#' f=function(x){
#'   2*x^3-4*x^2+3*x-6
#' }
#' fgrad=function(x){
#'   6*x^2-8*x+3
#' }
#' dic.vs.new(f,fgrad,xrange=c(0,3),x0=1,eps=1e-5,it_max=1000)
#' }
#' @importFrom mvtnorm dmvnorm
#' @import Rcpp
#' @import microbenchmark
#' @import MASS
#' @import mvtnorm
#' @import nloptr
#' @import xtable
#' @import bootstrap
#' @import boot
#' @import DAAG
#' @import RANN
#' @import energy
#' @import Ball
#' @import devtools
#' @import roxygen2
#' @export
dic.vs.new <- function(f,fgrad,xrange,x0,eps,it_max){
  #dichotomy
  start=xrange[1];end=xrange[2]
  dichotomy<-function(f,xrange,start,end,eps=1e-5,it_max=1000){
    i=0;a=xrange[1];b=xrange[2]
    while(i<=it_max){
      if(f(a)*f(b)>0)
        list(fail="find root is fail!")
      else{
        repeat{
          if(abs(b-a)<eps) break;
          x<-(a+b)/2
          if(f(a)*f(x)<0) b<-x else a<-x
          i=i+1
        }
      }
      return(list(x,i) )#root and times of iterations
    }
  }
  #Newton
  funs=function(f,fgrad,x){
    f=x
    J=f(x)/fgrad(x)
    list(f=f,J=J);
  }
  Newtons=function(fun,f,fgrad,x,eps=1e-5,it_max=1000){
    index=0;k=1
    while(k<=it_max){
      norm=abs(fun(f,fgrad,x)$J)
      x=x-fun(f,fgrad,x)$J
      if(norm<eps){ index=1;break}
      k=k+1
    }
    if(index==0) list(fail="find root is fail!")
    return(c(x,k))#root and times of iterations
  }
  #Compare the two algorithms
  result_dic=dichotomy(f,xrange,eps,it_max)
  result_new=Newtons(funs,f,fgrad,x0,eps,it_max)
  result=cbind(dichotomy=result_dic,Newtons=result_new)
  rownames(result)=c("root","iters")
  return(result)
}
