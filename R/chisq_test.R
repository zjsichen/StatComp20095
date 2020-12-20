#' @title Make a faster version of chisq.test()
#' @name chisq_test
#' @description The chi-square statistic of the cross-contingency table is calculated to verify the independence between two variables.
#' @param x A contingency table matrix.
#' @return Chi-square test gives chi-square values, degrees of freedom and P values.
#' @examples
#' \dontrun{
#' M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
#' dimnames(M) <- list(gender = c("F", "M"),party = c("Democrat","Independent", "Republican"))
#' (M)
#' (Xsq <- chisq.test(M))
#'   chisq_test(M)
#' print(microbenchmark::microbenchmark(chisq_test(M),chisq.test(M)))
#' }
#' @import microbenchmark
#' @importFrom stats dchisq
#' @export
chisq_test <- function(x){
  m <- nrow(x);  n <- ncol(x); N <- sum(x)
  E <- matrix(0,m,n)
  rowsums <- unlist(lapply(1:m, function(i) sum(x[i,]))) # which is used to computed pi.
  colsums <- unlist(lapply(1:n, function(j) sum(x[,j])))
  for (i in 1:m){
    for (j in 1:n) {
      E[i,j] <- rowsums[i]*colsums[j]/N#Compute the number if the calculation is independent.
    }
  }
  # degree of freedom
  df <- (m-1) * (n-1)
  #Calculate the chi-square statistic
  chi_sqr <- sum((x-E)^2/E) 
  p_value <- dchisq(chi_sqr, df = df)
  (test <- list(chi_sqr = chi_sqr, df = df, p_value = p_value))
}



