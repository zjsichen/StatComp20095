% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chisq_test.R
\name{chisq_test}
\alias{chisq_test}
\title{Make a faster version of chisq.test()}
\usage{
chisq_test(x)
}
\arguments{
\item{x}{A contingency table matrix.}
}
\value{
Chi-square test gives chi-square values, degrees of freedom and P values.
}
\description{
The chi-square statistic of the cross-contingency table is calculated to verify the independence between two variables.
}
\examples{
\dontrun{
M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
dimnames(M) <- list(gender = c("F", "M"),party = c("Democrat","Independent", "Republican"))
(M)
(Xsq <- chisq.test(M))
  chisq_test(M)
print(microbenchmark::microbenchmark(chisq_test(M),chisq.test(M)))
}
}
