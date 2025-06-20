% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ckmultacdoe.R
\name{ckmultacdoe}
\alias{ckmultacdoe}
\title{To identify Outliers in Multi-Response Experiments with Auto-Correlated Errors(1)}
\usage{
ckmultacdoe(trt, Rep, resps, noutli)
}
\arguments{
\item{trt}{Numeric or complex vector containing the treatment levels}

\item{Rep}{Numeric or complex vector containing the Replication}

\item{resps}{data frame containing response variables}

\item{noutli}{a  number indicating the number of outliers}
}
\value{
The output contains Multivariate Normality test and BoxM test for
Residuals of the of the model and Cook's distance for each treatment or a
combination of treatments in Multi - Response Experiments with Auto-Correlated Errors (1)
}
\description{
To identify Outliers in Multi-Response Experiments with Auto-Correlated Errors(1)
}
\examples{
# example code
data(ex2)
ckmultacdoe(ex3$trt,ex3$rep,ex3[,3:4],2)

}
