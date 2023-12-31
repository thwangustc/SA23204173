% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/l0tfABESSR.R
\name{l0tfABESSR}
\alias{l0tfABESSR}
\title{Sequential L0 Trend Filtering based on ABESS with R}
\usage{
l0tfABESSR(y = y, kmax = kmax, q = q)
}
\arguments{
\item{y}{A numeric vector of response variables.}

\item{kmax}{The maximum number of variables to include in the model.}

\item{q}{The order of the difference matrix used in the trend filtering.}
}
\value{
A list containing the following elements:
\itemize{
  \item{kopt}{The optimal number of variables chosen by the algorithm.}
  \item{beta}{Coefficients of the model for each iteration up to \code{kmax}.}
  \item{A.all}{A list of indices of variables included in the model for each iteration up to \code{kmax}.}
  \item{A}{The final set of indices of variables included in the model.}
  \item{yhat}{Predicted values using the final model.}
}
}
\description{
This function implements a sequential L0 trend filtering algorithm for sparse signal recovery. It iteratively selects the most significant variables based on a given criterion, such as the Bayesian Information Criterion (BIC).
}
\examples{
\dontrun{
y <- rnorm(100)
result <- l0tfABESSR(y = y, kmax = 10, q = 1)
plot(y)
lines(result$yhat, col = "blue")
}
}
