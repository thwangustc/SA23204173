% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/l0tfABESSR.R
\name{DiffMatrix}
\alias{DiffMatrix}
\title{Create a Difference Matrix using the ABESS approach}
\usage{
DiffMatrix(n, q)
}
\arguments{
\item{n}{The number of rows and columns of the matrix.}

\item{q}{The order of the difference (0 or 1).}
}
\value{
A matrix of order \code{n} representing the difference of order \code{q}.
}
\description{
This function creates a difference matrix of a specified size and order based on ABESS idea.
}
\examples{
\dontrun{
DiffMatrix(4, 0)
DiffMatrix(10, 1)
}
}
