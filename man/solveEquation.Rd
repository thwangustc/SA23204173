% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{solveEquation}
\alias{solveEquation}
\title{Solving Linear Equations using Rcpp and Armadillo}
\usage{
solveEquation(X, Alnew, y)
}
\arguments{
\item{X}{A matrix representing the coefficient matrix.}

\item{Alnew}{A vector of unsigned integers indicating the indices of columns of X to be used.}

\item{y}{A vector representing the right-hand side of the equation.}
}
\value{
A vector that is the solution to the linear system.
}
\description{
Solves a linear system of equations by extracting certain columns from matrix X, then applying transformations and solving Ax = b.
}
\examples{
# Assuming X is a matrix, Alnew is a vector of indices, and y is a vector
# solveEquation(X, Alnew, y) solves the system of equations and returns the solution
}
