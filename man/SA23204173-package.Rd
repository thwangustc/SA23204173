\name{SA23204173-package}
\alias{SA23204173-package}
\docType{package}
\title{Final Project in Statistical Computing: Implementation of L0 Trend Filtering Based on adaptive best subset selection (ABESS) and alternating minimization induced active set (AMIAS) search Methods}
\description{The L0 trend filtering model has advantages in knot detection and trend estimation compared to traditional trend filtering models. This R package addresses the L0 trend filtering problem from two perspectives: one approach transforms it into a regression problem with an L0 penalty term, processed based on the concept of ABESS, and the other approach starts directly from the dual problem, replicating the AMIAS algorithm. Both the ABESS and AMIAS algorithms involve extensive iterative steps and matrix operations. Here, we use Rcpp functions to optimize the matrix computation part of each iteration step and compare the computational efficiency of the above methods.}
\author{
Tianhao Wang
Maintainer: Tianhao Wang <tianhaowang@mail.ustc.edu.cn>
}
\keyword{package}


