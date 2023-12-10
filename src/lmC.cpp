#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace arma;
using namespace Rcpp;

arma::mat subMatrixCpp(const arma::mat& X, const arma::uvec& indices) {
  return X.cols(indices);
}

//' @title Solving Linear Equations using Rcpp and Armadillo
//' @description Solves a linear system of equations by extracting certain columns from matrix X, then applying transformations and solving Ax = b.
//' @param X A matrix representing the coefficient matrix.
//' @param Alnew A vector of unsigned integers indicating the indices of columns of X to be used.
//' @param y A vector representing the right-hand side of the equation.
//' @return A vector that is the solution to the linear system.
//' @examples
//' # Assuming X is a matrix, Alnew is a vector of indices, and y is a vector
//' # solveEquation(X, Alnew, y) solves the system of equations and returns the solution
//' @export
// [[Rcpp::export]]
 arma::vec solveEquation(const arma::mat& X, const arma::uvec& Alnew, const arma::vec& y) {
   arma::mat subMat = subMatrixCpp(X, Alnew);
   arma::mat transSubMat = trans(subMat);
   arma::mat product = transSubMat * subMat;
   arma::vec rightSide = transSubMat * y;
   arma::vec solution = solve(product, rightSide);
   return solution;
 }

