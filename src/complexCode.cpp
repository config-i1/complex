#include <RcppArmadillo.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/* # Function calculates inverse of matrix of complex variables */
//' Function calculates inverse of matrix of complex variables
//'
//' The function accepts a square complex matrix and returns inverse of it.
//'
//' @param x The square matrix of complex variables.
//'
//' @template author
//'
//' @return The function returns a matrix of the same size as the original
//' matrix \code{x}
//'
//' @seealso \link[base]{solve}
//'
//' @examples
//'
//' \dontrun{invert(matrix(complex(real=c(1,2), imaginary=c(1.1,2.1)), 2, 2))}
//'
//' @useDynLib complex
//' @export
// [[Rcpp::export]]
arma::cx_mat invert(arma::cx_mat x) {
  return x.i();
}
