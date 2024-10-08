% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{invert}
\alias{invert}
\title{Function calculates inverse of matrix of complex variables}
\usage{
invert(x)
}
\arguments{
\item{x}{The square matrix of complex variables.}
}
\value{
The function returns a matrix of the same size as the original
matrix \code{x}
}
\description{
The function accepts a square complex matrix and returns inverse of it.
}
\examples{

invert(matrix(complex(real=c(1,2), imaginary=c(1.1,2.1)), 2, 2))

}
\seealso{
\link[base]{solve}
}
\author{
Ivan Svetunkov, \email{ivan@svetunkov.com}
}
