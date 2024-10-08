% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cplot.R
\name{cplot}
\alias{cplot}
\title{Scatterplots for complex variables}
\usage{
cplot(x, y, which = 1, ...)
}
\arguments{
\item{x}{vector of a complex variable.}

\item{y}{second vector of a complex variable.}

\item{which}{defines, what type of plot to produce. \code{which=1} will produce
six scatterplots, while \code{which=2} will produce a scatterplot of data after
multidimensional scaling (creating projections of complex variables to x and y axes).}

\item{...}{otehr parameters passed to plot method. Works only for \code{which=2}.}
}
\value{
The function produces a plot and does not return any value
}
\description{
Function produces six scatterplots to show relations between the two complex variables x and y.
}
\details{
The plots are positioned to satisfy two rules:
1. When a scatterplot for a c.r.v. is produced, the real part should be in x-axis,
while the imaginary should be in the y-axis.
2. When parts of variables x and y are compared, the part for $x$ should be in x-axis,
while the part for y should be in y-axis, which should the reflect the idea that x could
 be an explanatory variable for y.
}
\examples{

# Generate random complex variables
x <- complex(real=rnorm(100,10,10), imaginary=rnorm(100,10,10))
y <- complex(real=rnorm(100,10,10), imaginary=rnorm(100,10,10))

cplot(x, y)

}
\references{
\itemize{
\item Svetunkov, S. & Svetunkov I. (2022) Complex Autoregressions. In Press.
}
}
\seealso{
\code{\link[complex]{ccor}}
}
\author{
Ivan Svetunkov, \email{ivan@svetunkov.com}
}
\keyword{univar}
