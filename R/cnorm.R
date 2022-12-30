#' Complex Normal Distribution
#'
#' Density, cumulative distribution, quantile functions and random number
#' generation for the Complex Normal distribution.
#'
#' Complex Normal distribution is a special case of a multivariate normal
#' distribution, which is parametrised using direct and conjugate variances
#' instead of the covariance matrix.
#'
#' These functions are just wrappers for the functions from the \code{mvtnorm}
#' package.
#'
#' Note that \code{sigma2} and \code{varsigma2} are the conjugate and direct
#' variances, not the standard deviations!
#'
#' Both \code{pcnorm} and \code{qcnorm} are returned for the lower
#' tail of the distribution.
#'
#' All the functions are defined for non-negative values only.
#'
#' @template author
#' @keywords distribution
#'
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. Should be a single number.
#' @param mu vector of location parameters (means).
#' @param sigma2 vector of conjugate variances.
#' @param varsigma2 vector of direct variances.
#' @param log if \code{TRUE}, then probabilities are returned in
#' logarithms.
#' @param ... Other parameters passed to the \code{mvtnorm} functions.
#'
#' @return Depending on the function, various things are returned
#' (usually either vector or scalar):
#' \itemize{
#' \item \code{dcnorm} returns the density function values for the
#' provided parameters, based on \link[mvtnorm](Mvnorm) function.
#' \item \code{pcnorm} returns the values of the cumulative function
#' for the provided parameters, based on \link[mvtnorm](pmvnorm) function.
#' \item \code{qcnorm} returns quantiles of the distribution,
#' based on \link[mvtnorm](qmvnorm) function.
#' \item \code{rcnorm} returns a vector of random variables
#' generated from the Complex Normal distribution,
#' based on \link[mvtnorm](Mvnorm) function.
#' }
#'
#' @examples
#' dcnorm(89+90i, 100+100i, 1, 1+1i)
#' pcnorm(89+90i, 100+100i, 1, 1+1i)
#' qcnorm(c(0.025,0.975), 100+100i, 1, 1+1i)
#' rcnorm(1000, 100+100i, 1, 1+1i)
#'
#' @rdname cnormal
#' @aliases cnormal dcnorm
#' @importFrom mvtnorm dmvnorm pmvnorm qmvnorm rmvnorm
#' @export dcnorm
dcnorm <- function(q, mu=0+0i, sigma=1, varsigma=0+0i, log=FALSE, ...){

    # Create a covariance matrix based on the provided variances
    Sigma <- matrix(c((sigma+Re(varsigma)),Im(varsigma),Im(varsigma),(sigma-Re(varsigma)))/2, 2, 2);

    cnormReturn <- dmvnorm(complex2vec(q), mean=complex2vec(mu), sigma=Sigma, log=log, ...);

    # Return logs if needed
    if(log){
        cnormReturn[] <- log(cnormReturn);
    }

    return(cnormReturn);
}

#' @rdname cnormal
#' @export pcnorm
#' @aliases pcnorm
pcnorm <- function(lower=-Inf, upper=Inf, mu=0+0i, sigma=1, varsigma=0+0i, ...){

    # Create a covariance matrix based on the provided variances
    Sigma <- matrix(c((sigma+Re(varsigma)),Im(varsigma),Im(varsigma),(sigma-Re(varsigma)))/2, 2, 2);

    cnormReturn <- pmvnorm(lower=lower, upper=upper, mean=complex2vec(mu), sigma=Sigma, ...);

    return(cnormReturn);
}

#' @rdname cnormal
#' @export qcnorm
#' @aliases qcnorm
qcnorm <- function(p, mu=0+0i, sigma=1, varsigma=0+0i, ...){

    # Create a covariance matrix based on the provided variances
    Sigma <- matrix(c((sigma+Re(varsigma)),Im(varsigma),Im(varsigma),(sigma-Re(varsigma)))/2, 2, 2);

    cnormReturn <- qmvnorm(p, mean=as.vector(complex2vec(mu)), sigma=Sigma, ...);

    return(cnormReturn);
}

#' @rdname cnormal
#' @export rcnorm
#' @aliases rcnorm
rcnorm <- function(n=1, mu=0+0i, sigma=1, varsigma=0+0i, ...){

    # Create a covariance matrix based on the provided variances
    Sigma <- matrix(c((sigma+Re(varsigma)),Im(varsigma),Im(varsigma),(sigma-Re(varsigma)))/2, 2, 2);

    cnormReturn <- rmvnorm(n, mean=complex2vec(mu), sigma=Sigma, ...)

    return(vec2complex(cnormReturn));
}
