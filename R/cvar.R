#' Correlation, Variance and Covariance (Matrices) for complex variables
#'
#' Functions \code{cvar()}, \code{ccov()} and \code{ccor()} return respectively
#' complex variance, covariance and correlation based on the
#' provided complex vector/matrix \code{x}. Function \code{covar()} returns the covariance
#' matrix based on a complex vector/matrix.
#'
#' Only the parametric correlation is supported by the function. If \code{x}
#' is matrix, then \code{y} is ignored.
#'
#' \code{covar()} function returns a covariance matrix calculated for the provided complex
#' vector or matrix \code{x}.
#'
#' @template references
#' @template author
#'
#' @keywords univar
#'
#' @param x vector or matrix of complex variables. If it is matrix then the
#' variable \code{y} is ignored.
#' @param y second vector to calculate covariance or correlations with.
#' @param kind kind of measure to calculate. \code{"conjugate"} means that it is based
#' on the multiplication by conjugate number. \code{"direct"} means the calculation
#' without the conjugate (i.e. "pseudo" moment). For \code{ccor} the variable \code{kind}
#' can also be "pearson", "kendall", or "spearman", defining what correlation coefficient
#' to use after MDS transformation of complex variable \code{x} and \code{y}.
#' @param ... parameters passed to \code{mean()} functions. For example, this can be
#' \code{na.rm=TRUE} to remove missing values or \code{trim} to define the trimming
#' in the mean (see \link[base]{mean}).
#'
#' @return A scalar or a matrix with resulting complex variables.
#'
#' @seealso \code{\link[stats]{cor}}
#'
#' @examples
#'
#' # Generate random complex variables
#' x <- complex(real=rnorm(100,10,10), imaginary=rnorm(100,10,10))
#' y <- complex(real=rnorm(100,10,10), imaginary=rnorm(100,10,10))
#'
#' # Create a matrix of complex variables
#' z <- cbind(x,y)
#'
#' # Calculate measures
#' cvar(x)
#' cvar(z)
#' ccor(x,y)
#' ccor(z)
#'
#' @rdname ccor
#' @export
cvar <- function(x, kind=c("direct","conjugate"),
                 df=NULL, ...){
    kind <- match.arg(kind);
    if(kind=="direct"){
        if(is.matrix(x)){
            if(is.null(df)){
                df <- nrow(x)-1;
            }
            xMeans <- matrix(colMeans(x, ...),nrow(x),ncol(x),byrow=TRUE);
            return(t(x-xMeans) %*% (x-xMeans) / df);
        }
        else{
            if(is.null(df)){
                df <- length(x)-1;
            }
            return(sum((x-mean(x, ...))^2, ...) / df);
        }
    }
    else{
        if(is.matrix(x)){
            if(is.null(df)){
                df <- nrow(x)-1;
            }
            xMeans <- matrix(colMeans(x, ...),nrow(x),ncol(x),byrow=TRUE);
            return(Conj(t(x-xMeans)) %*% (x-xMeans) / df);
        }
        else{
            if(is.null(df)){
                df <- length(x)-1;
            }
            return(sum(Conj(x-mean(x, ...)) * (x-mean(x, ...)), ...) / df);
        }
    }
}

#' @rdname ccor
#' @export
ccov <- function(x, y, kind=c("direct","conjugate"),
                 df=NULL, ...){
    kind <- match.arg(kind);
    if(kind=="direct"){
        if(is.matrix(x)){
            return(cvar(x, kind=kind, ...));
        }
        else{
            if(is.null(df)){
                df <- length(x)-1;
            }
            return(sum((x-mean(x, ...))*(y-mean(y, ...)), ...) / df);
        }
    }
    else{
        if(is.matrix(x)){
            return(cvar(x, kind=kind, ...));
        }
        else{
            if(is.null(df)){
                df <- length(x)-1;
            }
            return(sum(Conj(x-mean(x, ...))*(y-mean(y, ...)), ...) / df);
        }
    }
}

#' @rdname ccor
#' @importFrom stats cor
#' @export
ccor <- function(x, y, kind=c("direct","conjugate","pearson","kendall", "spearman"),
                 ...){
    kind <- match.arg(kind);
    if(any(kind==c("direct","conjugate"))){
        if(is.matrix(x)){
            ccov2cor(cvar(x, kind=kind, ...));
        }
        else{
            return(sqrt((ccov(x, y, kind=kind, ...) * ccov(y, x, kind=kind, ...)) / (cvar(x, kind=kind, ...)*cvar(y, kind=kind, ...))));
        }
    }
    else{
        xScaled <- cmdscale(dist(complex2vec(x)), k=1)
        yScaled <- cmdscale(dist(complex2vec(y)), k=1)
        cor(xScaled, yScaled, method=kind, ...)
    }
}

#' @param V complex (pseudo)covariance matrix.
#' @rdname ccor
#' @export
ccov2cor <- function(V){
    return(V / sqrt(diag(V) %*% t(diag(V))));
}


#' @param df Number of degrees of freedom to use in the calculation of the statistics.
#' @rdname ccor
#' @export
covar <- function(x, df=NULL){
    if(is.complex(x)){
        x <- complex2vec(x);
    }

    # Centre the variable
    if(is.matrix(x)){
        obs <- nrow(x);
        if(is.null(df)){
            df <- obs-1;
        }
        x[] <- x - matrix(colMeans(x), nrow=obs, ncol=ncol(x), byrow=TRUE);
    }
    else{
        if(is.null(df)){
            df <- length(x)-1;
        }
        x[] <- x - mean(x);
    }

    return(t(x) %*% x / df);
}
