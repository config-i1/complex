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
#' @param type type of measure to calculate. \code{"conjugate"} means that it is based
#' on the multiplication by conjugate number. \code{"direct"} means the calculation
#' without the conjugate (i.e. "pseudo" moment).
#' @param y second vector to calculate pseudo-covariance or pseudo-correlations with.
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
cvar <- function(x, type=c("direct","conjugate"),
                 ...){
    type <- match.arg(type);
    if(type=="direct"){
        if(is.matrix(x)){
            xMeans <- matrix(colMeans(x, ...),nrow(x),ncol(x),byrow=TRUE);
            return(t(x-xMeans) %*% (x-xMeans));
        }
        else{
            return(mean((x-mean(x, ...))^2, ...));
        }
    }
    else{
        if(is.matrix(x)){
            xMeans <- matrix(colMeans(x, ...),nrow(x),ncol(x),byrow=TRUE);
            return(t(x-xMeans) %*% Conj(x-xMeans));
        }
        else{
            return(mean((x-mean(x, ...)) * Conj(x-mean(x, ...)), ...));
        }
    }
}

#' @rdname ccor
#' @export
ccov <- function(x, y, type=c("direct","conjugate"),
                 ...){
    type <- match.arg(type);
    if(type=="direct"){
        if(is.matrix(x)){
            return(cvar(x, type=type, ...));
        }
        else{
            return(mean((x-mean(x, ...))*(y-mean(y, ...)), ...));
        }
    }
    else{
        if(is.matrix(x)){
            return(cvar(x, type=type, ...));
        }
        else{
            return(mean((x-mean(x, ...))*Conj(y-mean(y, ...)), ...));
        }
    }
}

#' @rdname ccor
#' @export
ccor <- function(x, y, type=c("direct","conjugate"),
                 ...){
    type <- match.arg(type);
    if(is.matrix(x)){
        ccov2cor(cvar(x, type=type, ...));
    }
    else{
        return(ccov(x, y, type=type, ...) / sqrt(cvar(x, type=type, ...)*cvar(y, type=type, ...)));
    }
}

#' @param V complex (pseudo)covariance matrix.
#' @rdname ccor
#' @export
ccov2cor <- function(V){
    return(V / sqrt(diag(V) %*% t(diag(V))));
}


#' @param df Number of degrees of freedom to use in the calculation of the covariance matrix.
#' @rdname ccor
#' @export
covar <- function(x, df=NULL){
    obs <- length(x);

    if(is.null(df)){
        df <- obs-1;
    }

    if(is.complex(x)){
        x <- complex2vec(x);
    }

    return(t(x) %*% x / df);
}
