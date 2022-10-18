#' Pseudo Correlation, Variance and Covariance (Matrices)
#'
#' Functions return (complex) pseudo-variance, pseudo-covariance and pseudo-correlation based on the
#' provided complex vector/matrix \code{x}. Function \code{covar()} returns the covariance
#' matrix based
#'
#' Functions calculate pseudo statistics for complex variables. Only the parametric
#' pseudo-correlation is supported by the function. If \code{x} is matrix, then \code{y}
#' is ignored.
#'
#' \code{covar()} function returns a covariance matrix calculated for the provided complex
#' vector or matrix.
#'
#' @author Ivan Svetunkov, \email{ivan@svetunkov.ru}
#' @keywords univar
#'
#' @param x vector or matrix of complex variables. If it is matrix then the
#' variable \code{y} is ignored.
#' @param y second vector to calculate pseudo-covariance or pseudo-correlations with.
#' @param ... parameters passed to \code{mean()} functions. For example, this can be
#' \code{na.rm=TRUE} to remove missing values or \code{trim} to define the trimming
#' in the mean (see \link[base]{mean}).
#'
#' @return A scalar or a matrix with resulting complex variables.
#'
#' #' @references \itemize{
#' \item Svetunkov, S. (2022) Complex Autoregressions. In Press.
#' }
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
cvar <- function(x, ...){
    if(is.matrix(x)){
        xMeans <- matrix(colMeans(x, ...),nrow(x),ncol(x),byrow=TRUE);
        return(t(x-xMeans) %*% (x-xMeans));
    }
    else{
        return(mean((x-mean(x, ...))^2, ...));
    }
}

#' @rdname ccor
#' @export
ccov <- function(x, y, ...){
    if(is.matrix(x)){
        return(cvar(x, ...));
    }
    else{
        return(mean((x-mean(x, ...))*(y-mean(y, ...)), ...))
    }
}

#' @rdname ccor
#' @export
ccor <- function(x, y, ...){
    if(is.matrix(x)){
        ccov2cor(cvar(x, ...));
    }
    else{
        return(ccov(x, y, ...) / sqrt(cvar(x, ...)*cvar(y, ...)));
    }
}

#' @param V complex pseudo-covariance matrix.
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
