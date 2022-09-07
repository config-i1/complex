#' Complex Correlation, Variance and Covariance (Matrices)
#'
#' Functions return complex variance, covariance and correlation based on the
#' provided complex vector/matrix \code{x}.
#'
#' Functions calculate complex statistics for complex variables. Only the parametric
#' correlation is supported by the function. If \code{x} is matrix, then \code{y}
#' is ignored.
#'
#' @author Ivan Svetunkov, \email{ivan@svetunkov.ru}
#' @keywords univar
#'
#' @param x Vector or matrix of complex variables. If it is matrix then the
#' variable \code{y} is ignored.
#' @param y Second vector to calculate covariance or correlations with.
#' @param V Complex covariance matrix.
#' @param ... Parameters passed to \code{mean()} functions. For example, this can be
#' \code{na.rm=TRUE} to remove missing values.
#'
#' @return A scalar or a matrix with resulting complex variables.
#' @seealso \code{\link[stats]{cor}}
#'
#' @examples
#'
#' # Generat random complex variables
#' x <- complex(real=rnorm(100,10,10), imaginary=rnorm(100,10,10))
#' y <- complex(real=rnorm(100,10,10), imaginary=rnorm(100,10,10))
#'
#' # Create a matrix of complex variables
#' z <- cbind(x,y)
#'
#' # Calculate measures
#' cvar(x,y)
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
        ccov2ccor(cvar(x, ...));
    }
    else{
        return(ccov(x, y, ...) / sqrt(cvar(x, ...)*cvar(y, ...)));
    }
}

#' @rdname ccor
#' @export
ccov2ccor <- function(V){
    return(V / sqrt(diag(V) %*% t(diag(V))));
}