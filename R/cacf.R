#' Complex Autocorrelation Function Estimation
#'
#' The function computes (and by default plots) estimates of the complex autocovariance
#' or complex autocorrelation function.
#'
#' For \code{type="correlation"} and \code{"covariance"}, the estimates are based
#' on the sample complex covariance and use complex correlation \link[complex]{ccor} and complex
#' covariance \link[complex]{ccov} respectively. Note that the function does not calculate values for
#' lag 0. Also, the function will automatically remove NAs. Finally, function does not have
#' \code{demean} parameter, because \code{ccov()} and \code{ccor()} do that automatically.
#'
#' The generic function plot has a method for objects of class "cacf".
#'
#' The lag is returned and plotted in units of time, and not numbers of observations.
#'
#' There is a print method for objects of class "cacf".
#'
#' @author Ivan Svetunkov, \email{ivan@svetunkov.ru}
#' @keywords univar
#'
#' @param x vector of complex variables.
#' @param lag.max maximum number of lags. See \link[stats]{acf} for more details.
#' @param type character string giving the type of cacf to be computed. Allowed values
#' are "correlation" (the default) and "covariance". Will be partially matched.
#' @param plot logical. If \code{TRUE} (the default) the cacf is plotted on complex plane
#' and as two linear graphs for real and imaginary parts.
#'
#' @return An object of class "cacf", which is a list with the following elements:
#' \itemize{
#' \item \code{lag} A three dimensional array containing the lags at which the cacf is estimated.
#' \item \code{acf} An array with the same dimensions as lag containing the estimated cacf.
#' \item \code{type} The type of correlation (same as the type argument).
#' \item \code{n.used} The number of observations in the time series.
#' \item \code{series} The name of the series x.
#' }
#'
#' #' @references \itemize{
#' \item Svetunkov, S. (2022) Complex Autoregressions. In Press.
#' }
#'
#' @seealso \link[stats]{acf}, \link[complex]{ccor}
#'
#' @examples
#'
#' # Generate random complex variables
#' x <- complex(real=rnorm(100,10,10), imaginary=rnorm(100,10,10))
#'
#' # Calculate CACF
#' cacf(x)
#'
#' @importFrom greybox xregExpander
#' @export
cacf <- function(x, lag.max=NULL, type=c("correlation", "covariance"),
                 plot=TRUE){
    # Function is based on acf() from stats

    type <- match.arg(type);
    series <- deparse1(substitute(x))

    obs <- length(x);

    if(is.null(lag.max)){
        lag.max <- floor(10*log10(obs));
    }
    lag.max <- as.integer(min(lag.max, obs - 1));
    if(is.na(lag.max) || lag.max < 1){
        stop("'lag.max' must be at least 1");
    }

    xLagged <- xregExpander(x, lags=c(-1:-lag.max), gaps="NAs");
    xACF <- vector("complex", lag.max);

    for(i in 1:lag.max){
        xACF[i] <- switch(type,
                          "covariance"=ccov(xLagged[,1], xLagged[,i+1], na.rm=TRUE),
                          "correlation"=ccor(xLagged[,1], xLagged[,i+1], na.rm=TRUE));
    }

    acf.out <- structure(list(acf=xACF, type=type, n.used=obs,
                              lag=c(1:lag.max), series=series), class="cacf");

    mainLabel <- switch(type,
                        "covariance"="Autocovariance function",
                        "correlation"="Autocorrelation function");

    if(plot){
        plot(acf.out, 1);
        invisible(acf.out);
    }
    else{
        return(acf.out);
    }
}

#' @importFrom stats setNames
#' @export
print.cacf <- function(x, ...){
    cat("Complex Autocorrelations of series", x$series, "by lag\n\n");
    x$acf |> setNames(x$lag) |> print();
}

#' @importFrom graphics layout text devAskNewPage
#' @export
plot.cacf <- function(x, which=c(1,2), ask=length(which)>1, ...){
    mainLabel <- switch(x$type,
                        "covariance"="Complex autocovariance function",
                        "correlation"="Complex autocorrelation function");

    # Define, whether to wait for the hit of "Enter"
    if(ask){
        devAskNewPage(TRUE);
        on.exit(devAskNewPage(FALSE));
    }

    for(i in which){
        if(i==1){
            parDefault <- par(no.readonly=TRUE);
            on.exit(par(parDefault),add=TRUE);

            layout(matrix(c(1,2,1,3),2,2))
            par(mar=c(4,4,4,2))
            plot(x$acf, type="l",
                 xlab="Real CACF", ylab="Imaginary CACF", main=mainLabel);
            points(x$acf);
            text(Re(x$acf), Im(x$acf), c(1:length(x$acf)), pos=3);

            par(mar=c(4,4,1,2))
            plot(Re(x$acf), type="h", xlab="Lag", ylab="Real CACF");
            abline(h=0);

            plot(Im(x$acf), type="h", xlab="Lag", ylab="Imaginary CACF");
            abline(h=0);
        }
        if(i==2){
            parDefault <- par(no.readonly=TRUE);
            on.exit(par(parDefault),add=TRUE);

            layout(matrix(c(1,2,1,3),2,2))
            par(mar=c(4,4,4,2))
            plot(abs(x$acf), Arg(x$acf), type="l",
                 xlab="Absolute of CACF", ylab="Argument of CACF", main=mainLabel);
            points(abs(x$acf), Arg(x$acf));
            text(abs(x$acf), Arg(x$acf), c(1:length(x$acf)), pos=3);

            par(mar=c(4,4,1,2))
            plot(abs(x$acf), type="h", xlab="Lag", ylab="Absolute of CACF");
            abline(h=0);

            plot(Arg(x$acf), type="h", xlab="Lag", ylab="Argument of CACF");
            abline(h=0);
        }
    }
}