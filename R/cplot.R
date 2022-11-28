#' Scatterplots for complex variables
#'
#' Function produces six scatterplots to show relations between the two complex variables x and y.
#'
#' The plots are positioned to satisfy two rules:
#' 1. When a scatterplot for a c.r.v. is produced, the real part should be in x-axis,
#' while the imaginary should be in the y-axis.
#' 2. When parts of variables x and y are compared, the part for $x$ should be in x-axis,
#' while the part for y should be in y-axis, which should the reflect the idea that x could
#'  be an explanatory variable for y.
#'
#' @template references
#' @template author
#'
#' @keywords univar
#'
#' @param x vector of a complex variable.
#' @param y second vector of a complex variable.
#'
#' @return The function produces a plot and does not return any value
#'
#' @seealso \code{\link[complex]{ccor}}
#'
#' @examples
#'
#' # Generate random complex variables
#' x <- complex(real=rnorm(100,10,10), imaginary=rnorm(100,10,10))
#' y <- complex(real=rnorm(100,10,10), imaginary=rnorm(100,10,10))
#'
#' cplot(x, y)
#'
#' @importFrom graphics axis box mtext
#' @importFrom grDevices grey palette
#' @export
cplot <- function(x, y){
    # The function creates scatterplots for two complex variables, x and y

    obs <- length(x)
    if(length(y)!=obs){
        stop("The length of x and y differs, cannot proceed", call.=FALSE);
    }

    ourData <- data.frame(x=x, y=y)
    ourData <- ourData[order(abs(ourData$x-complex(real=min(Re(ourData$x)), imaginary=min(Im(ourData$x))))),]

    parDefault <- par(no.readonly=TRUE);
    on.exit(par(parDefault));

    nColours <- min(1000,obs);
    nColoursTimes <- floor(obs/nColours);
    colours <- rep(grey(seq(0.9, 0.1, len=nColours)), each=nColoursTimes);
    omaValues <- c(5,5,2,2);

    par(mfcol=c(3,3), mar=rep(0,4), oma=omaValues, xaxt="s",yaxt="s",cex.main=1.5);

    plot(Re(ourData$x),Im(ourData$x), col=colours, xlab="Re(x)", ylab="Im(x)", axes=F)
    axis(side=2)
    box(col="darkred", lwd=2)
    mtext("Im(x)", side=2, line=3, at=2.5/3, outer=TRUE)
    plot(Re(ourData$x),Re(ourData$y), col=colours, xlab="Re(x)", ylab="Re(y)", axes=F)
    axis(side=2)
    box(col="darkblue", lwd=2)
    mtext("Re(y)", side=2, line=3, outer=TRUE)
    plot(Re(ourData$x),Im(ourData$y), col=colours, xlab="Re(x)", ylab="Im(y)", axes=F)
    axis(side=2)
    axis(side=1)
    box(col="darkgreen", lwd=2)
    mtext("Im(y)", side=2, line=3, at=0.5/3, outer=TRUE)
    mtext("Re(x)", side=1, line=3, at=0.5/3, outer=TRUE)

    plot(0,0, col="white", xlab="", ylab="", axes=F)
    # box(col="black", lwd=1)
    plot(Im(ourData$x),Re(ourData$y), col=colours, ylab="Re(y)", xlab="Im(x)", axes=F)
    box(col="darkgreen", lwd=2)
    plot(Im(ourData$x),Im(ourData$y), col=colours, xlab="Im(x)", ylab="Im(y)", axes=F)
    axis(side=1)
    box(col="darkblue", lwd=2)
    mtext("Im(x)", side=1, line=3, outer=TRUE)

    plot(0,0, col="white", xlab="", ylab="", axes=F)
    plot(0,0, col="white", xlab="", ylab="", axes=F)
    plot(Re(ourData$y),Im(ourData$y), col=colours, xlab="Re(y)", ylab="Im(y)", axes=F)
    axis(side=1)
    box(col="darkred", lwd=2)
    mtext("Re(y)", side=1, line=3, at=2.5/3, outer=TRUE)
}
