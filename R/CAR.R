#' Complex Autoregression
#'
#' Function constructs Complex Autoregression model as a restricted VAR for a
#' provided pair time series.
#'
#' The CAR(p) model can be written as a VAR(p) model with the following
#' restrictions on the proportionality matrices:
#' \deqn{
#' A_j = \begin{matrix} a_{j,1} & -a_{j,2} \\ a_{j,2} & a_{j,1} \end{matrix}
#' }{
#' A_j =
#' [ a_{j,1} | -a_{j,2} ]
#' [ a_{j,2} |  a_{j,1} ]
#' }
#'
#' @author Ivan Svetunkov, \email{ivan@svetunkov.ru}
#' @keywords univar ts models regression
#'
#' @param y The pair of the time series (the matrix) in columns.
#' @param order The order p of the model. If \code{NULL}, then it is selected
#' automatically based on the selected information criterion.
#' @param h Length of forecasting horizon.
#' @param holdout If \code{TRUE}, holdout sample of size \code{h} is taken from
#' the end of the data.
#' nothing happens.
#' @param silent If \code{TRUE}, then the plot is produced. Otherwise nothing happens.
#' @param ic The information criterion used in the model selection procedure.
#' @param restrict If \code{FALSE}, then VAR(p) is constructed instead of CAR(p).
#'
#' @return Object of class "vsmooth" is returned. It contains the following list of
#' values:
#' \itemize{
#' \item \code{model} - The name of the fitted model;
#' \item \code{timeElapsed} - The time elapsed for the construction of the model;
#' \item \code{coefficients} - The vector of all the estimated coefficients;
#' \item \code{nParam} - The number of estimated parameters;
#' \item \code{y} - The matrix with the original data;
#' \item \code{fitted} - The matrix of the fitted values;
#' \item \code{holdout} - The matrix with the holdout values (if \code{holdout=TRUE} in
#' the estimation);
#' \item \code{residuals} - The matrix of the residuals of the model;
#' \item \code{forecast} - The matrix of point forecasts;
#' \item \code{Sigma} - The covariance matrix of the errors (estimated with the correction
#' for the number of degrees of freedom);
#' \item \code{ICs} - The values of the information criteria for the selected model;
#' \item \code{logLik} - The log-likelihood function;
#' \item \code{lossValue} - The value of the loss function;
#' \item \code{loss} - The type of the used loss function;
#' \item \code{accuracy} - the values of the error measures. Currently not available.
#' \item \code{ICsAll} - In case of the order selection, the values of information
#' criteria for the tested models.
#' }
#' @seealso \code{\link[legion]{ves}}
#'
#' @examples
#'
#' Y <- ts(cbind(rnorm(100,100,10),rnorm(100,75,8)),frequency=4)
#'
#' # The simplest model applied to the data with the default values
#' CARModel <- CAR(Y,h=10,holdout=TRUE)
#'
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom nloptr nloptr
#' @importFrom greybox AICc BICc measures
#' @importFrom graphics abline lines par plot points
#' @importFrom stats AIC BIC deltat frequency start time ts fitted
#' @import smooth
#' @import legion
#'
#' @export
CAR <- function(y, order=NULL, h=13, holdout=TRUE, silent=TRUE,
                ic=c("AIC","AICc","BIC","BICc"), restrict=TRUE){
    # The function fits the restricted VAR, which corresponds to CAR. The likelihood is multivariate normal.
    # restrict defines whether the CAR models is constructed or a basic VAR

    # Start measuring the time of calculations
    startTime <- Sys.time();

    IC <- switch(ic[1],
                 "AIC"=AIC,
                 "AICc"=AICc,
                 "BIC"=BIC,
                 "BICc"=BICc);


    if(!is.null(ncol(y)) && (ncol(y)>2 || ncol(y)==1)){
        stop("This function can only work with two series.", call.=FALSE);
    }

    # Prepare the sample
    obsInSample <- nrow(y) - h * holdout;
    obsAll <- nrow(y) + h * (!holdout);
    nSeries <- 2;

    # Prepare the data and the matrices
    yInSample <- y[1:obsInSample,];
    yFitted <- matrix(0,obsInSample,nSeries);
    yForecast <- matrix(NA,h,nSeries);
    errors <- matrix(NA,obsInSample,nSeries);

    # Write down the names and ts structure
    dataNames <- colnames(y);
    dataFreq <- frequency(y);
    dataDeltat <- deltat(y);
    dataStart <- start(y);
    yForecastStart <- time(y)[obsInSample]+deltat(y);
    colnames(yFitted) <- colnames(yForecast) <- colnames(errors) <- dataNames;

    ##### If order needs to be selected #####
    if(is.null(order)){
        CARmodels <- vector("list",dataFreq);
        if(!silent){
            cat("Estimation progress:    ");
        }

        for(i in 1:dataFreq){
            if(!silent){
                if(i==1){
                    cat("\b");
                }
                cat(paste0(rep("\b",nchar(round((i-1)/dataFreq,2)*100)+1),collapse=""));
                cat(paste0(round(i/dataFreq,2)*100,"%"));
            }
            CARmodels[[i]] <- CAR(y=y, order=i, h=h, holdout=holdout, restrict=TRUE);
        }

        if(!silent){
            cat("... Done! \n");
        }

        CARICs <- sapply(CARmodels,IC);
        names(CARICs) <- paste0("CAR(",c(1:dataFreq),")");
        i <- which.min(CARICs);
        CARmodels[[i]]$ICsAll <- CARICs;

        yFitted <- fitted(CARmodels[[i]]);
        yForecast <- CARmodels[[i]]$forecast;

        model <- CARmodels[[i]];
        modelname <- model$model;
    }
    ##### If a specific order is provided #####
    else{
        # A is the matrix of 2 x order with the last row for the
        CARFitter <- function(A, yInSample, order){
            yFitted[] <- 0;
            for(i in (order+1):obsInSample){
                for(j in 1:order){
                    yFitted[i,1] <- yFitted[i,1] + yInSample[i-j,] %*% A[,(j-1)*nSeries+1];
                    yFitted[i,2] <- yFitted[i,2] + yInSample[i-j,] %*% A[,j*nSeries];
                }
            }
            yFitted[,1] <- yFitted[,1] + A[1,order*nSeries+1]
            yFitted[,2] <- yFitted[,2] + A[2,order*nSeries+1]

            return(yFitted);
        }

        CARCF <- function(A, yInSample, order){
            if(restrict){
                A <- matrix(A,nrow=nSeries);
                B <- matrix(0, nSeries, order*nSeries);
                for(i in 1:order){
                    B[,(i-1)*nSeries+1] <- A[,i];
                    B[,i*nSeries] <- rev(A[,i]) * c(-1, 1);
                }
                A <- cbind(B, A[,order+1]);
            }
            else{
                A <- matrix(A,nrow=nSeries);
            }
            yFitted[] <- CARFitter(A,yInSample,order);
            sigma <- t(yInSample - yFitted) %*% (yInSample - yFitted);
            logLikReturned <- sum(dmvnorm(yInSample - yFitted, sigma=sigma, log=TRUE));

            return(-logLikReturned);
        }

        CARForecaster <- function(A, yInSample, order){
            yBuffer <- matrix(0, h+order, nSeries);
            yBuffer[1:order,] <- yInSample[obsInSample-c(1:order)+1,]
            for(i in 1:h){
                for(j in 1:order){
                    yBuffer[i+order,1] <- yBuffer[i+order-j,] %*% A[,(j-1)*nSeries+1];
                    yBuffer[i+order,2] <- yBuffer[i+order-j,] %*% A[,j*nSeries];
                }
            }
            yBuffer[,1] <- yBuffer[,1] + A[1,order*nSeries+1]
            yBuffer[,2] <- yBuffer[,2] + A[2,order*nSeries+1]

            return(yBuffer[-c(1:order),]);
        }

        ##### Define the initials and start the optimisation #####
        # AR terms and constant
        if(restrict){
            A <- rep(0.5, order*nSeries);
        }
        else{
            A <- rep(0.5, order*nSeries^2);
        }
        # Define constants
        A <- c(A, colMeans(yInSample));

        res <- nloptr(A, CARCF,# lb=AList$ALower, ub=AList$AUpper,
                      opts=list("algorithm"="NLOPT_LN_BOBYQA", "xtol_rel"=1e-8, "maxeval"=1000, print_level=0),
                      yInSample=yInSample, order=order);
        A <- res$solution;

        if(restrict){
            A <- matrix(A,nrow=nSeries);
            B <- matrix(0, nSeries, order*nSeries);
            for(i in 1:order){
                B[,(i-1)*nSeries+1] <- A[,i];
                B[,i*nSeries] <- rev(A[,i]) * c(-1, 1);
            }
            A <- cbind(B, A[,order+1]);
        }
        else{
            A <- matrix(A,nrow=nSeries);
        }

        logLik <- -res$objective;

        #### Produce fitted and forecasts ####
        yFitted[] <- CARFitter(A,yInSample,order);
        yFitted[1:order,] <- NA;
        errors[] <- yInSample - yFitted;
        yForecast[] <- CARForecaster(A,yInSample,order);

        yForecast <- ts(yForecast,start=time(y)[obsInSample] + dataDeltat,frequency=dataFreq);
        yFitted <- ts(yFitted,start=dataStart,frequency=dataFreq);
        yInSample <- ts(y,start=dataStart,frequency=dataFreq);
        errors <- ts(errors,start=dataStart,frequency=dataFreq);
        if(holdout){
            yHoldout <- ts(y[-c(1:obsInSample),],start=time(y)[obsInSample] + dataDeltat,frequency=dataFreq);
        }
        else{
            yHoldout <- NA;
        }

        # Number of parameters to estimate / provided
        parametersNumber <- matrix(0,2,4,
                                   dimnames=list(c("Estimated","Provided"),
                                                 c("nParamInternal","nParamXreg",
                                                   "nParamIntermittent","nParamAll")));
        # orders + constant + cov matrix
        parametersNumber[1,4] <- parametersNumber[1,1] <- 2*order+2+3;

        modelname <- paste0("CAR(",order,")");

        ##### Now let's deal with the holdout #####
        if(holdout){
            measureFirst <- measures(yHoldout[,1],yForecast[,1],yInSample[1,]);
            errorMeasures <- matrix(NA,nSeries,length(measureFirst));
            rownames(errorMeasures) <- dataNames;
            colnames(errorMeasures) <- names(measureFirst);
            errorMeasures[1,] <- measureFirst;
            for(i in 2:nSeries){
                errorMeasures[i,] <- measures(yHoldout[,i],yForecast[,i],yInSample[i,]);
            }
        }
        else{
            errorMeasures <- NA;
        }

        model <- list(model=modelname, timeElapsed=Sys.time()-startTime,
                      y=yInSample, fitted=yFitted, coefficients=A, residuals=errors, forecast=yForecast,
                      nParam=parametersNumber, logLik=logLik, holdout=yHoldout, PI=NA, loss="likelihood",
                      lossValue=-logLik, logLik=logLik, Sigma=1/obsInSample * sum(t(errors) %*% errors),
                      accuracy=errorMeasures);
        model <- structure(model,class=c("vsmooth","smooth"));
        model$ICs <- c(AIC(model),AICc(model),BIC(model),BICc(model));
        names(model$ICs) <- c("AIC","AICc","BIC","BICc");
    }

    #### Produce a graph if needed ####
    if(!silent){
        pages <- ceiling(nSeries / 5);
        perPage <- ceiling(nSeries / pages);
        packs <- seq(1, nSeries+1, perPage);
        if(packs[length(packs)]<nSeries+1){
            packs <- c(packs,nSeries+1);
        }
        parDefault <- par(no.readonly=TRUE);
        for(j in 1:pages){
            par(mar=c(4,4,2,1),mfcol=c(perPage,1));
            for(i in packs[j]:(packs[j+1]-1)){
                # if(any(intervalType==c("u","i"))){
                #     plotRange <- range(min(y[,i],yForecast[,i],yFitted[,i],PI[,i*2-1]),
                #                        max(y[,i],yForecast[,i],yFitted[,i],PI[,i*2]));
                # }
                # else{
                plotRange <- range(min(y[,i],yForecast[,i],yFitted[,i],na.rm=TRUE),
                                   max(y[,i],yForecast[,i],yFitted[,i],na.rm=TRUE));
                # }
                plot(y[,i],main=paste0(modelname," ",dataNames[i]),ylab="Y",
                     ylim=plotRange, xlim=range(time(y[,i])[1],time(yForecast)[max(h,1)]),
                     type="l");
                lines(yFitted[,i],col="purple",lwd=2,lty=2);
                if(h>1){
                    # if(any(intervalType==c("u","i"))){
                    #     lines(PI[,i*2-1],col="darkgrey",lwd=3,lty=2);
                    #     lines(PI[,i*2],col="darkgrey",lwd=3,lty=2);
                    #
                    #     polygon(c(seq(dataDeltat*(yForecastStart[2]-1)+yForecastStart[1],
                    #                   dataDeltat*(end(yForecast)[2]-1)+end(yForecast)[1],dataDeltat),
                    #               rev(seq(dataDeltat*(yForecastStart[2]-1)+yForecastStart[1],
                    #                       dataDeltat*(end(yForecast)[2]-1)+end(yForecast)[1],dataDeltat))),
                    #             c(as.vector(PI[,i*2]), rev(as.vector(PI[,i*2-1]))), col="lightgray",
                    #             border=NA, density=10);
                    # }
                    lines(yForecast[,i],col="blue",lwd=2);
                }
                else{
                    # if(any(intervalType==c("u","i"))){
                    #     points(PI[,i*2-1],col="darkgrey",lwd=3,pch=4);
                    #     points(PI[,i*2],col="darkgrey",lwd=3,pch=4);
                    # }
                    points(yForecast[,i],col="blue",lwd=2,pch=4);
                }
                abline(v=dataDeltat*(yForecastStart[2]-2)+yForecastStart[1],col="red",lwd=2);
            }
        }
        par(parDefault);
    }
    return(model);
}
