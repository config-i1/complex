#' Complex Linear Model
#'
#' Function estimates complex variables model
#'
#' This is a function, similar to \link[stats]{lm}, but supporting several estimation
#' techniques for complex variables regression.
#'
#' @template author
#' @template keywords
#' @template references
#'
#' @param formula an object of class "formula" (or one that can be coerced to
#' that class): a symbolic description of the model to be fitted. Can also include
#' \code{trend}, which would add the global trend.
#' @param data a data frame or a matrix, containing the variables in the model.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param na.action	a function which indicates what should happen when the
#' data contain NAs. The default is set by the na.action setting of
#' \link[base]{options}, and is \link[stats]{na.fail} if that is unset. The
#' factory-fresh default is \link[stats]{na.omit}. Another possible value
#' is NULL, no action. Value \link[stats]{na.exclude} can be useful.
#' @param loss The type of Loss Function used in optimization. \code{loss} can
#' be:
#' \itemize{
#' \item \code{OLS} - Ordinary Least Squares method, relying on the minimisation of
#' the conjoint variance of the error term;
#' \item \code{CLS} - Complex Least Squares method, relying on the minimisation of
#' the complex variance of the error term;
#' \item \code{likelihood} - the model is estimated via the maximisation of the
#' likelihood of the complex Normal distribution;
#' \item \code{MSE} (Mean Squared Error),
#' \item \code{MAE} (Mean Absolute Error),
#' \item \code{HAM} (Half Absolute Moment),
#' }
#'
#' A user can also provide their own function here as well, making sure
#' that it accepts parameters \code{actual}, \code{fitted} and \code{B}. Here is an
#' example:
#'
#' \code{lossFunction <- function(actual, fitted, B, xreg) return(mean(abs(actual-fitted)))}
#' \code{loss=lossFunction}
#'
#' @param parameters vector of parameters of the linear model. When \code{NULL}, it
#' is estimated.
#' @param fast if \code{TRUE}, then the function won't check whether
#' the data has variability and whether the regressors are correlated. Might
#' cause trouble, especially in cases of multicollinearity.
#' @param ... additional parameters to pass to distribution functions. This
#' includes:
#' \itemize{
#' \item \code{FI=TRUE} will make the function also produce Fisher Information
#' matrix, which then can be used to calculated variances of smoothing parameters
#' and initial states of the model. This is used in the \link[stats]{vcov} method;
#' }
#'
#' You can also pass parameters to the optimiser:
#' \enumerate{
#' \item \code{B} - the vector of starting values of parameters for the optimiser,
#' should correspond to the explanatory variables. If formula for scale was provided,
#' the parameters for that part should follow the parameters for location;
#' \item \code{algorithm} - the algorithm to use in optimisation
#' (\code{"NLOPT_LN_SBPLX"} by default);
#' \item \code{maxeval} - maximum number of evaluations to carry out. Default is 40 per
#' estimated parameter. In case of LASSO / RIDGE the default is 80 per estimated parameter;
#' \item \code{maxtime} - stop, when the optimisation time (in seconds) exceeds this;
#' \item \code{xtol_rel} - the precision of the optimiser (the default is 1E-6);
#' \item \code{xtol_abs} - the absolute precision of the optimiser (the default is 1E-8);
#' \item \code{ftol_rel} - the stopping criterion in case of the relative change in the loss
#' function (the default is 1E-4);
#' \item \code{ftol_abs} - the stopping criterion in case of the absolute change in the loss
#' function (the default is 0 - not used);
#' \item \code{print_level} - the level of output for the optimiser (0 by default).
#' If equal to 41, then the detailed results of the optimisation are returned.
#' }
#' You can read more about these parameters by running the function
#' \link[nloptr]{nloptr.print.options}.
#'
#' @return Function returns \code{model} - the final model of the class
#' "clm", which contains:
#' \itemize{
#' \item coefficients - estimated parameters of the model,
#' \item FI - Fisher Information of parameters of the model. Returned only when \code{FI=TRUE},
#' \item fitted - fitted values,
#' \item residuals - residuals of the model,
#' \item mu - the estimated location parameter of the distribution,
#' \item scale - the estimated scale parameter of the distribution. If a formula was provided for
#' scale, then an object of class "scale" will be returned.
#' \item logLik - log-likelihood of the model. Only returned, when \code{loss="likelihood"}
#' and in a special case of complex least squares.
#' \item loss - the type of the loss function used in the estimation,
#' \item lossFunction - the loss function, if the custom is provided by the user,
#' \item lossValue - the value of the loss function,
#' \item df.residual - number of degrees of freedom of the residuals of the model,
#' \item df - number of degrees of freedom of the model,
#' \item call - how the model was called,
#' \item rank - rank of the model,
#' \item data - data used for the model construction,
#' \item terms - terms of the data. Needed for some additional methods to work,
#' \item B - the value of the optimised parameters. Typically, this is a duplicate of coefficients,
#' \item other - the list of all the other parameters either passed to the
#' function or estimated in the process, but not included in the standard output
#' (e.g. \code{alpha} for Asymmetric Laplace),
#' \item timeElapsed - the time elapsed for the estimation of the model.
#' }
#'
#' @seealso \code{\link[greybox]{alm}}
#'
#' @examples
#'
#' ### An example with mtcars data and factors
#' x <- complex(real=rnorm(1000,10,10), imaginary=rnorm(1000,10,10))
#' a0 <- 10 + 15i
#' a1 <- 2-1.5i
#' y <- a0 + a1 * x + 1.5*complex(real=rnorm(length(x),0,1), imaginary=rnorm(length(x),0,1))
#'
#' complexData <- cbind(y=y,x=x)
#' complexModel <- clm(y~x, complexData)
#' summary(complexModel)
#'
#' @importFrom nloptr nloptr
#' @importFrom stats model.frame sd terms model.matrix update.formula as.formula
#' @importFrom stats formula residuals sigma
#' @rdname clm
#' @export clm
clm <- function(formula, data, subset, na.action,
                loss=c("OLS","CLS","likelihood","MSE","MAE","HAM"),
                parameters=NULL, fast=FALSE, ...){
    # Start measuring the time of calculations
    startTime <- Sys.time();

    # Create substitute and remove the original data
    dataSubstitute <- substitute(data);

    cl <- match.call();
    # This is needed in order to have a reasonable formula saved, so that there are no issues with it
    cl$formula <- eval(cl$formula);
    if(is.function(loss)){
        lossFunction <- loss;
        loss <- "custom";
    }
    else{
        lossFunction <- NULL;
        loss <- match.arg(loss);
    }

    #### Functions used in the estimation ####
    ifelseFast <- function(condition, yes, no){
        if(condition){
            return(yes);
        }
        else{
            return(no);
        }
    }

    # Function that calculates mean. Works faster than mean()
    meanFast <- function(x){
        return(sum(x) / length(x));
    }

    # Basic fitter for non-dynamic models
    fitter <- function(B, y, matrixXreg){
        mu[] <- matrixXreg %*% B

        # Get the scale value
        if(loss=="CLS"){
            scale <- sqrt(sum((y-mu)^2)/obsInsample)
        }
        else{
            scale <- Re(sqrt(sum(Conj((y-mu))*(y-mu))/obsInsample));
        }

        return(list(mu=mu,scale=scale));
    }

    ### Fitted values in the scale of the original variable
    extractorFitted <- function(mu, scale){
        return(mu);
    }

    ### Error term in the transformed scale
    extractorResiduals <- function(mu, yFitted){
        return(y - mu);
    }

    CF <- function(B, loss, y, matrixXreg){

        if(loss=="likelihood"){
            B <- complex(real=B[1:(nVariables/2)],imaginary=B[(nVariables/2+1):nVariables]);
            fitterReturn <- fitter(B, y, matrixXreg);
            errors <- complex2vec(y - fitterReturn$mu);
            sigmaMat <- covar(errors, df=obsInsample);
            # # Concentrated logLik
            # CFValue <- obsInsample*(log(2*pi) + 1 + 0.5*log(det(sigmaMat)));
            # This is the correct one - concentrated does not work for whatever reason...
            CFValue <- obsInsample*(log(2*pi) + 0.5*log(det(sigmaMat))) +
                0.5*sum(diag(errors %*% Re(invert(sigmaMat)) %*% t(errors)));
            ### Likelihood using dcnorm
            # errors <- y - fitterReturn$mu;
            # sigma2 <- mean(errors * Conj(errors));
            # varsigma2 <- mean(errors^2);
            # CFValue <- -sum(dcnorm(errors, mu=0, sigma2=sigma2, varsigma2=varsigma2, log=TRUE));
            ### Likelihood using MVNorm
            # CFValue <- -sum(dmvnorm(errors, mean=c(0,0), sigma=sigmaMat, log=TRUE));
        }
        else if(loss=="CLS"){
            fitterReturn <- fitter(B, y, matrixXreg);
            CFValue <- sum((y - fitterReturn$mu)^2);
        }
        else if(loss=="OLS"){
            fitterReturn <- fitter(B, y, matrixXreg);
            CFValue <- sum(abs(y - fitterReturn$mu)^2);
        }
        else{
            B <- complex(real=B[1:(nVariables/2)],imaginary=B[(nVariables/2+1):nVariables]);
            fitterReturn <- fitter(B, y, matrixXreg);
            yFitted[] <- extractorFitted(fitterReturn$mu, fitterReturn$scale);

            if(loss=="MSE"){
                CFValue <- meanFast((y-yFitted)^2);
            }
            else if(loss=="MAE"){
                CFValue <- meanFast(abs(y-yFitted));
            }
            else if(loss=="HAM"){
                CFValue <- meanFast(sqrt(abs(y-yFitted)));
            }
            else if(loss=="custom"){
                CFValue <- lossFunction(actual=y,fitted=yFitted,B=B,xreg=matrixXreg);
            }
        }

        if(is.nan(CFValue) || is.na(CFValue) || is.infinite(CFValue)){
            CFValue[] <- 1E+300;
        }

        return(CFValue);
    }

    #### Define the rest of parameters ####
    ellipsis <- list(...);
    # ellipsis <- match.call(expand.dots = FALSE)$`...`;

    # Fisher Information
    if(is.null(ellipsis$FI)){
        FI <- FALSE;
    }
    else{
        FI <- ellipsis$FI;
    }

    # Starting values for the optimiser
    if(is.null(ellipsis$B)){
        B <- NULL;
    }
    else{
        B <- ellipsis$B;
    }
    # Parameters for the nloptr from the ellipsis
    if(is.null(ellipsis$xtol_rel)){
        xtol_rel <- 1E-6;
    }
    else{
        xtol_rel <- ellipsis$xtol_rel;
    }
    if(is.null(ellipsis$algorithm)){
        # if(recursiveModel){
        # algorithm <- "NLOPT_LN_BOBYQA";
        # }
        # else{
        algorithm <- "NLOPT_LN_SBPLX";
        # }
    }
    else{
        algorithm <- ellipsis$algorithm;
    }
    if(is.null(ellipsis$maxtime)){
        maxtime <- -1;
    }
    else{
        maxtime <- ellipsis$maxtime;
    }
    if(is.null(ellipsis$xtol_abs)){
        xtol_abs <- 1E-8;
    }
    else{
        xtol_abs <- ellipsis$xtol_abs;
    }
    if(is.null(ellipsis$ftol_rel)){
        ftol_rel <- 1E-4;
    }
    else{
        ftol_rel <- ellipsis$ftol_rel;
    }
    if(is.null(ellipsis$ftol_abs)){
        ftol_abs <- 0;
    }
    else{
        ftol_abs <- ellipsis$ftol_abs;
    }
    if(is.null(ellipsis$print_level)){
        print_level <- 0;
    }
    else{
        print_level <- ellipsis$print_level;
    }
    if(is.null(ellipsis$stepSize)){
        stepSize <- .Machine$double.eps^(1/4);
    }
    else{
        stepSize <- ellipsis$stepSize;
    }

    #### Form the necessary matrices -- needs to be rewritten for complex variables ####
    # Call similar to lm in order to form appropriate data.frame
    mf <- match.call(expand.dots = FALSE);
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L);
    mf <- mf[c(1L, m)];
    mf$drop.unused.levels <- TRUE;
    mf[[1L]] <- quote(stats::model.frame);

    # If data is provided explicitly, check it
    if(exists("data",inherits=FALSE,mode="numeric") || exists("data",inherits=FALSE,mode="complex") ||
       exists("data",inherits=FALSE,mode="list")){
        if(!is.data.frame(data)){
            data <- as.data.frame(data);
        }
        else{
            dataOrders <- unlist(lapply(data,is.ordered));
            # If there is an ordered factor, remove the bloody ordering!
            if(any(dataOrders)){
                data[dataOrders] <- lapply(data[dataOrders],function(x) factor(x, levels=levels(x), ordered=FALSE));
            }
        }
        mf$data <- data;
        rm(data);

        # If there are NaN values, remove the respective observations
        if(any(sapply(mf$data,is.nan))){
            warning("There are NaN values in the data. This might cause problems. Removing these observations.", call.=FALSE);
            NonNaNValues <- !apply(sapply(mf$data,is.nan),1,any);
            # If subset was not provided, change it
            if(is.null(mf$subset)){
                mf$subset <- NonNaNValues
            }
            else{
                mf$subset <- NonNaNValues & mf$subset;
            }
            dataContainsNaNs <- TRUE;
        }
        else{
            dataContainsNaNs <- FALSE;
        }

        # If there are spaces in names, give a warning
        if(any(grepl("[^A-Za-z0-9,;._-]", all.vars(formula))) ||
           # If the names only contain numbers
           any(suppressWarnings(!is.na(as.numeric(all.vars(formula)))))){
            warning("The names of your variables contain special characters ",
                    "(such as numbers, spaces, comas, brackets etc). clm() might not work properly. ",
                    "It is recommended to use `make.names()` function to fix the names of variables.",
                    call.=FALSE);
            formula <- as.formula(paste0(gsub(paste0("`",all.vars(formula)[1],"`"),
                                              make.names(all.vars(formula)[1]),
                                              all.vars(formula)[1]),"~",
                                         paste0(mapply(gsub, paste0("`",all.vars(formula)[-1],"`"),
                                                       make.names(all.vars(formula)[-1]),
                                                       labels(terms(formula))),
                                                collapse="+")));
            mf$formula <- formula;
        }
        # Fix names of variables. Switch this off to avoid conflicts between formula and data, when numbers are used
        colnames(mf$data) <- make.names(colnames(mf$data), unique=TRUE);

        # If the data is a matrix / data.frame, get the nrows
        if(!is.null(dim(mf$data))){
            obsAll <- nrow(mf$data);
        }
        else{
            obsAll <- length(mf$data);
        }

        # If the user asked for trend, but it's not in the data, add it
        if(any(all.vars(formula)=="trend") && all(colnames(mf$data)!="trend")){
            mf$data <- cbind(mf$data,trend=c(1:obsAll));
        }
    }
    else{
        dataContainsNaNs <- FALSE;
    }

    responseName <- all.vars(formula)[1];

    dataWork <- eval(mf, parent.frame());
    #### Temporary solution for model.matrix() ####
    complexVariables <- apply(dataWork, 2, is.complex);
    complexVariablesNames <- names(complexVariables)[complexVariables];
    responseIsComplex <- any(responseName==complexVariablesNames);
    complexVariablesNames <- complexVariablesNames[complexVariablesNames!=responseName];
    # Save the original data to come back to it
    originalData <- dataWork;
    dataWork[,] <- sapply(dataWork[,complexVariables], Re);

    dataTerms <- terms(dataWork);
    cl$formula <- formula(dataTerms);

    interceptIsNeeded <- attr(dataTerms,"intercept")!=0;
    # Create a model from the provided stuff. This way we can work with factors
    dataWork <- model.matrix(dataWork,data=dataWork);
    dataWork[,complexVariablesNames] <- as.matrix(originalData[,complexVariablesNames]);
    if(responseIsComplex){
        y <- originalData[,responseName]
    }
    obsInsample <- nrow(dataWork);

    # matrixXreg should not contain 1 for the further checks
    if(interceptIsNeeded){
        variablesNames <- colnames(dataWork)[-1];
        matrixXreg <- as.matrix(dataWork[,-1,drop=FALSE]);
        # Include response to the data
        # dataWork <- cbind(y,dataWork[,-1,drop=FALSE]);
    }
    else{
        variablesNames <- colnames(dataWork);
        matrixXreg <- dataWork;
        # Include response to the data
        # dataWork <- cbind(y,dataWork);
        warning("You have asked not to include intercept in the model. We will try to fit the model, ",
                "but this is a very naughty thing to do, and we cannot guarantee that it will work...", call.=FALSE);
    }
    # colnames(dataWork) <- c(responseName, variablesNames);
    rm(dataWork, originalData);

    nVariables <- length(variablesNames);
    colnames(matrixXreg) <- variablesNames;

    # Record the subset used in the model
    if(is.null(mf$subset)){
        subset <- rep(TRUE, obsInsample);
    }
    else{
        if(dataContainsNaNs){
            subset <- mf$subset[NonNaNValues];
        }
        else{
            subset <- mf$subset;
        }
    }

    mu <- vector("numeric", obsInsample);
    yFitted <- vector("numeric", obsInsample);
    errors <- vector("numeric", obsInsample);
    ot <- vector("logical", obsInsample);

    if(!fast){
        #### Checks of the exogenous variables ####
        # Remove the data for which sd=0
        noVariability <- vector("logical",nVariables);
        noVariability[] <- apply((matrixXreg==matrix(matrixXreg[1,],obsInsample,nVariables,byrow=TRUE)),2,all);
        if(any(noVariability)){
            if(all(noVariability)){
                warning("None of exogenous variables has variability. Fitting the straight line.",
                        call.=FALSE);
                matrixXreg <- matrix(1,obsInsample,1);
                nVariables <- 1;
                variablesNames <- "(Intercept)";
            }
            else{
                warning("Some exogenous variables did not have any variability. We dropped them out.",
                        call.=FALSE);
                matrixXreg <- matrixXreg[,!noVariability,drop=FALSE];
                nVariables <- ncol(matrixXreg);
                variablesNames <- variablesNames[!noVariability];
            }
        }

        #### Not clear yet how to check multicollinearity for complex models ####
        # Check the multicollinearity.
        # corThreshold <- 0.999;
        # if(nVariables>1){
        #     # Check perfectly correlated cases
        #     corMatrix <- cor(matrixXreg,use="pairwise.complete.obs");
        #     corHigh <- upper.tri(corMatrix) & abs(corMatrix)>=corThreshold;
        #     if(any(corHigh)){
        #         removexreg <- unique(which(corHigh,arr.ind=TRUE)[,1]);
        #         matrixXreg <- matrixXreg[,-removexreg,drop=FALSE];
        #         nVariables <- ncol(matrixXreg);
        #         variablesNames <- colnames(matrixXreg);
        #         if(!occurrenceModel && !CDF){
        #             warning("Some exogenous variables were perfectly correlated. We've dropped them out.",
        #                     call.=FALSE);
        #         }
        #     }
        # }
        #
        # # Do these checks only when intercept is needed. Otherwise in case of dummies this might cause chaos
        # if(nVariables>1 & interceptIsNeeded){
        #     # Check dummy variables trap
        #     detHigh <- suppressWarnings(determination(matrixXreg))>=corThreshold;
        #     if(any(detHigh)){
        #         while(any(detHigh)){
        #             removexreg <- which(detHigh>=corThreshold)[1];
        #             matrixXreg <- matrixXreg[,-removexreg,drop=FALSE];
        #             nVariables <- ncol(matrixXreg);
        #             variablesNames <- colnames(matrixXreg);
        #
        #             detHigh <- suppressWarnings(determination(matrixXreg))>=corThreshold;
        #         }
        #         if(!occurrenceModel){
        #             warning("Some combinations of exogenous variables were perfectly correlated. We've dropped them out.",
        #                     call.=FALSE);
        #         }
        #     }
        # }

        #### Finish forming the matrix of exogenous variables ####
        # Remove the redundant dummies, if there are any
        varsToLeave <- apply(matrixXreg,2,cvar)!=0;
        matrixXreg <- matrixXreg[,varsToLeave,drop=FALSE];
        variablesNames <- variablesNames[varsToLeave];
        nVariables <- length(variablesNames);
    }

    # Reintroduce intercept
    if(interceptIsNeeded){
        matrixXreg <- cbind(1,matrixXreg);
        variablesNames <- c("(Intercept)",variablesNames);
        colnames(matrixXreg) <- variablesNames;

        #### Not clear yet how to check multicollinearity for complex models ####
        # Check, if redundant dummies are left. Remove the first if this is the case
        # Don't do the check for LASSO / RIDGE
        # if(is.null(parameters) && !fast){
        #     determValues <- suppressWarnings(determination(matrixXreg[, -1, drop=FALSE]));
        #     determValues[is.nan(determValues)] <- 0;
        #     if(any(determValues==1)){
        #         matrixXreg <- matrixXreg[,-(which(determValues==1)[1]+1),drop=FALSE];
        #         variablesNames <- colnames(matrixXreg);
        #     }
        # }
        nVariables <- length(variablesNames);
    }
    variablesNamesAll <- variablesNames;
    # The number of exogenous variables
    nVariablesExo <- nVariables;

    #### Estimate parameters of the model ####
    if(is.null(parameters)){
        if(loss=="CLS"){
            B <- as.vector(invert(t(matrixXreg) %*% matrixXreg) %*% t(matrixXreg) %*% y);
            CFValue <- CF(B, loss, y, matrixXreg);
        }
        else if(loss=="OLS"){
            B <- as.vector(invert(t(Conj(matrixXreg)) %*% matrixXreg) %*% t(Conj(matrixXreg)) %*% y);
            CFValue <- CF(B, loss, y, matrixXreg);
        }
        else{
            # Set bounds for B as NULL. Then amend if needed
            BLower <- NULL;
            BUpper <- NULL;
            if(is.null(B)){
                B <- as.vector(invert(t(Conj(matrixXreg)) %*% matrixXreg) %*% t(Conj(matrixXreg)) %*% y);
                B <- c(Re(B),Im(B));
            }
            nVariables <- length(B);

            BLower <- rep(-Inf,nVariables);
            BUpper <- rep(Inf,nVariables);

            print_level_hidden <- print_level;
            if(print_level==41){
                print_level[] <- 0;
            }

            #### Define what to do with the maxeval ####
            if(is.null(ellipsis$maxeval)){
                maxeval <- nVariables * 40;
            }
            else{
                maxeval <- ellipsis$maxeval;
            }

            # Although this is not needed in case of distribution="dnorm", we do that in a way, for the code consistency purposes
            res <- nloptr(B, CF,
                          opts=list(algorithm=algorithm, xtol_rel=xtol_rel, maxeval=maxeval, print_level=print_level,
                                    maxtime=maxtime, xtol_abs=xtol_abs, ftol_rel=ftol_rel, ftol_abs=ftol_abs),
                          # lb=BLower, ub=BUpper,
                          loss=loss, y=y, matrixXreg=matrixXreg);
            B[] <- res$solution;
            B <- complex(real=B[1:(nVariables/2)],imaginary=B[(nVariables/2+1):nVariables]);
            nVariables <- length(B);
            CFValue <- res$objective;

            if(print_level_hidden>0){
                print(res);
            }
        }
    }
    # If the parameters are provided
    else{
        B <- parameters;
        nVariables <- length(B);
        if(!is.null(names(B))){
            variablesNames <- names(B);
        }
        else{
            names(B) <- variablesNames;
            names(parameters) <- variablesNames;
        }
        variablesNamesAll <- colnames(matrixXreg);
        CFValue <- CF(B, loss, y, matrixXreg);
    }

    fitterReturn <- fitter(B, y, matrixXreg);
    mu[] <- fitterReturn$mu;
    scale <- fitterReturn$scale;

    #### Produce Fisher Information ####
    # if(FI){
    #     # Only vcov is needed, no point in redoing the occurrenceModel
    #     FI <- hessian(CF, B, h=stepSize,
    #                   distribution=distribution, loss=loss, y=y, matrixXreg=matrixXreg,
    #                   recursiveModel=recursiveModel, denominator=denominator);
    #
    #     if(any(is.nan(FI))){
    #         warning("Something went wrong and we failed to produce the covariance matrix of the parameters.\n",
    #                 "Obviously, it's not our fault. Probably Russians have hacked your computer...\n",
    #                 "Try a different distribution maybe?", call.=FALSE);
    #         FI <- diag(1e+100,nVariables);
    #     }
    #     dimnames(FI) <- list(variablesNames,variablesNames);
    # }

    # Give names to additional parameters
    if(is.null(parameters)){
        parameters <- B;
        names(parameters) <- variablesNames;
        names(B) <- variablesNames;
    }

    ### Fitted values in the scale of the original variable
    yFitted[] <- extractorFitted(mu, scale);

    ### Error term in the transformed scale
    errors[] <- extractorResiduals(mu, yFitted);

    nParam <- nVariables + (loss=="likelihood")*3/2;

    if(interceptIsNeeded){
        # This shit is needed, because R has habit of converting everything into vectors...
        dataWork <- cbind(y,matrixXreg[,-1,drop=FALSE]);
        variablesUsed <- variablesNamesAll[variablesNamesAll!="(Intercept)"];
    }
    else{
        dataWork <- cbind(y,matrixXreg);
        variablesUsed <- variablesNamesAll;
    }
    colnames(dataWork) <- c(responseName, variablesUsed);

    # Return LogLik, depending on the used loss
    if(loss=="likelihood"){
        logLik <- -CFValue;
    }
    # else if(loss=="CLS"){
    #     logLik <- -CF(B, loss="likelihood", y, matrixXreg);
    # }
    else{
        logLik <- NA;
    }

    finalModel <- structure(list(coefficients=parameters, FI=FI, fitted=yFitted, residuals=as.vector(errors),
                                 mu=mu, scale=scale, logLik=logLik,
                                 loss=loss, lossFunction=lossFunction, lossValue=CFValue,
                                 df.residual=obsInsample-nParam, df=nParam, call=cl, rank=nParam,
                                 data=dataWork, terms=dataTerms,
                                 subset=subset, other=ellipsis, B=B,
                                 timeElapsed=Sys.time()-startTime),
                            class=c("clm","greybox"));

    return(finalModel);
}

#' @importFrom greybox actuals
#' @export
actuals.clm <- function(object, all=TRUE, ...){
    return(object$data[,1]);
}

#' @importFrom stats nobs
#' @export
nobs.clm <- function(object, all=TRUE, ...){
    return(length(actuals(object, ...)));
}

#' @importFrom greybox nparam
#' @export
nparam.clm <- function(object, all=TRUE, ...){
    return(object$rank);
    # x <- c(Re(coef(object)),Im(coef(object)));
    # x <- x[x!=0];
    # # Divide by two, because we calculate df per series
    # return(length(x)/2);
}

#' @importFrom stats logLik
#' @export
logLik.clm <- function(object, ...){
    return(structure(object$logLik,nobs=nobs(object),df=nparam(object),class="logLik"));
}

#' @rdname clm
#' @param object Model estimated using \code{clm()} function.
#' @param type Type of sigma to return. This is calculated based on the residuals
#' of the estimated model and can be \code{"direct"}, based on the direct variance,
#' \code{"conjugate"}, based on the conjugate variance and \code{"matrix"}, returning
#' covariance matrix for the complex error. If \code{NULL} then will return value based
#' on the loss used in the estimation: OLS -> "conjugate", CLS -> "direct", likelihood ->
#' "matrix".
#' @param ... Other parameters passed to internal functions.
#' @importFrom stats sigma
#' @export
sigma.clm <- function(object, type=NULL, ...){

    # See the type
    if(!is.null(type)){
        type <- match.arg(type, c("direct","conjugate","matrix"));
    }
    else{
        # Default values
        type <- switch(object$loss,
                       "CLS"="direct",
                       "OLS"="conjugate",
                       "likelihood"="matrix")
    }

    errors <- resid(object);
    sigmaValue <- switch(type,
                         "direct"=sqrt(sum(errors^2, ...)/(nobs(object) - nparam(object))),
                         "conjugate"=sqrt(sum(errors * Conj(errors), ...)/(nobs(object) - nparam(object))),
                         covar(errors, df=nobs(object)-nparam(object)));

    if(type=="matrix"){
        colnames(sigmaValue) <- rownames(sigmaValue) <- c("e_r", "e_i");
    }

    return(sigmaValue);
}

#' @rdname clm
#' @importFrom stats vcov
#' @export
vcov.clm <- function(object, type=NULL, ...){

    # See the type
    if(!is.null(type)){
        type <- match.arg(type, c("direct","conjugate","matrix"));
    }
    else{
        # Default values
        type <- switch(object$loss,
                       "CLS"="direct",
                       "OLS"="conjugate",
                       "likelihood"="matrix")
    }

    nVariables <- length(coef(object));
    variablesNames <- names(coef(object));
    interceptIsNeeded <- any(variablesNames=="(Intercept)");
    ellipsis <- list(...);

    matrixXreg <- object$data;
    if(interceptIsNeeded){
        matrixXreg[,1] <- 1;
        colnames(matrixXreg)[1] <- "(Intercept)";
    }
    else{
        matrixXreg <- matrixXreg[,-1,drop=FALSE];
    }

    sigmaValue <- sigma(object, type);

    if(type=="conjugate"){
        # Get variance
        sigmaValue[] <- sigmaValue^2;
        if(object$loss=="CLS"){
            # Conjugate matrix
            matrixXregConj <- Conj(matrixXreg);
            # Transposed original
            matrixXregTrans <- t(matrixXreg);
            # Conjugate transposed matrix
            matrixXregConjTrans <- t(Conj(matrixXreg));

            # Calculate the (X'X)^{-1} X'X~ (X~' X~)^{-1'}
            vcov <-
                invert(matrixXregTrans %*% matrixXreg) %*%
                (matrixXregTrans %*% matrixXregConj) %*%
                t(invert(t(matrixXregConj) %*% matrixXregConj));
        }
        # OLS and likelihood
        else {
            # Conjugate matrix
            matrixXregConj <- Conj(matrixXreg);
            # Transposed conjugate original
            matrixXregConjTrans <- t(Conj(matrixXreg));
            # Transposed matrix
            matrixXregTrans <- t(matrixXreg);

            # Calculate the (X'X)^{-1} X'X~ (X~' X~)^{-1'}
            vcov <-
                invert(matrixXregConjTrans %*% matrixXreg) %*%
                (matrixXregConjTrans %*% matrixXreg) %*%
                t(Conj(invert(matrixXregTrans %*% matrixXregConj)));
        }
    }
    else if(type=="direct"){
        # Get variance
        sigmaValue[] <- sigmaValue^2;
        if(object$loss=="CLS"){
            # Simple transposition of the matrix
            matrixXregTrans <- t(matrixXreg);
            matrixXreg <- matrixXreg;
        }
        # OLS and likelihood
        else {
            # Simple transposition of the matrix
            matrixXregTrans <- t(Conj(matrixXreg));
        }

        # Calculate the (X' X)^{-1}
        vcov <- invert(matrixXregTrans %*% matrixXreg);
    }
    else{
        # Transform the complex matrix to be a matrix
        matrixXreg <- complex2mat(matrixXreg);
        # Conjugate transposition
        matrixXregTrans <- t(matrixXreg);

        # Calculate the (X'X)^{-1}
        # Re() is needed to drop the 0i
        vcov <- Re(invert(matrixXregTrans %*% matrixXreg));
    }
    ndimVcov <- ncol(vcov);

    if(any(type==c("direct","conjugate"))){
        vcov[] <- vcov * sigmaValue;
        rownames(vcov) <- colnames(vcov) <- names(coef(object));
    }
    # Likelihood estimate
    else{
        for(i in 1:(ndimVcov/2)){
            for(j in 1:(ndimVcov/2)){
                vcov[(1:2)+(i-1)*2,(1:2)+(j-1)*2] <- vcov[(1:2)+(i-1)*2,(1:2)+(j-1)*2] * sigmaValue;
            }
        }
        rownames(vcov) <- colnames(vcov) <- names(complex2mat(coef(object))[,1]);
    }

    return(vcov);
}

#' @importFrom stats resid qt
#' @export
confint.clm <- function(object, parm, level = 0.95, ...){

    confintNames <- c(paste0((1-level)/2*100,"%"),
                      paste0((1+level)/2*100,"%"));

    # Extract coefficients
    parameters <- complex2mat(coef(object))[,1];
    parametersNames <- names(parameters);
    parametersLength <- length(parametersNames);

    if(object$loss=="likelihood"){
        # Get covariance matrix
        parametersSE <- sqrt(diag(vcov(object, type="matrix")));
    }
    else{
        parametersSEDir <- Re(diag(vcov(object, type="direct")));
        parametersSEConj <- Re(diag(vcov(object, type="conjugate")));
        parametersSE <- vector("numeric", parametersLength);
        # Real values
        parametersSE[1:(parametersLength/2)*2-1] <- sqrt((parametersSEConj + parametersSEDir)/2);
        # Imaginary values
        parametersSE[1:(parametersLength/2)*2] <- sqrt((parametersSEConj - parametersSEDir)/2);
    }
    varLength <- length(parametersSE);

    # Define quantiles using Student distribution
    paramQuantiles <- qt((1+level)/2,df=object$df.residual);

    # We can use normal distribution, because of the asymptotics of MLE
    confintValues <- cbind(parameters-paramQuantiles*parametersSE,
                           parameters+paramQuantiles*parametersSE);
    colnames(confintValues) <- confintNames;

    # Return S.E. as well, so not to repeat the thing twice...
    confintValues <- cbind(parametersSE, confintValues);
    # Give the name to the first column
    colnames(confintValues)[1] <- "S.E.";
    rownames(confintValues) <- parametersNames;

    return(confintValues);
}

#' @importFrom stats coef confint
#' @importFrom greybox nparam
#' @rdname clm
#' @param object Object of class "clm" estimated via \code{clm()} function.
#' @param level What confidence level to use for the parameters of the model.
#' @export
summary.clm <- function(object, level=0.95, ...){
    bootstrap <- FALSE;
    errors <- residuals(object);
    obs <- nobs(object, all=TRUE);

    # Collect parameters and their standard errors
    parametersConfint <- confint(object, level=level, bootstrap=bootstrap, ...);
    parameters <- complex2mat(coef(object))[,1];
    parametersLength <- length(parameters);
    parametersTable <- cbind(parameters,parametersConfint);
    parametersTableColnames <- c("Estimate","Std. Error",
                                 paste0("Lower ",(1-level)/2*100,"%"),
                                 paste0("Upper ",(1+level)/2*100,"%"));
    parametersTableRownames <- names(parameters)
    rownames(parametersTable) <- parametersTableRownames;
    colnames(parametersTable) <- parametersTableColnames;

    ourReturn <- list(coefficients=parametersTable);

    # Mark those that are significant on the selected level
    ourReturn$significance <- !(parametersTable[,3]<=0 & parametersTable[,4]>=0);

    # If there is a likelihood, then produce ICs
    if(!is.na(logLik(object))){
        ICs <- c(AIC(object),AICc(object),BIC(object),BICc(object));
        names(ICs) <- c("AIC","AICc","BIC","BICc");
        ourReturn$ICs <- ICs;
    }
    ourReturn$loss <- object$loss;
    ourReturn$other <- object$other;
    ourReturn$responseName <- formula(object)[[2]];

    # Table with degrees of freedom
    dfTable <- c(obs,nparam(object),obs-nparam(object));
    names(dfTable) <- c("n","k","df");

    ourReturn$r.squared <- 1 - sum(errors^2) / sum((actuals(object)-mean(actuals(object)))^2);
    ourReturn$adj.r.squared <- 1 - (1 - ourReturn$r.squared) * (obs - 1) / (dfTable[3]);

    ourReturn$dfTable <- dfTable;
    ourReturn$s2 <- sigma(object, type="matrix");

    ourReturn <- structure(ourReturn,class=c("summary.clm","summary.greybox"));
    return(ourReturn);
}

#' @export
print.summary.clm <- function(x, ...){
    ellipsis <- list(...);
    if(!any(names(ellipsis)=="digits")){
        digits <- 4;
    }
    else{
        digits <- ellipsis$digits;
    }

    cat("Complex Linear Regression estimated via clm()\n");
    cat(paste0("Response variable: ", paste0(x$responseName,collapse="")));
    cat(paste0("\nLoss function used in estimation: ",x$loss));

    cat("\nCoefficients:\n");
    stars <- setNames(vector("character",length(x$significance)),
                      names(x$significance));
    stars[x$significance] <- "*";
    print(data.frame(round(x$coefficients,digits),stars,
                     check.names=FALSE,fix.empty.names=FALSE));

    cat("\nError covariance matrix:\n"); print(round(x$s2,digits));
    cat("\nSample size: "); cat(x$dfTable[1]);
    cat("\nNumber of estimated parameters: "); cat(x$dfTable[2]);
    cat("\nNumber of degrees of freedom: "); cat(x$dfTable[3]);
    if(!is.null(x$ICs)){
        cat("\nInformation criteria:\n");
        print(round(x$ICs,digits));
    }
    cat("\n");
}

#' @export
plot.clm <- function(x, which=c(1,2,4,6), ...){
    # Amend the object to make it work with legion plots
    x$Sigma <- sigma(x, type="matrix");
    x$data <- complex2vec(actuals(x));
    x$fitted <- complex2vec(fitted(x));
    x$residuals <- complex2vec(residuals(x));
    x$forecast <- matrix(NA,1,2);
    x$model <- "VES(ANN)";
    class(x) <- "legion";

    if(any(which>=12)){
        which <- which[which<12];
    }

    plot(x, which=which, ...);
}

#' @export
predict.clm <- function(object, newdata=NULL, interval=c("none", "confidence", "prediction"),
                            level=0.95, side=c("both","upper","lower"), ...){
    interval <- match.arg(interval);
    side <- match.arg(side);

    parameters <- coef(object);
    parametersNames <- names(parameters);

    nLevels <- length(level);
    levelLow <- levelUp <- vector("numeric",nLevels);
    if(side=="upper"){
        levelLow[] <- 0;
        levelUp[] <- level;
    }
    else if(side=="lower"){
        levelLow[] <- 1-level;
        levelUp[] <- 1;
    }
    else{
        levelLow[] <- (1-level) / 2;
        levelUp[] <- (1+level) / 2;
    }
    paramQuantiles <- qt(c(levelLow, levelUp),df=object$df.residual);

    if(is.null(newdata)){
        matrixOfxreg <- object$data;
        newdataProvided <- FALSE;
        # The first column is the response variable. Either substitute it by ones or remove it.
        if(any(parametersNames=="(Intercept)")){
            matrixOfxreg[,1] <- 1;
        }
        else{
            matrixOfxreg <- matrixOfxreg[,-1,drop=FALSE];
        }
    }
    else{
        newdataProvided <- TRUE;

        if(!is.data.frame(newdata)){
            if(is.vector(newdata)){
                newdataNames <- names(newdata);
                newdata <- matrix(newdata, nrow=1, dimnames=list(NULL, newdataNames));
            }
            newdata <- as.data.frame(newdata);
        }
        else{
            dataOrders <- unlist(lapply(newdata,is.ordered));
            # If there is an ordered factor, remove the bloody ordering!
            if(any(dataOrders)){
                newdata[dataOrders] <- lapply(newdata[dataOrders],function(x) factor(x, levels=levels(x), ordered=FALSE));
            }
        }

        # The gsub is needed in order to remove accidental special characters
        colnames(newdata) <- make.names(colnames(newdata), unique=TRUE);

        # Extract the formula and get rid of the response variable
        testFormula <- formula(object);

        # If the user asked for trend, but it's not in the data, add it
        if(any(all.vars(testFormula)=="trend") && all(colnames(newdata)!="trend")){
            newdata <- cbind(newdata,trend=nobs(object)+c(1:nrow(newdata)));
        }

        testFormula[[2]] <- NULL;
        # Expand the data frame
        newdataExpanded <- model.frame(testFormula, newdata);
        interceptIsNeeded <- attr(terms(newdataExpanded),"intercept")!=0;

        #### Temporary solution for model.matrix() ####
        responseName <- all.vars(formula(object))[[1]];
        complexVariables <- apply(newdataExpanded, 2, is.complex);
        complexVariablesNames <- names(complexVariables)[complexVariables];
        complexVariablesNames <- complexVariablesNames[complexVariablesNames!=responseName];
        # Save the original data to come back to it
        originalData <- newdataExpanded;
        newdataExpanded[,] <- sapply(newdataExpanded[,complexVariables], Re);

        # Create a model from the provided stuff. This way we can work with factors
        matrixOfxreg <- model.matrix(newdataExpanded,data=newdataExpanded);
        matrixOfxreg[,complexVariablesNames] <- originalData[,complexVariablesNames];
        matrixOfxreg <- matrixOfxreg[,parametersNames,drop=FALSE];
    }

    h <- nrow(matrixOfxreg);

    if(!is.matrix(matrixOfxreg)){
        matrixOfxreg <- as.matrix(matrixOfxreg);
        h <- nrow(matrixOfxreg);
    }

    if(h==1){
        matrixOfxreg <- matrix(matrixOfxreg, nrow=1);
    }

    ourForecast <- as.vector(matrixOfxreg %*% parameters);
    vectorOfVariances <- NULL;

    if(interval!="none"){
        matrixOfxreg <- complex2mat(matrixOfxreg);
        ourVcov <- vcov(object, ...);
        # abs is needed for some cases, when the likelihood was not fully optimised
        vectorOfVariances <- complex2vec(diag(mat2complex(matrixOfxreg %*% ourVcov %*% t(matrixOfxreg))));

        #### Confidence works somehow. Prediction doesn't.

        yUpper <- yLower <- matrix(NA, h, nLevels);
        if(interval=="confidence"){
            for(i in 1:nLevels){
                yLower[,i] <- ourForecast + paramQuantiles[i] * vec2complex(sqrt(abs(vectorOfVariances)));
                yUpper[,i] <- ourForecast + paramQuantiles[i+nLevels] * vec2complex(sqrt(abs(vectorOfVariances)));
            }
        }
        else if(interval=="prediction"){
            sigmaValues <- sigma(object);
            vectorOfVariances[] <- vectorOfVariances + matrix(diag(sigmaValues),h,2,byrow=TRUE);
            for(i in 1:nLevels){
                yLower[,i] <- ourForecast + paramQuantiles[i] * vec2complex(sqrt(vectorOfVariances));
                yUpper[,i] <- ourForecast + paramQuantiles[i+nLevels] * vec2complex(sqrt(vectorOfVariances));
            }
        }

        colnames(yLower) <- switch(side,
                                   "both"=,
                                   "lower"=paste0("Lower bound (",levelLow*100,"%)"),
                                   "upper"=rep("Lower 0%",nLevels));

        colnames(yUpper) <- switch(side,
                                   "both"=,
                                   "upper"=paste0("Upper bound (",levelUp*100,"%)"),
                                   "lower"=rep("Upper 100%",nLevels));
    }
    else{
        yLower <- NULL;
        yUpper <- NULL;
    }

    ourModel <- list(model=object, mean=ourForecast, lower=yLower, upper=yUpper, level=c(levelLow, levelUp), newdata=newdata,
                     variances=vectorOfVariances, newdataProvided=newdataProvided);
    return(structure(ourModel,class=c("predict.clm","predict.greybox")));
}

#' @export
plot.predict.clm <- function(x, ...){
    # Amend the object to make it work with legion plots
    x$data <- complex2vec(actuals(x$model));
    x$fitted <- complex2vec(fitted(x$model));
    x$forecast <- complex2vec(x$mean);
    if(!is.null(x$lower)){
        x$PI <- complex2vec(cbind(x$lower,x$upper))[,c(1,3,2,4)];
    }
    x$model <- "CLM";
    class(x) <- "legion";

    plot(x, which=7, ...);
}
