#' Complex Linear Model
#'
#' Function estimates complex variables model
#'
#' This is a function, similar to \link[stats]{lm}, but supporting several estimation
#' techniques for complex variables regression.
#'
#' @template author
#' @template keywords
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
#' \item \code{LS} - least squares method for complex variables;
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
#' @importFrom stats .lm.fit formula residuals sigma
#' @export clm
clm <- function(formula, data, subset, na.action,
                loss=c("LS","likelihood","MSE","MAE","HAM"),
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
        scale <- sqrt(sum((y-mu)^2)/obsInsample)

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
        fitterReturn <- fitter(B, y, matrixXreg);

        if(loss=="likelihood"){
            # The original log-likelilhood
            #### !!! Needs to be modified for complex case !!! ####
            # CFValue <- -sum("dnorm" = dnorm(y, mean=fitterReturn$mu, sd=fitterReturn$scale, log=TRUE));
            CFValue <- NA;
        }
        else if(loss=="LS"){
            CFValue <- sum((y - fitterReturn$mu)^2);
        }
        else{
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
        # LASSO / RIDGE need more accurate estimation
        if(any(loss==c("LASSO","RIDGE"))){
            ftol_rel <- 1e-8;
        }
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
    dataWork[,complexVariablesNames] <- originalData[,complexVariablesNames];
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
        varsToLeave <- apply(matrixXreg,2,pvar)!=0;
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
        if(loss=="LS"){
            B <- as.vector(invert(t(matrixXreg) %*% matrixXreg) %*% t(matrixXreg) %*% y);
            CFValue <- CF(B, loss, y, matrixXreg);
        }
        else{
            # Set bounds for B as NULL. Then amend if needed
            BLower <- NULL;
            BUpper <- NULL;
            if(is.null(B)){
                B <- .lm.fit(matrixXreg,y)$coefficients;
                BLower <- -Inf;
                BUpper <- Inf;
            }

            BLower <- rep(-Inf,length(B));
            BUpper <- rep(Inf,length(B));

            print_level_hidden <- print_level;
            if(print_level==41){
                print_level[] <- 0;
            }

            #### Define what to do with the maxeval ####
            if(is.null(ellipsis$maxeval)){
                maxeval <- length(B) * 40;
            }
            else{
                maxeval <- ellipsis$maxeval;
            }

            # Although this is not needed in case of distribution="dnorm", we do that in a way, for the code consistency purposes
            res <- nloptr(B, CF,
                          opts=list(algorithm=algorithm, xtol_rel=xtol_rel, maxeval=maxeval, print_level=print_level,
                                    maxtime=maxtime, xtol_abs=xtol_abs, ftol_rel=ftol_rel, ftol_abs=ftol_abs),
                          lb=BLower, ub=BUpper,
                          loss=loss, y=y, matrixXreg=matrixXreg);
            B[] <- res$solution;
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

    nParam <- nVariables + (loss=="likelihood")*1;

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
    # else if(loss=="LS"){
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

#' @importFrom stats logLik
#' @export
logLik.clm <- function(object, ...){
    return(structure(object$logLik,nobs=nobs(object),df=nparam(object),class="logLik"));
}

#' @importFrom stats vcov
#' @export
vcov.clm <- function(object, ...){
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

    vcov <- invert(t(matrixXreg) %*% matrixXreg) * sigma(object)^2;
    rownames(vcov) <- colnames(vcov) <- variablesNames;

    return(vcov);
}

#' @importFrom stats coef confint
#' @importFrom greybox nparam
#' @export
summary.clm <- function(object, level=0.95, bootstrap=FALSE, ...){
    errors <- residuals(object);
    obs <- nobs(object, all=TRUE);

    # Collect parameters and their standard errors
    parametersConfint <- confint(object, level=level, bootstrap=bootstrap, ...);
    parameters <- coef(object);
    parametersTable <- cbind(parameters,sqrt(diag(vcov(object))),parametersConfint);
    rownames(parametersTable) <- names(parameters);
    colnames(parametersTable) <- c("Estimate","Std. Error",
                                   paste0("Lower ",(1-level)/2*100,"%"),
                                   paste0("Upper ",(1+level)/2*100,"%"));
    ourReturn <- list(coefficients=parametersTable);
    # Mark those that are significant on the selected level
    ourReturn$significance <- !(abs(parametersTable[,3])<=0 & abs(parametersTable[,4])>=0);

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
    ourReturn$s2 <- sigma(object)^2;

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

    cat(paste0("Response variable: ", paste0(x$responseName,collapse="")));
    cat(paste0("\nLoss function used in estimation: ",x$loss));

    cat("\nCoefficients:\n");
    stars <- setNames(vector("character",length(x$significance)),
                      names(x$significance));
    stars[x$significance] <- "*";
    print(data.frame(round(x$coefficients,digits),stars,
                     check.names=FALSE,fix.empty.names=FALSE));

    cat("\nError standard deviation: "); cat(round(sqrt(x$s2),digits));
    cat("\nSample size: "); cat(x$dfTable[1]);
    cat("\nNumber of estimated parameters: "); cat(x$dfTable[2]);
    cat("\nNumber of degrees of freedom: "); cat(x$dfTable[3]);
    if(!is.null(x$ICs)){
        cat("\nInformation criteria:\n");
        print(round(x$ICs,digits));
    }
    cat("\n");
}
