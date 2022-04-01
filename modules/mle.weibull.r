mle.weibull <- function(data, data.censored=NA, param='auto', param.control=2, plots=FALSE, debug=FALSE) {
    
    ## weibull distribution
    ## MLE (Maximum Likelihood Estimate) fit to determine parameters
    ## LR (Likelihood Ratio) appraoch to find tolerance limit

    ## input: data  = vector of data
    ##        param = initial guess for fit parameters for: shape and scale
    ##                if type is also provided, it will not be used
    ##              = 'auto' (default) uses tolerance::exttol.int() for initial guess of parameters
    ##                and will try using list(shape=1, scale=1) if that fails

    ## based on approach found here:
    ## https://www.r-bloggers.com/2019/08/maximum-likelihood-estimation-from-scratch/
    ## https://personal.psu.edu/abs12/stat504/Lecture/lec3_4up.pdf
    x <- data
    scale <- NA
    shape <- NA

    ## make sure Weibull is possible
    if (min(x, na.rm=TRUE) < 0) {
        cat('####################################################\n')
        cat('               ERROR: MINIMUM x < 0                 \n')
        cat('####################################################\n')
        return(list(scale=scale, shape=shape))
    }
  
    if (is.data.frame(data.censored[1])) {
        ## censored data also provided
        xcen <- data.frame(x.low = data.censored[[1]], x.high = data.censored[[2]])
    } else {
        xcen <- NA
    }
    
    ## if (isTRUE(plots)) par(mfrow=c(1,2))
    out <- NULL

    ##-----------------------------------------------------------------------------
    ## let R figure out which Johnson distribution fits best
    tol_out_weib <- tolerance::exttol.int(x)
    ## params.compare <- as.data.frame(t(unlist(parms.extol)))
    shape   <- tol_out_weib$'shape.1'
    scale   <- tol_out_weib$'shape.2'
    param <- list(shape=shape, scale=scale)
    params.compare <- as.data.frame(param)
    params.compare$description <- 'tolerance::exttol.int(x)'

    ##-----------------------------------------------------------------------------
    ## set initial guess for MLE fit
    if (param[1] == 'auto') {
        if (is.na(parms$shape) | is.na(parms$scale)) {
            ## exttol failed to converge 
            param <- list(shape = 1,
                          scale = 1)
        }
    }
    temp <- as.data.frame(t(unlist(param)))
    temp$description <- 'Initial guess for MLE'
    params.compare <- fastmerge(params.compare, temp)
    
    ##-----------------------------------------------------------------------------
    ## determine best fit using nll
    param.fix <- function(param, param.control) {
        if (param.control == 1) {
            ## keep param positive so subsequent functions are defined
            param <- max(1E-15, param)
        } else if (param.control == 2) {
            ## keep param positive so subsequent functions are defined
            param <- abs(param)
            if (param == 0) param <- 1E-15
        }
        return(param)
    }
    nll <- function(data, data.censored=NA, param, param.control=0, debug=FALSE){
        ## calculate nll (negative log likelihhod) for distribution
        x      <- data
        xcen   <- data.censored
        shape  <- param[[1]]  
        scale  <- param[[2]]
        scale  <- param.fix(scale, param.control)
        pdf <- shape / scale^shape * x^(shape-1) * exp(-(x/scale)^shape)
        ## the above is equivalent
        ## pdf <- stats::dweibull(x, shape, scale)
        ## cdf <- pnorm(z)
        ## plot(z, cdf)
        if (is.data.frame(xcen)) {
            xcen$F.low  <- 1 - exp(- (xcen$x.low  / scale)^shape) # CDF; could use ExtDist::pWeibull(xcen$x.low, shape, scale)
            xcen$F.high <- 1 - exp(- (xcen$x.high / scale)^shape)
            ## if low CDF is NA, set to 0
            xcen$F.low[is.na(xcen$F.low)]   <- 0
            ## if high CDF is NA, set to 1
            xcen$F.high[is.na(xcen$F.high)] <- 1
            ## calculate probability for the censored interval
            xcen$probability <- xcen$F.high - xcen$F.low
            nll     <- -sum(log(pdf), log(xcen$probability))
        } else {
            nll     <- -sum(log(pdf))
        }
        if (isTRUE(debug)) cat('shape=', signif(shape,11),
                               'scale=', signif(scale,11),
                               'nll=', signif(nll,11), "\n")
        return(nll)
    }        
    ## print('Attempting MLE fit on regular parameters')
    tryCatch({
        out.bestfit <- optim(par     = param, 
                             fn      = nll, 
                             data    = x,
                             data.censored = xcen,
                             debug   = debug,
                             param.control = param.control,
                             control = list(trace=TRUE),
                             method  = "BFGS")
        nll.max.bestfit <- out.bestfit$value
        parms.mle <- as.list(out.bestfit$par)
    }, error = function(e) {
        ## what to do if error
        cat('WARNING: CONVERGENCE FAILURE IN mle.weibull()\n')
        parms.mle <- list(shape=NA, scale=NA)
    })
    ## the following is needed if use param.control because the
    ## parameter returned by optim() is the input pamameter to the function
    ## not the output parameter. So if the input parameter is negative,
    ## the optim() routine thinks that was OK as the answer even though
    ## it really only ever used abs(param) in the optimization.
    parms.mle$scale  <- param.fix(parms.mle$scale,  param.control)
        
    ## add MLE parameters to params.compare dataframe
    temp <- as.data.frame(parms.mle)
    temp$description <- 'MLE'
    params.compare <- fastmerge(params.compare, temp)

    ## remove factor levels from params.compare
    params.compare <- droplevels.all(params.compare)

    if (isTRUE(plots)) {
        out.hist <- hist(x, plot=FALSE)
        shape <- parms.mle$shape
        scale <- parms.mle$scale
        curve.points <- stats::dweibull(x, shape, scale)
        hist(x, freq=FALSE, ylim=range(out.hist$density, curve.points))
        curve(stats::dweibull(x, shape, scale), min(x), max(x), add=TRUE)
        qqplot_nwj(x, type='w', jfit=parms.mle)
    }

    return(list(parms=parms.mle, parms.compare=params.compare))
}



mle.weibull.test <- function() {
    source("F:\\Documents\\01_Dave's Stuff\\Programs\\GitHub_home\\R-setup\\setup.r")
    source('setup.r')
    
    fit.compare <- function(x, xcen, main=NULL) {
        ## plot histogram with no censored data and fit
        hist(x, freq=FALSE, border='black', main=main)
        out.fit0 <- mle.weibull(x, data.censored=NA, plots=FALSE)
        parms0 <- out.fit0$parms
        curve(stats::dweibull(x, parms0$shape, parms0$scale), min(x), max(x), col='black', add=TRUE)
        ## plot histogram with all data treated as known fit at numeric value
        xcen.avg <- rowMeans(xcen, na.rm=TRUE) # use the average for interval data
        x.all <- c(x, xcen.avg)
        hist(x.all, freq=FALSE, border='red', add=TRUE)
        out.fit1 <- mle.weibull(x.all, data.censored=NA, plots=FALSE)
        parms1 <- out.fit1$parms
        curve(stats::dweibull(x, parms1$shape, parms1$scale), min(x), max(x), col='red', add=TRUE)
        ## plot fit if treat xcen as censored
        out.fit2 <- mle.weibull(x, data.censored=xcen, plots=FALSE)
        parms2 <- out.fit2$parms
        curve(stats::dweibull(x, parms2$shape, parms2$scale), min(x), max(x), col='blue', type='p', add=TRUE)
        ## add legend
        legend('topright', 
               legend=c('known', 'all',   'censored'),
               col   =c('black', 'red', 'blue'),
               lty   =c( 1     ,  1   ,  NA),
               pch   =c( NA    ,  NA  ,  1))
    }

    plotspace(2,2)
    x <- iris$Sepal.Width
    xnum <- 10
    x.low  <- NA
    x.high <- 2.2
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit.compare(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    
    x.low  <- 3.7
    x.high <- NA
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit.compare(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    
    x.low  <- 2.2
    x.high <- 3.7
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit.compare(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    
    x.low  <- 2.2
    x.high <- NA
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit.compare(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
}

