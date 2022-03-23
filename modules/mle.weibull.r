mle.weibull <- function(data, data.censored=NA, param='auto', param.control=2, plots=FALSE, debug=FALSE) {
    
    ## weibull distribution
    ## MLE (Maximum Likelihood Estimate) fit to determine parameters
    ## LR (Likelihood Ratio) appraoch to find tolerance limit

    ## input: data  = vector of data
    ##        param = initial guess for fit parameters for: gamma, delta, xi, and lambda
    ##                if type is also provided, it will not be used
    ##              = 'auto' uses tolerance package to determine initial guess

    ## based on approach found here:
    ## https://www.r-bloggers.com/2019/08/maximum-likelihood-estimation-from-scratch/
    ## https://personal.psu.edu/abs12/stat504/Lecture/lec3_4up.pdf
    x <- data
  
    if (is.data.frame(data.censored[1])) {
        ## censored data also provided
        xcen <- data.frame(x.low = data.censored[[1]], x.high = data.censored[[2]])
    } else {
        xcen <- NA
    }
    
    ## if (isTRUE(plots)) par(mfrow=c(1,2))
    out <- NULL

    ##-----------------------------------------------------------------------------
    ## initial guesses for Johnson parameters
    if (param[1] == 'auto') {
        ## use R to guess fit parameters
        tol_out_weib <-  tolerance::exttol.int(x, alpha=alpha, P=P, side=1, dist="Weibull")
        shape   <- tol_out_weib$'shape.1'
        scale   <- tol_out_weib$'shape.2'
        param <- list(shape=shape, scale=scale)
        params.compare <- as.data.frame(param)
        params.compare$description <- 'tolerance::extol.int(x)'
    } else {
        params.compare <- as.data.frame(param)
        params.compare$description <- 'input parameters'
    }
    params.guess <- param
    
    ##-----------------------------------------------------------------------------
    ## determine best fit using nll
    param.fix <- function(param, param.control) {
        if (param.control == 1) {
            ## keep param positive so subsequent functions are defined
            param <- param(1E-15, param)
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
        shape  <- param[[1]]  
        scale  <- param[[2]]
        scale  <- param.fix(scale, param.control)
        pdf <- shape / scale^shape * x^(shape-1) * exp(-(x/scale)^shape)
        ## the above is equivalent
        ## pdf <- dWeibull(x, shape, scale)
        if (is.data.frame(xcen)) {
            xcen$F.low  <- 1 - exp(- (xcen$x.low  / scale)^shape) # CDF
            xcen$F.high <- 1 - exp(- (xcen$x.high / scale)^shape) # CDF
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
        if (isTRUE(debug)) cat('shape=', signif(shape,11), 'scale=', signif(scale,11), 'nll=', signif(nll,11), "\n")
        return(nll)
    }        
    ## print('Attempting MLE fit on regular parameters')
    out.bestfit <- stats::optim(par     = param, 
                                fn      = nll, 
                                data    = x,
                                data.censored = xcen,
                                debug   = debug,
                                control = list(trace=TRUE),
                                method  = "BFGS")
    nll.max.bestfit <- out.bestfit$value
    params.mle <- as.list(out.bestfit$par)
    params.mle$scale <- param.fix(params$scale, param.control)
    shape  <- params.mle$shape
    scale  <- params.mle$scale
    print(as.data.frame(params.mle))
    cat('\n')


    ## add MLE parameters to params.compare dataframe
    temp <- as.data.frame(params.mle)
    temp$description <- 'MLE'
    params.compare <- fastmerge(params.compare, temp)

    ## remove factor levels from params.compare
    params.compare <- droplevels.all(params.compare)

    if (isTRUE(plots)) {
        out.hist <- hist(x, plot=FALSE)
        curve.points <- stats::dweibull(x, shape, scale)
        hist(x, freq=FALSE, ylim=range(out.hist$density, curve.points))
        curve(stats::dweibull(x, shape, scale), min(x), max(x), add=TRUE)
        qqplot_nwj(x, type='w', jfit=params.mle)
    }

    
    return(list(parms=params.mle, parms.compare=params.compare))
}

