mle.johnsonsu <- function(data, data.censored=NA, param='auto', param.control=2, plots=FALSE, debug=FALSE) {
    
    ## johnsonsu distribution
    ## MLE (Maximum Likelihood Estimate) fit to determine parameters
    ## LR (Likelihood Ratio) appraoch to find tolerance limit

    ## input: data  = vector of data
    ##        param = initial guess for fit parameters for: gamma, delta, xi, and lambda
    ##                if type is also provided, it will not be used
    ##              = 'auto' (default) uses ExtDist::eJohnsonSU() for initial guess of parameters,
    ##                will switch to SuppDists::JohnsonFit(x) if that fails, and will switch to
    ##                try using list(gamma=xx, delta=xx, xi=xx, lambda=xx) if that fails

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
    ## if (param[1] == 'SuppDists') {
    ##     ## let R figure out which Johnson distribution fits best
    ##     jparms <- SuppDists::JohnsonFit(x)
    ##     if ('type' != 'SU') {
    ##         cat('\n')
    ##         cat('################################################\n')
    ##         cat('                  WARNING:\n')
    ##         cat('SuppDists::JohnsonFit() did not return type="SU"\n')
    ##         cat('Try param="ExtDist" or user defined list.\n')
    ##         print(data.frame(t(unlist(jparms))))
    ##         cat('################################################\n')
    ##         cat('\n')
    ##     }
    ## } else if (param[1] == 'ExtDist') {

    ##-----------------------------------------------------------------------------
    ## let R figure out which Johnson distribution fits best
    jparms.SuppDists <- SuppDists::JohnsonFit(x)
    ## params.compare <- as.data.frame(t(unlist(jparms.SuppDists)))
    params.compare <- as.data.frame(jparms.SuppDists)
    params.compare$description <- 'SuppDists::JohnsonFit(x)'

    ##-----------------------------------------------------------------------------
    ## force the Johnson SU distribution
    jparms.out <- ExtDist::eJohnsonSU(x)
    jparms.ExtDist <- list(gamma   = jparms.out$gamma,
                           delta   = jparms.out$delta,
                           xi      = jparms.out$xi,
                           lambda  = jparms.out$lambda,
                           type    = 'SU')
    temp <- as.data.frame(t(unlist(jparms.ExtDist)))
    temp$description <- 'ExtDist::eJohnsonSU(x)'
    params.compare <- fastmerge(params.compare, temp)

    ##-----------------------------------------------------------------------------
    ## set initial guess for MLE fit
    if (param[1] == 'auto') {

        ## force the Johnson SU distribution
        jparms <- jparms.ExtDist

        if (jparms$gamma == -0.5 & jparms$delta == 2.0 & jparms$xi == -0.5 & jparms$lambda == 2.0) {
            ## ExtDist failed to converge so try SuppDists
            jparms <- jparms.SuppDists
            if (jparms$type != 'SU') {
                ## SuppDists::JohnsonFit(x) failed to converge or did not return type=SU
                ## so arbitrarily set jparms to something as a starting point
                jparms <- list(gamma   = 1,
                               delta   = 1,
                               xi      = 1,
                               lambda  = 1,
                               type    = 'SU')
            }
         }
        
    } else {
        ## use Johnson parameters specified in param
        ## needs to be in same list format as created by SuppDists::JohnsonFit
        jparms <- param
    }
    temp <- as.data.frame(t(unlist(jparms)))
    temp$description <- 'Initial guess for MLE'
    params.compare <- fastmerge(params.compare, temp)

    ## strip off type if that is provided in param list since not needed for MLE fit
    param <- jparms[1:4]   # list(gamma, delta, xi, lambda)
    
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
        gamma  <- param[[1]]  
        delta  <- param[[2]]
        xi     <- param[[3]]
        lambda <- param[[4]]
        delta  <- param.fix(delta,  param.control)
        lambda <- param.fix(lambda, param.control)
        ## PDF for Johnson SU
        z <- gamma + delta * asinh( (x-xi)/lambda )
        pdf <- delta / (lambda*sqrt(2*pi)) / sqrt(1 +((x-xi)/lambda)^2) * exp(-0.5*z^2)
        ## the above is equivalent
        ## pdf <- ExtDist::dJohnsonSU(x, parms=c(gamma, delta, xi, lambda))
        ## cdf <- pnorm(z)
        ## plot(z, cdf)
        if (is.data.frame(xcen)) {
            xcen$F.low  <- pnorm(xcen$x.low)
            xcen$F.high <- pnorm(xcen$x.high)
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
        if (isTRUE(debug)) cat('gamma=', signif(gamma,11),
                               'delta=', signif(delta,11),
                               'xi   =', signif(xi, 11),
                               'lambda=', signif(lambda,11),
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
        jparms.mle <- as.list(out.bestfit$par)
        jparms.mle$type <- 'SU'
    }, error = function(e) {
        ## what to do if error
        cat('WARNING: CONVERGENCE FAILURE IN mle.johnsonsu()\n')
        jparms.mle <- list(gamma=NA, delta=NA, xi=NA, lambda=NA, type=NA)
    })
    ## the following is needed if use param.control because the
    ## parameter returned by optim() is the input pamaeter to the function
    ## not the output parameter. So if the input lambda is negative,
    ## the optim() routine thinks that was OK as the answer even though
    ## it really only ever used abs(lambda) and abs(delta) in the optimization.
    jparms.mle$delta  <- param.fix(jparms.mle$delta,  param.control)
    jparms.mle$lambda <- param.fix(jparms.mle$lambda, param.control)

    ## add MLE parameters to params.compare dataframe
    temp <- as.data.frame(jparms.mle)
    temp$description <- 'MLE'
    params.compare <- fastmerge(params.compare, temp)

    ## remove factor levels from params.compare
    params.compare <- droplevels.all(params.compare)

    if (isTRUE(plots)) {
        out.hist <- hist(x, plot=FALSE)
        curve.points <- ExtDist::dJohnsonSU(x, params=jparms.mle)
        hist(x, freq=FALSE, ylim=range(out.hist$density, curve.points))
        curve(ExtDist::dJohnsonSU(x, params=jparms.mle), min(x), max(x), add=TRUE)
        qqplot_nwj(x, type='j', jfit=jparms.mle)
    }

    
    return(list(jparms=jparms.mle, jparms.compare=params.compare))
}



mle.johnsonsu.test <- function() {
    source("F:\\Documents\\01_Dave's Stuff\\Programs\\GitHub_home\\R-setup\\setup.r")
    x <- iris$Sepal.Width
    plotspace(3,2)
    xcen <- data.frame(x.low=c(NA, 2, 2.7), x.high=c(0.5,2.5,NA))
    out.fit1 <- mle.johnsonsu(x, data.censored=NA,   plots=TRUE)
    out.fit2 <- mle.johnsonsu(x, data.censored=xcen, plots=TRUE)
    xcen <- data.frame(x.low=c(NA, NA, NA), x.high=c(2.2, 2.3, 2.4))
    out.fit3 <- mle.johnsonsu(x, data.censored=xcen, plots=TRUE)
    print(out.fit1$jparms.compare)
    print(out.fit2$jparms.compare)
    print(out.fit3$jparms.compare)

    ## little impact seen from many data at low end? 
    ## makes sense if survival data but not for my purpose
    xcen <- data.frame(x.low=rep(NA,99), x.high=rep(2.2,99))
    out.fit1 <- mle.johnsonsu(x, data.censored=xcen, plots=TRUE)
    print(out.fit1$jparms.compare)

    ## little impact from many data with huge uncertainty in timing
    xcen <- data.frame(x.low=rep(NA,99), x.high=rep(3.7,99))
    out.fit2 <- mle.johnsonsu(x, data.censored=xcen, plots=TRUE)
    print(out.fit2$jparms.compare)

    ## little impact from many data at high end????????
    xcen <- data.frame(x.low=rep(3.7,99), x.high=rep(NA,99))
    out.fit3 <- mle.johnsonsu(x, data.censored=xcen, plots=TRUE)
    print(out.fit3$jparms.compare)
    ## huge impact if add those data as real data    
    x <- c(x, xcen[[1]])
    out.fit4 <- mle.johnsonsu(x, data.censored=NA, plots=TRUE)
    print(out.fit4$jparms.compare)
    ## what if add those data with small window?
    x <- iris$Sepal.Width
    xcen <- data.frame(x.low=rep(3.7,99), x.high=rep(3.7001,99))
    out.fit5 <- mle.johnsonsu(x, data.censored=xcen, plots=TRUE)
    print(out.fit5$jparms.compare)
    
}

