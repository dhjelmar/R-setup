mle.johnsonsu <- function(data, param='auto', lambda.control=2, plots=FALSE, debug=FALSE) {
    
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
    temp$description <- 'ExtDist'
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
    lambda.fix <- function(lambda, lambda.control) {
        if (lambda.control == 1) {
            ## keep lambda positive so subsequent functions are defined
            lambda <- max(1E-15, lambda)
        } else if (lambda.control == 2) {
            ## keep lambda positive so subsequent functions are defined
            lambda <- abs(lambda)
            if (lambda == 0) lambda <- 1E-15
        }
        return(lambda)
    }
    nll <- function(data, param, lambda.control=0, debug=FALSE){
        ## calculate nll (negative log likelihhod) for distribution
        x      <- data
        gamma  <- param[[1]]  
        delta  <- param[[2]]
        xi     <- param[[3]]
        lambda <- param[[4]]
        lambda <- lambda.fix(lambda, lambda.control)
        ## PDF for Johnson SU
        pdf <- delta /( lambda * sqrt(2 * pi)   ) *
            1 / sqrt(1 +            ( (x-xi)/lambda )^2)  *
            exp( -0.5*(gamma + delta * asinh( (x-xi)/lambda ))^2 )
        ## the above is equivalent
        ## pdf <- ExtDist::dJohnsonSU(x, parms=c(gamma, delta, xi, lambda))
        nll     <- -sum(log(pdf))
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
                             debug   = debug,
                             lambda.control = lambda.control,
                             control = list(trace=TRUE),
                             method  = "BFGS")
        nll.max.bestfit <- out.bestfit$value
        jparms.mle <- as.list(out.bestfit$par)
        jparms.mle$type <- 'SU'
    }, error = function(e) {
        ## what to do if error
        cat('WARNING: CONVERGENCE FAILURE IN mle.johnsonsu()\n')
        jparms <- list(gamma=NA, delta=NA, xi=NA, lambda=NA, type=NA)
    })
    ## the following is needed if use lambda.control because the
    ## parameter returned by optim() is the input pamaeter to the function
    ## not the output parameter. So if the input lambda is negative,
    ## the optim() routine thinks that was OK as the answer even though
    ## it really only ever used abs(lambda) in the optimization.
    jparms.mle$lambda <- lambda.fix(jparms.mle$lambda, lambda.control)

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
    source('setup.r')
    x <- iris$Sepal.Width
    out <- mle.johnsonsu(x, 'auto', plots=TRUE, debug=FALSE)
    print(out$jparms)
    print(out$jparms.compare)
}
