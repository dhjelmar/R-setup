loglik.johnsonsu <- function(x=NA, xcen=NA, param=c(gamma, delta, xi, lambda), debug=FALSE){
    ## calculate log likelihhod for distribution and given parameters
    if (is.data.frame(x)) x <- x[1] # convert to vector
    gamma  <- param[[1]]  
    delta  <- param[[2]]
    xi     <- param[[3]]
    lambda <- param[[4]]
    pdf <- 1
    if (!is.na(x[1])) {
        ## PDF for Johnson SU
        z <- gamma + delta * asinh( (x-xi)/lambda )
        pdf <- delta / (lambda*sqrt(2*pi)) / sqrt(1 +((x-xi)/lambda)^2) * exp(-0.5*z^2)
        ## the above is equivalent
        ## pdf <- ExtDist::dJohnsonSU(x, param=c(gamma, delta, xi, lambda))
        ## cdf <- pnorm(z)
        ## plot(z, cdf)
    }
    probability <- 1
    if (is.data.frame(xcen)) {
        ## following is equivalent to ExtDist::pJohnsonSU(xcen$x.low,  gamma, delta, xi, lambda)
        xcen$F.low  <- stats::pnorm(gamma+delta*asinh((xcen$x.low -xi)/lambda))
        xcen$F.high <- stats::pnorm(gamma+delta*asinh((xcen$x.high-xi)/lambda))
        ## if low CDF is NA, set to 0
        xcen$F.low[is.na(xcen$F.low)]   <- 0
        ## if high CDF is NA, set to 1
        xcen$F.high[is.na(xcen$F.high)] <- 1
        ## calculate probability for the censored interval
        probability <- xcen$F.high - xcen$F.low
    }
    loglik <- sum(log(pdf), log(probability))
    if (isTRUE(debug)) cat('gamma=', signif(gamma,11),
                           'delta=', signif(delta,11),
                           'xi   =', signif(xi, 11),
                           'lambda=', signif(lambda,11),
                           'loglik=', signif(loglik,11), "\n")
    return(loglik)
}        
        
mle.johnsonsu <- function(x=NA, xcen=NA, param='auto', plots=FALSE, debug=FALSE) {
    
    ## johnsonsu distribution
    ## MLE (Maximum Likelihood Estimate) fit to determine parameters
    ## LR (Likelihood Ratio) appraoch to find tolerance limit

    ## input: x     = vector of known data
    ##        xcen  = dataframe of censored data
    ##                (1st column = low value or NA; 2nd column = high value or NA)
    ##        param = initial guess for fit parameters for: gamma, delta, xi, and lambda
    ##                if type is also provided, it will not be used
    ##              = 'auto' (default) uses ExtDist::eJohnsonSU() for initial guess of parameters,
    ##                will switch to SuppDists::JohnsonFit(x) if that fails, and will switch to
    ##                try using list(gamma=xx, delta=xx, xi=xx, lambda=xx) if that fails

    ## based on approach found here:
    ## https://www.r-bloggers.com/2019/08/maximum-likelihood-estimation-from-scratch/
    ## https://personal.psu.edu/abs12/stat504/Lecture/lec3_4up.pdf

    ## alternately, could probably use bbmle::mle2() to determine coefficients
    ## an example is here: https://stackoverflow.com/questions/45208176/the-weibull-distribution-in-r-extdist

    if (is.data.frame(x)) x <- x[1] # convert to vector
  
    if (is.data.frame(xcen)) {
        ## censored data also provided (only reason for following is if
        ## x.low and x.high were not the names of the two columns of data)
        xcen <- data.frame(x.low = xcen[[1]], x.high = xcen[[2]])

        ## calculate average value, ignoring NA, xcen to for use in estimating parameters
        ## from packages that do not have censor capability
        xcen.avg <- rowMeans(xcen, na.rm=TRUE)

    } else {
        xcen.avg <- NA
    }
    x.avg <- as.numeric( na.omit( c(x, xcen.avg) ) )
    
    
    ## if (isTRUE(plots)) par(mfrow=c(1,2))
    out <- NULL
 
    ##-----------------------------------------------------------------------------
    ## let R figure out which Johnson distribution fits best
    parms.SuppDists <- SuppDists::JohnsonFit(x.avg)
    ## params.compare <- as.data.frame(t(unlist(parms.SuppDists)))
    params.compare <- as.data.frame(parms.SuppDists)
    params.compare$description <- 'SuppDists::JohnsonFit'

    ##-----------------------------------------------------------------------------
    ## force the Johnson SU distribution
    parms.out <- ExtDist::eJohnsonSU(x.avg)
    parms.ExtDist <- list(gamma   = parms.out$gamma,
                           delta   = parms.out$delta,
                           xi      = parms.out$xi,
                           lambda  = parms.out$lambda,
                           type    = 'SU')
    temp <- as.data.frame(t(unlist(parms.ExtDist)))
    temp$description <- 'ExtDist::eJohnsonSU'
    params.compare <- fastmerge(params.compare, temp)

    ##-----------------------------------------------------------------------------
    ## set initial guess for MLE fit
    if (param[1] == 'auto') {

        ## force the Johnson SU distribution
        parms <- parms.ExtDist

        if (is.na(parms$gamma) |
            (parms$gamma == -0.5 & parms$delta == 2.0 & parms$xi == -0.5 & parms$lambda == 2.0)) {
            ## ExtDist failed to converge so try SuppDists
            parms <- parms.SuppDists
            if (is.na(parms$gamma) | parms$type != 'SU') {
                ## SuppDists::JohnsonFit(x) failed to converge or did not return type=SU
                ## so arbitrarily set parms to something as a starting point
                parms <- list(gamma   = 1,
                               delta   = 1,
                               xi      = 1,
                               lambda  = 1,
                               type    = 'SU')
            }
         }
        
    } else {
        ## use Johnson parameters specified in param
        ## needs to be in same list format as created by SuppDists::JohnsonFit
        parms <- param
    }
    temp <- as.data.frame(t(unlist(parms)))
    temp$description <- 'Initial guess for MLE'
    params.compare <- fastmerge(params.compare, temp)

    ## strip off type if that is provided in param list since not needed for MLE fit
    param <- parms[1:4]   # list(gamma, delta, xi, lambda)
    
    ##-----------------------------------------------------------------------------
    ## determine best fit using nll
    ## print('Attempting MLE fit on regular parameters')
    ## constraints for gamma, delta, xi, and/or lambda
    ## A %*% param + B > 0
    ## row 1: delta  > 0
    ## row 2: lambda > 0
    A <- matrix(c(0,1,0,0,  0,0,0,1), 2, 4, byrow=TRUE)
    B <- matrix(c(0,0),               2, 1)
    constraints <- list(ineqA=A, ineqB=B)
    out <- NA
    out <- maxLik::maxLik(loglik.johnsonsu,
                          start = unlist(param),
                          x     = x,
                          xcen  = xcen,
                          debug = debug,
                          constraints = constraints,
                          iterlim = 2000)
    parms.mle      <- as.list(out$estimate)
    parms.mle$type <- 'SU'
    loglik         <- out$maximum
    convergence <- if (out$message == 'successful convergence ') {'successful'}
                   else {out$message}
    if (convergence != 'successful') {
        cat('###############################################\n')
        cat('WARNING: CONVERGENCE FAILURE IN mle.johnsonsu()\n')
        cat('###############################################\n')
    }
    
    ## add MLE parameters to params.compare dataframe
    temp <- as.data.frame(parms.mle)
    temp$description <- 'MLE'
    params.compare <- fastmerge(params.compare, temp)

    ## remove factor levels from params.compare
    params.compare <- droplevels.all(params.compare)

    if (isTRUE(plots)) {
        out.hist <- hist(x, plot=FALSE)
        curve.points <- ExtDist::dJohnsonSU(x, params=parms.mle)
        hist(x, freq=FALSE, ylim=range(out.hist$density, curve.points))
        curve(ExtDist::dJohnsonSU(x, params=parms.mle), min(x), max(x), add=TRUE)
        qqplot_nwj(x, type='j', jfit=parms.mle)
    }
    
    return(list(parms=parms.mle, parms.compare=params.compare, loglik=loglik, convergence=convergence))
}



mle.johnsonsu.test <- function() {
    source("F:\\Documents\\01_Dave's Stuff\\Programs\\GitHub_home\\R-setup\\setup.r")
    source('setup.r')
    
    fit.compare <- function(x, xcen, main=NULL) {
        ## plot histogram with no censored data and fit
        hist(x, freq=FALSE, border='black', main=main)
        out.fit0 <- mle.johnsonsu(x, xcen=NA, plots=FALSE)
        parms0 <- out.fit0$parms
        curve(ExtDist::dJohnsonSU(x, params = out.fit0$parms), min(x), max(x), col='black', add=TRUE)
        ## plot histogram with all data treated as known fit at numeric value
        xcen.avg <- rowMeans(xcen, na.rm=TRUE) # use the average for interval data
        x.all <- c(x, xcen.avg)
        hist(x.all, freq=FALSE, border='red', add=TRUE)
        out.fit1 <- mle.johnsonsu(x.all, xcen=NA, plots=FALSE)
        parms1 <- out.fit1$parms
        curve(ExtDist::dJohnsonSU(x, params = out.fit1$parms), min(x), max(x), col='red', add=TRUE)
        ## plot fit if treat xcen as censored
        out.fit2 <- mle.johnsonsu(x, xcen=xcen, plots=FALSE)
        parms2 <- out.fit2$parms
        curve(ExtDist::dJohnsonSU(x, params = out.fit2$parms), min(x), max(x), col='blue', type='p', add=TRUE)
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

