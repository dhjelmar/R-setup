loglik.weibull <- function(x=NA, xcen=NA, param=c(shape, scale), debug=FALSE){
    ## calculate log likelihhod for distribution
    if (is.data.frame(x)) x <- x[1] # convert to vector
    if (is.data.frame(xcen)) {
        ## if known data are inside xcen, move them to x and keep remainder in xcen
        xcen.na  <- xcen[ is.na(rowSums(xcen)),]                  # censored rows with NA, if any
        xcen.val <- xcen[!is.na(rowSums(xcen)),]                  # rows w/o NA, if any
        x.add   <- xcen.val[xcen.val[[1]] == xcen.val[[2]],][[1]] # known values from xcen, if any
        x <- as.numeric(na.omit(c(x, x.add)))                     # new set of known values
        xcen.lowhigh <- xcen.val[xcen.val[[1]] != xcen.val[[2]],] # censored rows with no NA
        xcen <- rbind(xcen.lowhigh, xcen.na)                      # new set ofcensored rows
    }
    ## set parameters
    shape  <- param[[1]]  
    scale  <- param[[2]]
    ## likelihood contribution from x
    pdf <- 1
    if (!is.na(x[1])) {
        pdf <- shape / scale^shape * x^(shape-1) * exp(-(x/scale)^shape)
        ## the above is equivalent
        ## pdf <- stats::dweibull(x, shape, scale)
        ## cdf <- pnorm(z)
        ## plot(z, cdf)
    }
    ## likelihood contribution from xcen
    probability <- 1
    if (is.data.frame(xcen)) {
        xcen$F.low  <- 1 - exp(- (xcen$x.low  / scale)^shape) # CDF; could use ExtDist::pWeibull(xcen$x.low, shape, scale)
        xcen$F.high <- 1 - exp(- (xcen$x.high / scale)^shape)
        ## if low CDF is NA, set to 0
        xcen$F.low[is.na(xcen$F.low)]   <- 0
        ## if high CDF is NA, set to 1
        xcen$F.high[is.na(xcen$F.high)] <- 1
        ## calculate probability for the censored interval
        probability <- xcen$F.high - xcen$F.low
    }
    loglik <- sum(log(pdf), log(probability))
    if (isTRUE(debug)) cat('shape=', signif(shape,11),
                           'scale=', signif(scale,11),
                           'loglik=', signif(loglik,11), "\n")
    return(loglik)
}        

mle.weibull <- function(x=NA, xcen=NA, param='auto', plots=FALSE, debug=FALSE) {
    
    ## weibull distribution
    ## MLE (Maximum Likelihood Estimate) fit to determine parameters
    ## LR (Likelihood Ratio) appraoch to find tolerance limit

    ## input: x     = vector of known data
    ##        xcen  = dataframe of censored data
    ##                (1st column = low value or NA; 2nd column = high value or NA)
    ##        param = initial guess for fit parameters for: shape and scale
    ##                if type is also provided, it will not be used
    ##              = 'auto' (default) uses tolerance::exttol.int() for initial guess of parameters
    ##                and will try using list(shape=1, scale=1) if that fails

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
    
    ## make sure Weibull is possible
    if (!is.na(x[1])        & min(x   , na.rm=TRUE) < 0) {
        if (!is.data.frame(xcen) | min(xcen, na.rm=TRUE) < 0) {
            cat('####################################################\n')
            cat('               ERROR: MINIMUM x < 0                 \n')
            cat('####################################################\n')
            return()
        }
    }
    
    ## if (isTRUE(plots)) par(mfrow=c(1,2))
    out <- NULL

    ##-----------------------------------------------------------------------------
    ## set initial guess for MLE fit
    if (param[1] == 'auto') {
        ## initial parameters not specified, so let R figure guess
        tol_out_weib <- tolerance::exttol.int(x.avg)
        shape   <- tol_out_weib$'shape.1'
        scale   <- tol_out_weib$'shape.2'
        param <- list(shape=shape, scale=scale)
        params.compare <- as.data.frame(param)
        params.compare$description <- 'tolerance::exttol.int'
        if (is.na(param$shape) | is.na(param$scale)) {
            ## exttol failed to converge 
            param <- list(shape = 1,
                          scale = 1)
        }
    } else {
        ## use provided initial parameters
        shape <- param$shape
        scale <- param$scale
        params.compare <- as.data.frame(param)
        params.compare$description <- 'user provided parameters'
    }
    
    temp <- as.data.frame(t(unlist(param)))
    temp$description <- 'initial guess for MLE'
    params.compare <- fastmerge(params.compare, temp)
    
    ##-----------------------------------------------------------------------------
    ## determine best fit using nll
    ## print('Attempting MLE fit on regular parameters')
    ## constraints for shape and/or scale
    ## A %*% param + B > 0
    ## row 1: shape > 0
    ## row 2: scale > 0
    A <- matrix(c(1,0, 0,1), 2, 2, byrow=TRUE)
    B <- matrix(c(0,0))
    constraints <- list(ineqA=A, ineqB=B)
    out <- NA
    out <- maxLik::maxLik(loglik.weibull,
                                  start = unlist(param),
                                  x     = x,
                                      xcen  = xcen,
                                      debug = debug,
                                      constraints = constraints,
                                      iterlim = 2000)
        parms.mle <- as.list(out$estimate)
        loglik <- out$maximum
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
        shape <- parms.mle$shape
        scale <- parms.mle$scale
        curve.points <- stats::dweibull(x, shape, scale)
        hist(x, freq=FALSE, ylim=range(out.hist$density, curve.points))
        curve(stats::dweibull(x, shape, scale), min(x), max(x), add=TRUE)
        qqplot_nwj(x, type='w', jfit=parms.mle)
    }
    
    return(list(parms=parms.mle, parms.compare=params.compare, loglik=loglik, convergence=convergence))
}



mle.weibull.test <- function() {
    source("F:\\Documents\\01_Dave's Stuff\\Programs\\GitHub_home\\R-setup\\setup.r")
    source('setup.r')
    
    fit.compare <- function(x, xcen, main=NULL) {
        ## plot histogram with no censored data and fit
        hist(x, freq=FALSE, border='black', main=main)
        out.fit0 <- mle.weibull(x, xcen=NA, plots=FALSE)
        parms0 <- out.fit0$parms
        curve(stats::dweibull(x, parms0$shape, parms0$scale), min(x), max(x), col='black', add=TRUE)
        ## plot histogram with all data treated as known fit at numeric value
        xcen.avg <- rowMeans(xcen, na.rm=TRUE) # use the average for interval data
        x.all <- c(x, xcen.avg)
        hist(x.all, freq=FALSE, border='red', add=TRUE)
        out.fit1 <- mle.weibull(x.all, xcen=NA, plots=FALSE)
        parms1 <- out.fit1$parms
        curve(stats::dweibull(x, parms1$shape, parms1$scale), min(x), max(x), col='red', add=TRUE)
        ## plot fit if treat xcen as censored
        out.fit2 <- mle.weibull(x, xcen=xcen, plots=FALSE)
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

