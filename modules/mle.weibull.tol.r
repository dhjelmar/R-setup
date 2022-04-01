mle.weibull.tol <- function(data, data.censored=NA, param='auto', param.control=2,
                              side.which='upper', sided=1, conf=0.99, alpha=NULL, P=0.99,
                              plots=FALSE, plots.nr=FALSE, debug=FALSE, main=NULL) {
    
    ## weibull distribution
    ## MLE (Maximum Likelihood Estimate) fit to determine parameters
    ## LR (Likelihood Ratio) appraoch to find tolerance limit

    ## input: data  = vector of data
    ##        param = initial guess for fit parameters for: shape and scale
    ##              = 'auto' (default) uses tolerance::exttol.int() for initial guess of parameters
    ##                and will try using list(shape=1, scale=1) if that fails
    ##        side.which = 'upper' or 'lower' (not used if sided = 2)
    ##        sided = 1 (default) means 1-sided tolerance limit
    ##              = 2 means 2-sided tolerance limit
    ##        P     = proportion (or coverage)
    ##              = any value if 1-sided
    ##                Note: lower 99/99 means  the lower bound on P=0.99
    ##                      if you really want the lower bound on P=0.01, then specify P=0.1
    ##              >= 0.5 if 2-sided
    ##        conf  = confidence used to determine chi-square
    ##        alpha = NULL (default) sets alpha = 2*(1-conf)/sided
    ##                e.g., if 1-sided conf = 99%:
    ##                          alpha = 2(1-0.99)/1 = 0.02
    ##                          chi-square = qchisq(1-alpha,1) = 5.411894
    ##                                     = 2-sided 98% confidence limit
    ##                                     = 1-sided 99% confidence limit <- which is needed
    ##              = # uses input value to overwrite value based on conf
    
    x <- data
    
    if (is.data.frame(data.censored[1])) {
        ## censored data also provided
        xcen <- data.frame(x.low = data.censored[[1]], x.high = data.censored[[2]])
    } else {
        xcen <- NA
    }

    if (is.null(alpha)) {
        ## set alpha level for chi-square for use in confidence limit calculation
        alpha <- 2*(1-conf)/sided
    } else {
        ## calculate confidence limit from alpha used in chi-square
        conf <- 1 - alpha * sided/2
    }
    
    ## if (isTRUE(plots) & sided == 2) par(mfrow=c(1,2))
    out <- NULL

    ##-----------------------------------------------------------------------------
    ## Weibull SU parameters
    out <- mle.weibull(x, xcen, plots=FALSE)
    params <- out$parms
    shape  <- params$shape
    scale  <- params$scale
    params.save <- params
    
    ##-----------------------------------------------------------------------------
    ## redefine nll function to fit on desired quantile to find standard error
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
    nll.q <- function(data, data.censored=NA, param, P, param.control, debug=FALSE) {
        ## calculate nll (negative log likelihhod) for distribution
        x       <- data
        xcen    <- data.censored
        quant  <- param[[1]]  # replaced gamma with quant as a parameter
        shape  <- param[[2]]
        scale  <- quant/(-log((1-P)))^(1/shape)
        scale  <- param.fix(scale, param.control)
        pdf <- shape / scale^shape * x^(shape-1) * exp(-(x/scale)^shape)
        ## above is equivalent to
        ## pdf <- stats::dweibull(x, shape, scale)
        if (is.data.frame(xcen)) {
            xcen$F.low  <- stats::pweibull(xcen$x.low,  shape, scale)
            xcen$F.high <- stats::pweibull(xcen$x.high, shape, scale)
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
        if (isTRUE(debug)) cat('quant=', signif(quant,11),
                               'shape=', signif(shape,11),
                               'scale=', signif(scale,11),
                               'nll=', signif(nll,11), "\n")
        return(nll)
    }
    
    ##-----------------------------------------------------------------------------
    ## find confidence limit at level alpha for requested coverage, P
    tol.limits <- NA
    tol.approx <- NA
    params.q.save <- NA
    k <- 0
    xyplot <- NA
    quant.P.save <- NA
    P.in <- P
    if (sided == 1) {
        P.values <- P
    } else {
        ## sided == 2
        P.values <- c((1-P)/2, 1-(1-P)/2)
        ## e.g., if P = 0.99, then will base tolerance limits on the
        ##       quantiles associated with 0.5% and 99.5% coverage
    }
    for (P in P.values) {

        k <- k+1  # counter

        ##----------------------
        ## fist set estimated parameters
        quant.P <- stats::qweibull(P, params$shape, params$scale)
        quant.P.save[k] <- quant.P
        quant.param <- c(quant=quant.P, params$shape)

        ##----------------------
        ## calculate confidence limits using LR (Likelihood Ratio)
        ## confidene limit is defined at likelihood that is lower than max by chi-squared
        ll.max <- -nll.q(x, xcen, quant.param, P, param.control=param.control)
        ll.tol <-  ll.max - qchisq(1 - alpha, 1)/2   # qchisq(1-0.02, 1) = 5.411894
        cat('MLE=', ll.max, '; tolerance limit at MLE=', ll.tol, '\n')

        ##----------------------
        ## refit on alternate parameters to determine standard error for fit on quantile
        cat('Attempting MLE fit on alternate parameters for P=', P, '\n')
        out.bestfit.q <- optim(par     = quant.param, 
                               fn      = nll.q, 
                               data    = x,
                               data.censored = xcen,
                               P       = P,
                               param.control = param.control,
                               debug   = debug,
                               control = list(trace=TRUE),
                               hessian = TRUE,
                               method  = "BFGS")
        nll.max.bestfit.q <- out.bestfit.q$value
        quant.P  <- out.bestfit.q$par[[1]]
        shape.P  <- out.bestfit.q$par[[2]]
        scale.P  <- quant.P/(-log((1-P)))^(1/shape.P)
        scale.P  <- param.fix(scale.P, param.control)
        params.q <- as.list(out.bestfit.q$par)
        params.q$scale <- scale.P
        params.q.save[k] <- list(params.q)
        print(as.data.frame(params.q))
        standard.error <- as.numeric( sqrt(diag(solve(out.bestfit.q$hessian))) )
        print(standard.error)
        cat('standard error:', standard.error, '\n')
        cat('\n')
        
        ##----------------------
        ## estimate confidence limit using standard error to serve as starting point for search
        ## estimates      <- as.numeric( out.bestfit.q$par )
        ## coef           <- data.frame(estimates, standard.error)
        ## rownames(coef) <- names(out.bestfit.q$par)
        dof    <- length(x) - 2    # 2 independent fitting parameters in Weibull
        student.t <- qt(1 - (1-conf)/sided, dof)  # 2.3326 for dof=598
        quant.P.alpha.l.guess <- quant.P - student.t * standard.error[1]
        quant.P.alpha.u.guess <- quant.P + student.t * standard.error[1]
        cat('Initial guesses for confidence interval for P=', P, '\n')
        cat(quant.P.alpha.l.guess, quant.P.alpha.u.guess, '\n\n')
        tol.approx <- c(tol.approx, quant.P.alpha.l.guess, quant.P.alpha.u.guess)


        ##----------------------
        if (isTRUE(plots)) {
            quant.dif <- quant.P.alpha.u.guess - quant.P
            xmin <- quant.P.alpha.l.guess - 1.1*quant.dif
            xmax <- quant.P.alpha.u.guess + 1.1*quant.dif
            xplot <- seq(xmin, xmax, length.out=101)
            yplot <- NA
            plot(quant.P.save[k], ll.max, col='blue', pch=16, cex=2,
                 xlab='quantile', ylab='log likelihood',
                 xlim=range(xmin, xmax),
                 ylim=range(ll.max, ll.tol),
                 main=main)
            abline(h=ll.tol)
            abline(v=c(quant.P.alpha.l.guess, quant.P.alpha.u.guess), col='red', lty=2)
        }
        
        ##----------------------
        ## function to fit on delta, xi and lambda
        nll.fixedq <- function(data, data.censored=NA, param, quant, P,
                               param.control=param.control, debug=FALSE) {
            ## calculate nll (negative log likelihhod) for distribution
            ## for specified quant (i.e., only fit delta, xi, and lambda)
            param <- list(quant=quant, shape=param[[1]])
            nll.q(data, data.censored, param, P, param.control, debug=FALSE)
        }
        ll.fixedq <- function(x0, data, data.censored=NA,
                              P, shape=shape.P,
                              param.control=param.control, debug=FALSE) {
            ## first determine best fit delta, xi, and lambda for given x0=quant (and P)
            fit <- NA
            tryCatch({
                fit <- optim(par     = shape, 
                             fn      = nll.fixedq, 
                             data    = data,
                             data.censored = data.censored,
                             quant   = x0,
                             P       = P,
                             param.control = param.control,
                             debug   = debug,
                             control = list(trace=TRUE,
                                            maxit=1e4),
                             hessian = TRUE,
                             method  = "BFGS")
                ll <- -fit$value
            }, error = function(e) {
                ## what to do if error
                cat('WARNING: CONVERGENCE FAILURE IN mle.weibull.tol.r MODULE\n')
                cat('         IN ll.fixedq OPTIMIZATION ROUTINE.\n\n')
                fit <- NA
                ll  <- NA
            })
            if (isFALSE(debug)) {
                return(ll)
            } else {
                return(list(fit=fit, ll=ll))
            }
        }

        ##----------------------
        if (isTRUE(plots)) {
            xplot <- seq(xmin, xmax, length.out=301)
            yplot <- NA
            for (ploti in 1:length(xplot)) {
                ## if (ploti == 3) browser()
                yplot[ploti] <- ll.fixedq(xplot[ploti], x, xcen, P, shape,
                                          param.control=param.control, debug=FALSE)
                cat(ploti, xplot[ploti], yplot[ploti], '\n')
            }
            xyplot <- data.frame(quantile       = xplot[1:length(yplot)],
                                 log.likelihood = yplot,
                                 ll.tol         = ll.tol)
            points(xyplot$quantile, xyplot$log.likelihood)
            ## converged <- which(!is.na(yplot))
            ## points(xplot[converged], yplot[converged],
            ##        xlab='quantile', ylab='log likelihood', col='red')
            plotit <- function(xplot=0.994, col='blue', pch=16, cex=2) {
                ## this is just for using interactively if debugging
                yplot <- ll.fixedq(x0     = xplot,
                                   data   = x,
                                   data.censored = xcen,
                                   P      = P,
                                   shape  = shape.P,
                                   param.control = param.control,
                                   debug  = FALSE)
                points(xplot, yplot, col=col, pch=16, cex=2)
                return(yplot)
            }
        }

        ## determine confidence bound as point where the likelihood ratio equals 11.tol
        quant.P.alpha.l <- NA
        quant.P.alpha.u <- NA
        if (side.which == 'lower' | (P < 0.5 & sided == 2)) {
            ## find lower tolerance limit
            out.nrl <- newton.raphson(f = ll.fixedq,
                                      xguess = quant.P.alpha.l.guess,
                                      ytarget = ll.tol,
                                      data   = x,
                                      data.censored = xcen,
                                      P      = P,
                                      shape  = shape,
                                      param.control = param.control,
                                      tol = 1e-5, n = 1000,
                                      plots=plots.nr)
            quant.P.alpha.l <- out.nrl$root
            if (is.null(quant.P.alpha.l)) { quant.P.alpha.l <- NA }
        } else {
            ## find upper tolerance limit
            out.nru <- newton.raphson(f = ll.fixedq,
                                      xguess = quant.P.alpha.u.guess,
                                      ytarget = ll.tol,
                                      data   = x,
                                      data.censored = xcen,
                                      P      = P,
                                      shape  = shape,
                                      param.control = param.control,
                                      tol = 1e-5, n = 1000,
                                      plots=plots.nr)
            quant.P.alpha.u <- out.nru$root
            if (is.null(quant.P.alpha.u)) { quant.P.alpha.u <- NA }
        }
        cat('Final confidence interval for P=', P, '\n')
        cat(quant.P.alpha.l, quant.P.alpha.u, '\n\n')

        ## collect tolerance limit calculations (i.e., for 1-P and P)
        tol.limits <- c(tol.limits, quant.P.alpha.l, quant.P.alpha.u)
        
        if (isTRUE(plots)) {
            ## plot intersection with log likelihood curve
            if (!is.na(quant.P.alpha.l)) points(quant.P.alpha.l, ll.tol, col='red', pch=16, cex=2)
            if (!is.na(quant.P.alpha.u)) points(quant.P.alpha.u, ll.tol, col='red', pch=16, cex=2)
        }

    }

    ## collect lower and upper tolerance limits
    tol.approx <- range(tol.approx, na.rm=TRUE)
    tol.limits <- range(tol.limits, na.rm=TRUE)
    if (sided == 2) {
        tol.lower <- tol.limits[1]
        tol.upper <- tol.limits[2]
    } else if (side.which == 'lower') {
        ## lower, 1-sided
        tol.lower <- tol.limits[1]
        tol.upper <- NA
    } else {
        ## upper, 1-sided
        tol.lower <- NA
        tol.upper <- tol.limits[2]
    }
    ## reset P to input value
    P <- P.in

    ## collect tolerance values in dataframe similar to extol.int for weibull
    tolerance <- data.frame(sided, alpha.chisq=alpha, conf, P, tol.lower, tol.upper)

    ## print final parameter comparison and tolerance limits
    print(as.data.frame(params.save))
    cat('\n')
    print(tolerance)
    cat('\n')
    
    return(list(ll.plot    = xyplot,
                ll.max     = ll.max,
                ll.tol     = ll.tol,
                quant.P    = quant.P.save,
                params     = params.save,
                tol.approx = tol.approx,
                tolerance  = tolerance))
}


mle.weibull.tol.test <- function() {
    source('setup.r')

    ##-------------------------------------------------------------------
    ## create Weibull dataset
    set.seed(1)
    parms <- list(shape=1, scale=1)
    parms <- list(shape=2, scale=1)
    x <- stats::rweibull(1000, parms$shape, parms$scale)

    ## determine weibull parameters
    plotspace(2,2)
    out.fit <- mle.weibull(x, plots=TRUE)
    parms.mle <- out.fit$parms
    print(out.fit$parms.compare)

    ## weibull tolerance limits
    out.tol <- mle.weibull.tol(x, plots=TRUE)

    ## incorporated into other functions
    plotspace(2,2)
    out.hist <- hist_nwj(x)
    out <- qqplot_nwj(x, type='n')
    out <- qqplot_nwj(x, type='w')
    out <- qqplot_nwj(x, type='j')

    ##-------------------------------------------------------------------
    ## consider the following dataset
    x <- iris$Sepal.Width
    plotspace(2,2)
    ## lower tolerance limit
    out.lower <- mle.weibull.tol(data=x, param='auto', param.control=2,
                                   side.which='lower', sided=1, conf=0.99, P=0.01, 
                                   plots=TRUE, plots.nr=FALSE, debug=FALSE, main='lower, 1-sided 99/1')
    ## using default parameters except for 'plots'
    out.upper <- mle.weibull.tol(data=x, param='auto', param.control=2,
                                   side.which='upper', sided=1, conf=0.99, P=0.99, 
                                   plots=TRUE, plots.nr=FALSE, debug=FALSE, main='upper, 1-sided 99/99')
    ## lower and upper 2-sided tolerance limits
    ## test that 1-sided 99/99 is the same as a 2-sided 98/98
    out.twosided <- mle.weibull.tol(data=x, param='auto', param.control=2,
                                      side.which='both', sided=2, conf=0.98, P=0.98, 
                                      plots=TRUE, plots.nr=FALSE, debug=FALSE, main='2-sided 98/98')

    
    ## test that upper and lower 1-sided 75/99 are the same as a 2-sided 50/98  
    plotspace(2,2)
    ## lower and upper 1-sided tolerance limits
    out.lower <- mle.weibull.tol(data=x, param='auto', param.control=2,
                                   side.which='lower', sided=1, conf=0.99, P=0.25, 
                                   plots=TRUE, plots.nr=FALSE, debug=FALSE, main='lower, 1-sided 99/25')
    out.upper <- mle.weibull.tol(data=x, param='auto', param.control=2,
                                   side.which='upper', sided=1, conf=0.99, P=0.75, 
                                   plots=TRUE, plots.nr=FALSE, debug=FALSE, main='upper, 1-sided 99/75')
    ## lower and upper 2-sided tolerance limits
    out.twosided <- mle.weibull.tol(data=x, param='auto', param.control=2,
                                      side.which='both', sided=2, conf=0.98, P=0.5, 
                                      plots=TRUE, plots.nr=FALSE, debug=FALSE, main='2-sided 98/50')


    ##-------------------------------------------------------------------
    ## TEST OUT CENSORED DATA FUNCTIONALITY
    ##-------------------------------------------------------------------
    fit.compare.cen <- function(x, xcen, main=NULL) {
        ## plot histogram with no censored data and fit
        hist(x, freq=FALSE, border='black', main=main, xlim=c(2,5.4))
        out.fit0 <- mle.weibull.tol(x, data.censored=NA, plots=FALSE)
        parms0 <- out.fit0$params
        curve(stats::dweibull(x, parms0$shape, params0$scale), min(x), max(x), col='black', add=TRUE)
        abline(v=out.fit0$tolerance$tol.upper, col='black', lty=2)
        ## plot histogram with all data treated as known fit at numeric value
        xcen.avg <- rowMeans(xcen, na.rm=TRUE) # use the average for interval data
        x.all <- c(x, xcen.avg)
        hist(x.all, freq=FALSE, border='red', add=TRUE)
        out.fit1 <- mle.weibull.tol(x.all, data.censored=NA, plots=FALSE)
        parms1 <- out.fit1$params
        curve(stats::dweibull(x, params1$shape, params1$scale), min(x), max(x), col='red', add=TRUE)
        abline(v=out.fit1$tolerance$tol.upper, col='red', lty=2)
        ## plot fit if treat xcen as censored
        out.fit2 <- mle.weibull.tol(x, data.censored=xcen, plots=FALSE)
        parms2 <- out.fit2$params
        curve(stats::dweibull(x, parms2$shape, parms2$scale), min(x), max(x), col='blue', type='p', add=TRUE)
        abline(v=out.fit2$tolerance$tol.upper, col='blue', lty=2, lwd=2)
        ## add legend
        legend('topright', 
               legend=c('known', 'all',   'censored'),
               col   =c('black', 'red', 'blue'),
               lty   =c( 1     ,  1   ,  NA),
               pch   =c( NA    ,  NA  ,  1))

        ## return table of results
        df0 <- as.data.frame(out.fit0$params)
        df0$tol.upper <- out.fit0$tolerance$tol.upper
        df1 <- as.data.frame(out.fit1$params)
        df1$tol.upper <- out.fit1$tolerance$tol.upper
        df2 <- as.data.frame(out.fit2$params)
        df2$tol.upper <- out.fit2$tolerance$tol.upper
        df  <- rbind(df0, df1, df2)
        return(df)
    }

    x <- iris$Sepal.Width
    plotspace(4,2)
    ## left censored data should reduce the upper 99/99
    ## following shows that did not happen (99/99 for fit with xcen is >> that for fit to x)
    xnum <- 10
    x.low  <- NA
    x.high <- 2.0
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit1a <- fit.compare.cen(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    print(fit1)
    ## it may be that the fit was bad because the following does pull the 99/99 in as I would expect
    ## out.qq <- qqplot_nwj(x, type='j', mainadder = 'x only')
    ## x.all <- na.omit(c(x, xcen$x.low, xcen$x.high))
    ## out.qq <- qqplot_nwj(x.all, type='j', mainadder='all')
    ## out.qq <- qqplot_nwj(x, xcen, type='j', mainadder='censored')
    x.all <- na.omit(c(x, xcen$x.low, xcen$x.high))
    out.qq <- qqplot_nwj(x.all, type='j', mainadder='all')
    
    x.high <- 2.14
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit1b <- fit.compare.cen(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    x.all <- na.omit(c(x, xcen$x.low, xcen$x.high))
    out.qq <- qqplot_nwj(x.all, type='j', mainadder='all')
    
    x.high <- 2.15
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit1c <- fit.compare.cen(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    out.qq <- qqplot_nwj(x, xcen, type='j', mainadder='censored')
    
    x.high <- 2.2
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit1d <- fit.compare.cen(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    out.qq <- qqplot_nwj(x, xcen, type='j', mainadder='censored')

    plotspace(1,1)
    ## data that are high should push out the upper 99/99
    ## following shows that is true (99/99 for fit with xcen is >> that for fit to x)
    x.low  <- 3.7
    x.high <- NA
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit2 <- fit.compare.cen(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    print(fit2)
    
    ## data that are unknown over the middle range of data should strengthen confidence in 99/99
    ## following shows that is true (99/99 for fit with xcen is < that for fit to x)
    x.low  <- 2.2
    x.high <- 3.7
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit3 <- fit.compare.cen(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    print(fit3)
    
    ## data that are unknown between almost the lowest value and Inf should have no real impact
    ## following shows that is true (99/99 for fit with xcen is about the same as that for fit to x)
    x.low  <- 2.2
    x.high <- NA
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit4 <- fit.compare.cen(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    print(fit4)
   
}

