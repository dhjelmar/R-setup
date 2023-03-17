loglik.weibull.q <- function(x=NA, xcen=NA, param=c(quant, shape), P, debug=FALSE){
    ## calculate log likelihhod for distribution and given parameters
    quant  <- param[[1]]  # replaced gamma with quant as a parameter
    shape  <- param[[2]]
    scale  <- quant/(-log((1-P)))^(1/shape)
    loglik.weibull(x, xcen, param=c(shape, scale), debug=FALSE)
}

loglik.weibull.q.set <- function(x=NA, xcen=NA, param=shape, quant, P, debug=FALSE){
    ## calculate log likelihhod for distribution and given parameters
    shape  <- param[[1]]
    scale  <- quant/(-log((1-P)))^(1/shape)
    loglik.weibull(x, xcen, param=c(shape, scale), debug=FALSE)
}

## loglik.weibull.q.set.optim <- function(shape, x=NA, xcen=NA, quant, P, debug=FALSE){
##     ## calculate log likelihhod for distribution and given parameters
##     ## first parameter is the optimization parmaeter for optimize()
##     ## browser()
##     ## shape  <- param[[1]]
##     scale  <- quant/(-log((1-P)))^(1/shape)
##     loglik.weibull(x, xcen, param=c(shape, scale), debug=FALSE)
## }

mle.weibull.tol <- function(x, xcen=NA, param='auto',
                              side.which='upper', sided=1, conf=0.99, alpha.eff=NULL, P=0.99,
                              plots=FALSE, plots.nr=FALSE, debug=FALSE, main=NULL) {
    
    ## weibull distribution
    ## MLE (Maximum Likelihood Estimate) fit to determine parameters
    ## LR (Likelihood Ratio) appraoch to find tolerance limit
    ## e.g., 95% upper confidence bound on the 99th percentile
    ## wikipedia     convention: "100×p%    / 100×(1−α) tolerance interval"
    ## more standard convention: "100×(1−α) / 100×p%    tolerance interval"

    ## input: x     = vector of known data
    ##        xcen  = dataframe of censored data
    ##                (1st column = low value or NA; 2nd column = high value or NA)
    ##        param = initial guess for fit parameters for: shape and scale
    ##              = 'auto' (default) uses tolerance::exttol.int() for initial guess of parameters
    ##                and will try using list(shape=1, scale=1) if that fails
    ##        side.which = 'upper' or 'lower' (not used if sided = 2)
    ##        sided = 1 (default) means 1-sided tolerance limit
    ##              = 2 means 2-sided tolerance limit
    ##        conf  = confidence used to determine chi-square
    ##        alpha.eff = NULL (default) sets effective alpha, alpha.eff = 2*(1-conf)/sided, for use with qchisq();
    ##                    Staticticians define alpha as 1 - confidence, but functions and tables require
    ##                    use of an effective alpha that may be different.
    ##                    e.g., if 1-sided conf = 99%:
    ##                          alpha.eff = 2(1-0.99)/1 = 0.02
    ##                          chi-square = qchisq(1-alpha.eff,1) = 5.411894
    ##                                     = 2-sided 98% confidence limit
    ##                                     = 1-sided 99% confidence limit <- which is needed
    ##              = # uses input value to overwrite value based on conf
    ##        P     = proportion (or coverage)
    ##              = any value if 1-sided
    ##                Note: lower 99/99 means  the lower bound on P=0.99
    ##                      if you really want the lower bound on P=0.01, then specify P=0.1
    ##              >= 0.5 if 2-sided
    
    if (is.data.frame(x)) x <- x[[1]] # convert to vector

    if (is.data.frame(xcen[1])) {
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

    if (is.null(alpha.eff)) {
        ## set alpha.eff level for chi-square for use in confidence limit calculation
        alpha.eff <- 2*(1-conf)/sided
    } else {
        ## calculate confidence limit from alpha.eff used in chi-square
        conf <- 1 - alpha.eff * sided/2
    }
    
    ## if (isTRUE(plots) & sided == 2) par(mfrow=c(1,2))
    out <- NULL

    ##-----------------------------------------------------------------------------
    ## Weibull SU parameters
    out <- mle.weibull(x, xcen, param=param, plots=FALSE)
    params <- out$parms
    shape  <- params$shape
    scale  <- params$scale
    params.save <- params
    params$type <- NULL
    out.all  <- as.data.frame(params)
    out.all$quant <- NA
    out.all$P     <- NA
    out.all$loglik <- out$loglik
    out.all$convergence <- out$convergence
    out.all$optimizer    <- 'maxLik'
    out.all$max.function <- 'loglik.weibull'
    out.all$fit.params   <- 'scale, shape'
    cat('convergence for best estimate parameters:', out.all$convergence, '\n')
    cat('\n')
   
    ##-----------------------------------------------------------------------------
    ## find confidence limit at level alpha.eff for requested coverage, P
    tol.limits <- NA
    tol.approx <- NA
    ## params.q.save <- NA
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
        quant.P.orig <- stats::qweibull(P, params$shape, params$scale)
        quant.P.save[k] <- quant.P.orig
        quant.param <- list(quant=quant.P.orig, shape=params$shape)

        ##----------------------
        ## calculate confidence limits using LR (Likelihood Ratio)
        ## confidene limit is defined at likelihood that is lower than max by chi-squared
        loglik.max <- loglik.weibull.q(x, xcen, quant.param, P)
        loglik.tol <- loglik.max - qchisq(1 - alpha.eff, 1)/2   # qchisq(1-0.02, 1) = 5.411894
        cat('P =', P, 'MLE =', loglik.max, '; tolerance limit at MLE =', loglik.tol, '\n')
        temp <- data.frame(shape=shape, scale=scale,
                           quant=quant.P.orig, P=P, loglik=loglik.max, convergence=NA,
                           optimizer    = NA,
                           max.function = NA,
                           fit.params   = 'quant calculated for P')
        out.all <- rbind(out.all, temp)
        cat('\n')

        ##----------------------
        ## refit on alternate parameters to determine standard error for fit on quantile
        cat('Attempting MLE fit on alternate parameters for P=', P, '\n')
        ## constraints for shape and/or scale
        ## A %*% param + B > 0
        ## row 1: shape > 0
        A <- matrix(c(0,1), 1, 2, byrow=TRUE)
        B <- 0
        constraints <- list(ineqA=A, ineqB=B)
        out.qs <- NA
        out.qs <- maxLik::maxLik(loglik.weibull.q,
                                   start = unlist(quant.param),  # quant, shape
                                   x     = x,
                                   xcen  = xcen,
                                   P     = P,
                                   debug = debug,
                                   constraints = constraints,
                                   method = 'BFGS',
                                   iterlim = 2000)
        print(summary(out.qs))
        convergence.qs <- if (out.qs$code == 0) {'successful'}
                          else                  {out.qs$message}
        if (convergence.qs != 'successful') {
            cat('####################################################################################\n')
            cat('WARNING: CONVERGENCE FAILURE IN mle.weibull.tol when maximizing loglik.weibull.q\n')
            cat('####################################################################################\n')
        }
        loglik.max.qs <- out.qs$maximum
        params.qs <- as.list(out.qs$estimate)
        quant.P  <- params.qs[[1]]
        shape.P  <- params.qs[[2]]
        scale.P  <- quant.P/(-log((1-P)))^(1/shape.P)
        params.qs$scale <- scale.P
        ## params.q.save[k] <- list(params.qs)
        print(as.data.frame(params.qs))
        ## following does not work for maxLik because hessian often has negative diagonals
        ## standard.error <- as.numeric( sqrt(diag(solve(out.qs$hessian))) )
        standard.error <- summary(out.qs)$estimate['quant', 'Std. error']
        temp <- data.frame(shape=shape.P, scale=scale.P,
                           quant=quant.P, P=P, loglik=loglik.max.qs, convergence=convergence.qs,
                           optimizer    = 'maxLik',
                           max.function = 'loglik.weibull.q',
                           fit.params   = 'quant, shape')
        out.all <- rbind(out.all, temp)
        cat('convergence for quant, shape:', convergence.qs, '\n')
        cat('\n')
        
        ##----------------------
        ## estimate confidence limit using standard error to serve as starting point for search
        ## estimates      <- as.numeric( out.qs$par )
        ## coef           <- data.frame(estimates, standard.error)
        ## rownames(coef) <- names(out.qs$par)
        dof    <- length(x.avg) - 2    # 2 independent fitting parameters in Weibull
        student.t <- qt(1 - (1-conf)/sided, dof)  # 2.3326 for dof=598
        if (is.nothing(standard.error) | standard.error == Inf) standard.error <- 0.001 * quant.P
        quant.P.alpha.eff.l.guess <- quant.P - student.t * standard.error
        quant.P.alpha.eff.u.guess <- quant.P + student.t * standard.error
        cat('Initial guesses for confidence interval for P=', P, '\n')
        cat(quant.P.alpha.eff.l.guess, quant.P.alpha.eff.u.guess, '\n\n')
        tol.approx <- c(tol.approx, quant.P.alpha.eff.l.guess, quant.P.alpha.eff.u.guess)


        ##----------------------
        if (isTRUE(plots)) {
            quant.dif <- quant.P.alpha.eff.u.guess - quant.P
            xmin <- quant.P.alpha.eff.l.guess - 1.1*quant.dif
            xmax <- quant.P.alpha.eff.u.guess + 1.1*quant.dif
            xplot <- seq(xmin, xmax, length.out=101)
            yplot <- NA
            plot(quant.P.save[k], loglik.max, col='blue', pch=16, cex=2,
                 xlab='quantile', ylab='log likelihood',
                 xlim=range(xmin, xmax),
                 ylim=range(loglik.max, loglik.tol),
                 main=main)
            abline(h=loglik.tol)
            abline(v=c(quant.P.alpha.eff.l.guess, quant.P.alpha.eff.u.guess), col='red', lty=2)
        }
        
        ##----------------------
        ## function to fit on shape
        loglik.fixedq <- function(quant, data, xcen, P, shape, debug=FALSE) {
            ## first determine best fit delta, xi, and lambda for given quant (and P)
            x <- data  # using data rather than x as a parameter is needed for newton.raphson()
            ## constraints for shape
            ## A %*% param + B > 0
            ## row 1: shape > 0
            A <- matrix(c(1), 1, 1, byrow=TRUE)
            B <- 0
            constraints <- list(ineqA=A, ineqB=B)
            out.q <- NA
            ## browser()
            ## maxLik
            out.s <- maxLik::maxLik(loglik.weibull.q.set,
                                    start = shape,
                                    x     = x,
                                    xcen  = xcen,
                                    quant = quant,
                                    P     = P,
                                    debug = debug,
                                    constraints = constraints,
                                    method = 'BFGS',
                                    iterlim = 2000)
            if (out.s$code == 0) {
                params.s      <- as.list(out.s$estimate)
                loglik.s      <- out.s$maximum
                convergence.s <- 'successful'
            } else {
                cat('####################################################################################\n')
                cat('WARNING: CONVERGENCE FAILURE IN mle.weibull when maximizing loglik.weibull.q.set    \n')
                cat('         for quant =', quant, 'and P =', P,                                        '\n')
                cat('####################################################################################\n')
                params.s      <- list(shape=NA)
                loglik.s      <- NA
                convergence.s <- out.smessage
            }


            ## ## optim
            ## out.q.optim <- stats::optim(par=c(shape=shape), # par holds optimization params
            ##                             loglik.weibull.q.set,
            ##                             x     = x,
            ##                             xcen  = xcen,
            ##                             quant = quant,
            ##                             P     = P,
            ##                             debug = debug)
            ## out.q.optim <- stats::optim(par=c(shape=shape), # par holds optimization params
            ##                             loglik.weibull.q.set,
            ##                             x     = x,
            ##                             xcen  = xcen,
            ##                             quant = quant,
            ##                             P     = P,
            ##                             debug = debug,
            ##                             method = 'Brent',
            ##                             lower = 1E-15,
            ##                             upper = max(x, xcen, na.rm=TRUE)*2)
            ## ## optimize
            ## out.q.optimize <- stats::optim(par=c(shape=shape), # par holds optimization params
            ##                                loglik.weibull.q.set,
            ##                                x     = x,
            ##                                xcen  = xcen,
            ##                                quant = quant,
            ##                                P     = P,
            ##                                debug = debug)
            ## out <- optimize(loglik.weibull.q.set.optim, interval=c(1E-15, max(x, xcen, na.rm=TRUE)),
            ##                 maximum = TRUE, x=x, xcen=xcen, quant=quant, P=P, debug=FALSE)

            if (isFALSE(debug)) {
                return(loglik.s)
            } else {
                return(list(out = out.s, loglik=loglik.s, params=params.s, convergence=convergence.s))
            }
        }


        ## test fit at quant.P (could comment out this block of code)
        cat('Attempting MLE fit on alternate parameters for P=', P, 'but with fixed quant=quant.P\n')
        out.s <- loglik.fixedq(quant.P.orig, x, xcen, P, shape.P, debug=TRUE)
        summary(out.s$out)
        temp <- data.frame(shape=out.s$params[[1]], scale=NA,
                           quant=quant.P.orig, P=P,
                           loglik=out.s$loglik, convergence=out.s$convergence,
                           optimizer    = 'maxLik',
                           max.function = 'loglik.weibull.q.set',
                           fit.params   = 'shape for given quant + P')
        out.all <- rbind(out.all, temp)
        cat('convergence for shape at quant.P:', out.s$convergence, '\n\n')
        
        
        ##----------------------
        if (isTRUE(plots)) {
            ## initial guess for 1st point
            shape.plot <- shape
            ## plot fit above quant.P
            xplot <- seq(quant.P, xmax, length.out=51)
            yplot <- NA
            for (ploti in 1:length(xplot)) {
                ## if (ploti == 3) browser()
                out <- loglik.fixedq(xplot[ploti], x, xcen, P, shape.plot, debug=TRUE)
                yplot[ploti] <- out$loglik
                points(xplot[ploti], yplot[ploti], col='black')
                if (yplot[ploti] < loglik.tol) break  # exit for loop
                ## better initial guess for next point
                ## shape.plot   <- out$params[[1]]
                ## cat(ploti, xplot[ploti], yplot[ploti], '\n')
            }
            xyplot.upper <- data.frame(quantile       = xplot[1:length(yplot)],
                                       log.likelihood = yplot,
                                       loglik.tol     = loglik.tol)
            ## plot fit below quant.P
            xplot <- seq(quant.P, xmin, length.out=51)
            yplot <- NA
            for (ploti in 1:length(xplot)) {
                ## if (ploti == 3) browser()
                out <- loglik.fixedq(xplot[ploti], x, xcen, P, shape.plot, debug=TRUE)
                yplot[ploti] <- out$loglik
                points(xplot[ploti], yplot[ploti], col='black')
                if (yplot[ploti] < loglik.tol) break  # exit for loop
            }
            xyplot.lower <- data.frame(quantile       = xplot[1:length(yplot)],
                                       log.likelihood = yplot,
                                       loglik.tol     = loglik.tol)
            xyplot <- rbind(xyplot.lower, xyplot.upper)
            points(xyplot$quantile, xyplot$log.likelihood)
            ## converged <- which(!is.na(yplot))
            ## points(xplot[converged], yplot[converged],
            ##        xlab='quantile', ylab='log likelihood', col='red')
            plotit <- function(xplot=0.994, col='blue', pch=16, cex=2) {
                ## this is just for using interactively if debugging
                yplot <- loglik.fixedq(xplot, x, xcen, P, shape.P, debug=FALSE)
                points(xplot, yplot, col=col, pch=16, cex=2)
                return(yplot)
            }
        }

        ## determine confidence bound as point where the likelihood ratio equals loglik.tol
        cat('Attempting MLE fit on alternate parameters for P=', P, 'to find quant where loglik corresponds to tolerance limit\n')
        quant.P.alpha.eff.l <- NA
        quant.P.alpha.eff.u <- NA
        if (side.which == 'lower' | (P < 0.5 & sided == 2)) {
            ## find lower tolerance limit
            out.nrl <- newton.raphson(f = loglik.fixedq,
                                      xguess = quant.P.alpha.eff.l.guess,
                                      ytarget = loglik.tol,
                                      data    = x,
                                      xcen = xcen,
                                      P      = P,
                                      shape  = shape,
                                      tol = 1e-5, 
                                      n = 1000,
                                      relax = 0.8,
                                      nrelax = 10,
                                      plots=plots.nr,
                                      plot.add=TRUE)
            quant.P.alpha.eff.l <- out.nrl$root
            if (is.null(quant.P.alpha.eff.l)) { quant.P.alpha.eff.l <- NA }
        } else {
            ## find upper tolerance limit
            out.nru <- newton.raphson(f = loglik.fixedq,
                                      xguess = quant.P.alpha.eff.u.guess,
                                      ytarget = loglik.tol,
                                      data   = x,
                                      xcen   = xcen,
                                      P      = P,
                                      shape  = shape,
                                      tol = 1e-5, 
                                      n = 1000,
                                      relax = 0.8,
                                      nrelax = 10,
                                      plots=plots.nr,
                                      plot.add=TRUE)
            quant.P.alpha.eff.u <- out.nru$root
            if (is.null(quant.P.alpha.eff.u)) { quant.P.alpha.eff.u <- NA }
        }
        cat('Final confidence interval for P=', P, '\n')
        cat(quant.P.alpha.eff.l, quant.P.alpha.eff.u, '\n\n')

        ## collect tolerance limit calculations (i.e., for 1-P and P)
        tol.limits <- c(tol.limits, quant.P.alpha.eff.l, quant.P.alpha.eff.u)
        
        if (isTRUE(plots)) {
            ## plot intersection with log likelihood curve
            if (!is.na(quant.P.alpha.eff.l)) points(quant.P.alpha.eff.l, loglik.tol, col='red', pch=16, cex=2)
            if (!is.na(quant.P.alpha.eff.u)) points(quant.P.alpha.eff.u, loglik.tol, col='red', pch=16, cex=2)
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
    tolerance <- data.frame(sided, alpha.eff=alpha.eff, conf, P, tol.lower, tol.upper)

    ## print final parameter comparison and tolerance limits
    print(as.data.frame(params.save))
    cat('\n')
    print(tolerance)
    cat('\n')
    
    return(list(loglik.plot    = xyplot,
                loglik.max     = loglik.max,
                loglik.tol     = loglik.tol,
                quant.P    = quant.P.save,
                params     = params.save,
                tol.approx = tol.approx,
                tolerance  = tolerance,
                out.all    = out.all))
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
    out.tol <- mle.weibull.tol(x, P=0.99, conf=0.95, plots=TRUE) # 1-sided, upper 99/95 = 
    tolerance::exttol.int(x)                                     # 1-sided, upper 99/95 = 2.0987

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
    out.lower <- mle.weibull.tol(x, param='auto',
                                 side.which='lower', sided=1, conf=0.99, P=0.01, 
                                 plots=TRUE, plots.nr=FALSE, debug=FALSE, main='lower, 1-sided 99/1')
    ## using default parameters except for 'plots'
    out.upper <- mle.weibull.tol(x, param='auto',
                                 side.which='upper', sided=1, conf=0.99, P=0.99, 
                                 plots=TRUE, plots.nr=FALSE, debug=FALSE, main='upper, 1-sided 99/99')
    ## lower and upper 2-sided tolerance limits
    ## test that 1-sided 99/99 is the same as a 2-sided 98/98
    out.twosided <- mle.weibull.tol(x, param='auto',
                                    side.which='both', sided=2, conf=0.98, P=0.98, 
                                    plots=TRUE, plots.nr=FALSE, debug=FALSE, main='2-sided 98/98')

    
    ## test that upper 1-sided 99/75 and lower 1-sided 99/25 are the same as a 2-sided 98/50  
    plotspace(2,2)
    ## lower and upper 1-sided tolerance limits
    out.lower <- mle.weibull.tol(x, param='auto',
                                 side.which='lower', sided=1, conf=0.99, P=0.25, 
                                 plots=TRUE, plots.nr=FALSE, debug=FALSE, main='lower, 1-sided 99/25')
    out.upper <- mle.weibull.tol(x, param='auto',
                                 side.which='upper', sided=1, conf=0.99, P=0.75, 
                                 plots=TRUE, plots.nr=FALSE, debug=FALSE, main='upper, 1-sided 99/75')
    ## lower and upper 2-sided tolerance limits
    out.twosided <- mle.weibull.tol(x, param='auto',
                                    side.which='both', sided=2, conf=0.98, P=0.5, 
                                    plots=TRUE, plots.nr=FALSE, debug=FALSE, main='2-sided 98/50')


    ##-------------------------------------------------------------------
    ## TEST OUT CENSORED DATA FUNCTIONALITY
    ##-------------------------------------------------------------------
    fit.compare.cen <- function(x, xcen, main=NULL) {
        ## plot histogram with no censored data and fit
        hist(x, freq=FALSE, border='black', main=main, xlim=c(2,5.4))
        out.fit0 <- mle.weibull.tol(x, xcen=NA, plots=FALSE)
        parms0 <- out.fit0$params
        curve(stats::dweibull(x, parms0$shape, parms0$scale), min(x), max(x), col='black', add=TRUE)
        abline(v=out.fit0$tolerance$tol.upper, col='black', lty=2)
        ## plot histogram with all data treated as known fit at numeric value
        xcen.avg <- rowMeans(xcen, na.rm=TRUE) # use the average for interval data
        x.all <- c(x, xcen.avg)
        hist(x.all, freq=FALSE, border='red', add=TRUE)
        out.fit1 <- mle.weibull.tol(x.all, xcen=NA, plots=FALSE)
        parms1 <- out.fit1$params
        curve(stats::dweibull(x, parms1$shape, parms1$scale), min(x), max(x), col='red', add=TRUE)
        abline(v=out.fit1$tolerance$tol.upper, col='red', lty=2)
        ## plot fit if treat xcen as censored
        out.fit2 <- mle.weibull.tol(x, xcen=xcen, plots=FALSE)
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
    print(fit1a)
    ## it may be that the fit was bad because the following does pull the 99/99 in as I would expect
    ## out.qq <- qqplot_nwj(x, type='j', mainadder = 'x only')
    ## x.all <- na.omit(c(x, xcen$x.low, xcen$x.high))
    ## out.qq <- qqplot_nwj(x.all, type='j', mainadder='all')
    ## out.qq <- qqplot_nwj(x, xcen, type='j', mainadder='censored')
    x.all <- na.omit(c(x, xcen$x.low, xcen$x.high))
    out.qq <- qqplot_nwj(x.all, type='j', mainadder='all')
    
    x.low  <- NA
    x.high <- 2.14
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit1b <- fit.compare.cen(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    x.all <- na.omit(c(x, xcen$x.low, xcen$x.high))
    out.qq <- qqplot_nwj(x.all, type='j', mainadder='all')
    
    x.low  <- NA
    x.high <- 2.15
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit1c <- fit.compare.cen(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))

    x.low  <- NA
    x.high <- 2.2
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit1d <- fit.compare.cen(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))

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

