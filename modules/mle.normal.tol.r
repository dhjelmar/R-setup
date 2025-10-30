loglik.normal.q <- function(x=NA, xcen=NA, param=c(quant, sdev), P, debug=FALSE){
    ## calculate log likelihhod for distribution and given parameters
    quant <- param[[1]]  # replaced xbar with quant as a parameter
    sdev  <- param[[2]]
    xbar  <- quant - sdev * qnorm(P)
    loglik.normal(x, xcen, param=c(xbar, sdev), debug=FALSE)
}

loglik.normal.q.set <- function(x=NA, xcen=NA, param=sdev, quant, P, debug=FALSE){
    ## calculate log likelihhod for distribution and given parameters
    sdev <- param[[1]]
    xbar <- quant - sdev * qnorm(P)
    loglik.normal(x, xcen, param=c(xbar, sdev), debug=FALSE)
}

mle.normal.tol <- function(x, xcen=NA, param='auto',
                           side.which='upper', sided=1, conf=0.99, alpha=NULL, P=0.99,
                           plots=FALSE, plots.nr=FALSE, debug=FALSE, main=NULL, old=TRUE) {
    
    ## normal distribution
    ## MLE (Maximum Likelihood Estimate) fit to determine parameters
    ## LR (Likelihood Ratio) appraoch to find tolerance limit
    ## e.g., 95% upper confidence bound on the 99th percentile
    ## wikipedia     convention: "100×p%    / 100×(1−α) tolerance interval"
    ## more standard convention: "100×(1−α) / 100×p%    tolerance interval"

    ## input: x     = vector of known data
    ##        xcen  = dataframe of censored data
    ##                (1st column = low value or NA; 2nd column = high value or NA)
    ##        param = initial guess for fit parameters for: sdev and xbar
    ##              = 'auto' (default) uses mean() and sd() for initial guess of parameters
    ##                and will try using list(sdev=1, xbar=1) if that fails
    ##        side.which = 'upper' or 'lower' (not used if sided = 2)
    ##        sided = 1 (default) means 1-sided tolerance limit
    ##              = 2 means 2-sided tolerance limit
    ##        conf  = 1 - alpha = confidence used to determine chi-square
    ##
    ##        old:
    ##            alpha.eff = NULL (default) sets effective alpha, alpha.eff = 2*(1-conf)/sided, for use with qchisq();
    ##
    ##        new:
    ##            alpha = NULL (default) sets effective alpha, alpha = 1-conf, for use with qchisq(1-2*(alpha/sided),1);
    ##
    ##                    e.g., if want 1-sided conf = 99%:
    ##                          sided = 1
    ##                          alpha = (1-0.99) = 0.01
    ##                          chi-square = qchisq(1-2*alpha/sided,1) = 5.411894
    ##
    ##                    e.g., if want 2-sided conf = 99%
    ##                          sided = 2
    ##                          alpha = 1-0.01 = 0.01
    ##                          chi-square = qchisq(1-2*alpha/sided,1) = 6.634897
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

    if (is.null(alpha)) {
        ## set alpha level for chi-square(1-2*alpha,1) for use in confidence limit calculation
        alpha <- (1-conf)*sided
    } else {
        ## calculate confidence limit from alpha
        conf <- 1 - alpha / sided
    }
    
    ## if (isTRUE(plots) & sided == 2) par(mfrow=c(1,2))
    out <- NULL

    ##-----------------------------------------------------------------------------
    ## Normal SU parameters
    out <- mle.normal(x, xcen, param=param, plots=FALSE)
    params <- out$parms
    sdev  <- params$sdev
    xbar  <- params$xbar
    params.save <- params
    params$type <- NULL
    out.all  <- as.data.frame(params)
    out.all$quant <- NA
    out.all$P     <- NA
    out.all$loglik <- out$loglik
    out.all$convergence <- out$convergence
    out.all$optimizer    <- 'maxLik'
    out.all$max.function <- 'loglik.normal'
    out.all$fit.params   <- 'xbar, sdev'
    cat('convergence for best estimate parameters:', out.all$convergence, '\n')
    cat('\n')
   
    ##-----------------------------------------------------------------------------
    ## find confidence limit at level alpha for requested coverage, P
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
        quant.P.orig <- stats::qnorm(P, mean=params$xbar, sd=params$sdev)
        quant.P.save[k] <- quant.P.orig
        quant.param <- list(quant=quant.P.orig, sdev=params$sdev)

        ##----------------------
        ## calculate confidence limits using LR (Likelihood Ratio)
        ## confidene limit is defined at likelihood that is lower than max by chi-squared
        loglik.max <- loglik.normal.q(x, xcen, quant.param, P)
        if (old) {
            alpha.eff <- 2*(1-conf)/sided
            loglik.tol <- loglik.max - qchisq(1 - alpha.eff, 1)/2   # qchisq(1-0.02, 1) = 5.411894
        } else {
            loglik.tol <- loglik.max - qchisq(1 - 2*(alpha/sided), 1)/2   # qchisq(1-0.02, 1) = 5.411894
        }
        cat('P =', P, 'MLE =', loglik.max, '; tolerance limit at MLE =', loglik.tol, '\n')
        temp <- data.frame(sdev=sdev, xbar=xbar,
                           quant=quant.P.orig, P=P, loglik=loglik.max, convergence=NA,
                           optimizer    = NA,
                           max.function = NA,
                           fit.params   = 'quant calculated for P')
        out.all <- rbind(out.all, temp)
        cat('\n')

        ##----------------------
        ## refit on alternate parameters to determine standard error for fit on quantile
        cat('Attempting MLE fit on alternate parameters for P=', P, '\n')
        ## constraints for sdev and/or xbar
        ## A %*% param + B > 0
        ## row 1: sdev > 0
        A <- matrix(c(0,1), 1, 2, byrow=TRUE)
        B <- 0
        constraints <- list(ineqA=A, ineqB=B)
        out.qs <- NA
        out.qs <- maxLik::maxLik(loglik.normal.q,
                                   start = unlist(quant.param),  # quant, sdev
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
            cat('WARNING: CONVERGENCE FAILURE IN mle.normal.tol when maximizing loglik.normal.q\n')
            cat('####################################################################################\n')
        }
        loglik.max.qs <- out.qs$maximum
        params.qs <- as.list(out.qs$estimate)
        quant.P  <- params.qs[[1]]
        sdev.P  <- params.qs[[2]]
        xbar.P  <- quant.P/(-log((1-P)))^(1/sdev.P)
        params.qs$xbar <- xbar.P
        ## params.q.save[k] <- list(params.qs)
        print(as.data.frame(params.qs))
        ## following does not work for maxLik because hessian often has negative diagonals
        ## standard.error <- as.numeric( sqrt(diag(solve(out.qs$hessian))) )
        standard.error <- summary(out.qs)$estimate['quant', 'Std. error']
        temp <- data.frame(sdev=sdev.P, xbar=xbar.P,
                           quant=quant.P, P=P, loglik=loglik.max.qs, convergence=convergence.qs,
                           optimizer    = 'maxLik',
                           max.function = 'loglik.normal.q',
                           fit.params   = 'quant, sdev')
        out.all <- rbind(out.all, temp)
        cat('convergence for quant, sdev:', convergence.qs, '\n')
        cat('\n')
        
        ##----------------------
        ## estimate confidence limit using standard error to serve as starting point for search
        ## estimates      <- as.numeric( out.qs$par )
        ## coef           <- data.frame(estimates, standard.error)
        ## rownames(coef) <- names(out.qs$par)
        dof    <- length(x.avg) - 2    # 2 independent fitting parameters in Normal
        student.t <- qt(1 - (1-conf)/sided, dof)  # 2.3326 for dof=598
        if (is.nothing(standard.error) | standard.error == Inf) standard.error <- 0.001 * quant.P
        quant.P.alpha.l.guess <- quant.P - student.t * standard.error
        quant.P.alpha.u.guess <- quant.P + student.t * standard.error
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
            plot(quant.P.save[k], loglik.max, col='blue', pch=16, cex=2,
                 xlab='quantile', ylab='log likelihood',
                 xlim=range(xmin, xmax),
                 ylim=range(loglik.max, loglik.tol),
                 main=main)
            abline(h=loglik.tol)
            abline(v=c(quant.P.alpha.l.guess, quant.P.alpha.u.guess), col='red', lty=2)
        }
        
        ##----------------------
        ## function to fit on sdev
        loglik.fixedq <- function(quant, data, xcen, P, sdev, debug=FALSE) {
            ## first determine best fit delta, xi, and lambda for given quant (and P)
            x <- data  # using data rather than x as a parameter is needed for newton.raphson()
            ## constraints for sdev
            ## A %*% param + B > 0
            ## row 1: sdev > 0
            A <- matrix(c(1), 1, 1, byrow=TRUE)
            B <- 0
            constraints <- list(ineqA=A, ineqB=B)
            out.q <- NA
            ## browser()
            ## maxLik
            out.s <- maxLik::maxLik(loglik.normal.q.set,
                                    start = sdev,
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
                cat('WARNING: CONVERGENCE FAILURE IN mle.normal when maximizing loglik.normal.q.set    \n')
                cat('         for quant =', quant, 'and P =', P,                                        '\n')
                cat('####################################################################################\n')
                params.s      <- list(sdev=NA)
                loglik.s      <- NA
                convergence.s <- out.smessage
            }

            if (isFALSE(debug)) {
                return(loglik.s)
            } else {
                return(list(out = out.s, loglik=loglik.s, params=params.s, convergence=convergence.s))
            }
        }


        ## test fit at quant.P (could comment out this block of code)
        cat('Attempting MLE fit on alternate parameters for P=', P, 'but with fixed quant=quant.P\n')
        out.s <- loglik.fixedq(quant.P.orig, x, xcen, P, sdev.P, debug=TRUE)
        summary(out.s$out)
        temp <- data.frame(sdev=out.s$params[[1]], xbar=NA,
                           quant=quant.P.orig, P=P,
                           loglik=out.s$loglik, convergence=out.s$convergence,
                           optimizer    = 'maxLik',
                           max.function = 'loglik.normal.q.set',
                           fit.params   = 'sdev for given quant + P')
        out.all <- rbind(out.all, temp)
        cat('convergence for sdev at quant.P:', out.s$convergence, '\n\n')
        
        
        ##----------------------
        if (isTRUE(plots)) {
            ## initial guess for 1st point
            sdev.plot <- sdev
            ## plot fit above quant.P
            xplot <- seq(quant.P, xmax, length.out=51)
            yplot <- NA
            for (ploti in 1:length(xplot)) {
                ## if (ploti == 3) browser()
                out <- loglik.fixedq(xplot[ploti], x, xcen, P, sdev.plot, debug=TRUE)
                yplot[ploti] <- out$loglik
                points(xplot[ploti], yplot[ploti], col='black')
                if (yplot[ploti] < loglik.tol) break  # exit for loop
                ## better initial guess for next point
                ## sdev.plot   <- out$params[[1]]
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
                out <- loglik.fixedq(xplot[ploti], x, xcen, P, sdev.plot, debug=TRUE)
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
                yplot <- loglik.fixedq(xplot, x, xcen, P, sdev.P, debug=FALSE)
                points(xplot, yplot, col=col, pch=16, cex=2)
                return(yplot)
            }
        }

        ## determine confidence bound as point where the likelihood ratio equals loglik.tol
        cat('Attempting MLE fit on alternate parameters for P=', P, 'to find quant where loglik corresponds to tolerance limit\n')
        quant.P.alpha.l <- NA
        quant.P.alpha.u <- NA
        if (side.which == 'lower' | (P < 0.5 & sided == 2)) {
            ## find lower tolerance limit
            out.nrl <- newton.raphson(f = loglik.fixedq,
                                      xguess = quant.P.alpha.l.guess,
                                      ytarget = loglik.tol,
                                      data    = x,
                                      xcen = xcen,
                                      P      = P,
                                      sdev  = sdev,
                                      tol = 1e-5, 
                                      n = 1000,
                                      relax = 0.8,
                                      nrelax = 10,
                                      plots=plots.nr,
                                      plot.add=TRUE)
            quant.P.alpha.l <- out.nrl$root
            if (is.null(quant.P.alpha.l)) { quant.P.alpha.l <- NA }
        } else {
            ## find upper tolerance limit
            out.nru <- newton.raphson(f = loglik.fixedq,
                                      xguess = quant.P.alpha.u.guess,
                                      ytarget = loglik.tol,
                                      data   = x,
                                      xcen   = xcen,
                                      P      = P,
                                      sdev  = sdev,
                                      tol = 1e-5, 
                                      n = 1000,
                                      relax = 0.8,
                                      nrelax = 10,
                                      plots=plots.nr,
                                      plot.add=TRUE)
            quant.P.alpha.u <- out.nru$root
            if (is.null(quant.P.alpha.u)) { quant.P.alpha.u <- NA }
        }
        cat('Final confidence interval for P=', P, '\n')
        cat(quant.P.alpha.l, quant.P.alpha.u, '\n\n')

        ## collect tolerance limit calculations (i.e., for 1-P and P)
        tol.limits <- c(tol.limits, quant.P.alpha.l, quant.P.alpha.u)
        
        if (isTRUE(plots)) {
            ## plot intersection with log likelihood curve
            if (!is.na(quant.P.alpha.l)) points(quant.P.alpha.l, loglik.tol, col='red', pch=16, cex=2)
            if (!is.na(quant.P.alpha.u)) points(quant.P.alpha.u, loglik.tol, col='red', pch=16, cex=2)
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

    ## collect tolerance values in dataframe similar to extol.int for normal
    tolerance <- data.frame(sided, alpha=alpha, conf, P, tol.lower, tol.upper)

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


mle.normal.tol.test <- function() {
    source('setup.r')

    ##-------------------------------------------------------------------
    ## create Normal dataset
    set.seed(1)
    parms <- list(xbar=10, sdev=1)
    parms <- list(xbar=10, sdev=2)
    x <- stats::rnorm(1000, mean=parms$xbar, sd=parms$sdev)

    ## determine normal parameters
    plotspace(1,2)
    out.fit <- mle.normal(x, plots=TRUE)
    parms.mle <- out.fit$parms
    print(out.fit$parms.compare)    # parameters and loglik are nearly identical

    ## normal tolerance limits from my module1
    out.tol <- mle.normal.tol(x, P=0.99, conf=0.95, plots=TRUE) # 1-sided, upper 99/95 = 15.00248
    
    ## normal tolerance limits from tolerance package
    out.tol.answer <- tolerance::normtol.int(x,P=0.99, alpha=0.1, side=1)
    out.tol.answer$`1-sided.upper`                              # 1-sided, upper 99/95 = 14.95875

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
    out.lower <- mle.normal.tol(x, param='auto',
                                 side.which='lower', sided=1, conf=0.99, P=0.01, 
                                 plots=TRUE, plots.nr=FALSE, debug=FALSE, main='lower, 1-sided 99/1')
    ## using default parameters except for 'plots'
    out.upper <- mle.normal.tol(x, param='auto',
                                 side.which='upper', sided=1, conf=0.99, P=0.99, 
                                 plots=TRUE, plots.nr=FALSE, debug=FALSE, main='upper, 1-sided 99/99')
    ## lower and upper 2-sided tolerance limits
    ## test that 1-sided 99/99 is the same as a 2-sided 98/98
    out.twosided <- mle.normal.tol(x, param='auto',
                                    side.which='both', sided=2, conf=0.98, P=0.98, 
                                    plots=TRUE, plots.nr=FALSE, debug=FALSE, main='2-sided 98/98')

    
    ## test that upper 1-sided 99/75 and lower 1-sided 99/25 are the same as a 2-sided 98/50  
    plotspace(2,2)
    ## lower and upper 1-sided tolerance limits
    out.lower <- mle.normal.tol(x, param='auto',
                                 side.which='lower', sided=1, conf=0.99, P=0.25, 
                                 plots=TRUE, plots.nr=FALSE, debug=FALSE, main='lower, 1-sided 99/25')
    out.upper <- mle.normal.tol(x, param='auto',
                                 side.which='upper', sided=1, conf=0.99, P=0.75, 
                                 plots=TRUE, plots.nr=FALSE, debug=FALSE, main='upper, 1-sided 99/75')
    ## lower and upper 2-sided tolerance limits
    out.twosided <- mle.normal.tol(x, param='auto',
                                    side.which='both', sided=2, conf=0.98, P=0.5, 
                                    plots=TRUE, plots.nr=FALSE, debug=FALSE, main='2-sided 98/50')


    ##-------------------------------------------------------------------
    ## TEST OUT CENSORED DATA FUNCTIONALITY
    ##-------------------------------------------------------------------
    fit.compare.cen <- function(x, xcen, main=NULL) {
        ## plot histogram with no censored data and fit
        hist(x, freq=FALSE, border='black', main=main, xlim=c(2,5.4))
        out.fit0 <- mle.normal.tol(x, xcen=NA, plots=FALSE)
        parms0 <- out.fit0$params
        curve(stats::dnorm(x, sd=parms0$sdev, mean=parms0$xbar), min(x), max(x), col='black', add=TRUE)
        abline(v=out.fit0$tolerance$tol.upper, col='black', lty=2)
        ## plot histogram with all data treated as known fit at numeric value
        xcen.avg <- rowMeans(xcen, na.rm=TRUE) # use the average for interval data
        x.all <- c(x, xcen.avg)
        hist(x.all, freq=FALSE, border='red', add=TRUE)
        out.fit1 <- mle.normal.tol(x.all, xcen=NA, plots=FALSE)
        parms1 <- out.fit1$params
        curve(stats::dnorm(x, sd=parms1$sdev, mean=parms1$xbar), min(x), max(x), col='red', add=TRUE)
        abline(v=out.fit1$tolerance$tol.upper, col='red', lty=2)
        ## plot fit if treat xcen as censored
        out.fit2 <- mle.normal.tol(x, xcen=xcen, plots=FALSE)
        parms2 <- out.fit2$params
        curve(stats::dnorm(x, sd=parms2$sdev, mean=parms2$xbar), min(x), max(x), col='blue', type='p', add=TRUE)
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

