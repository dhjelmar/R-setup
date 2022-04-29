loglik.johnsonsu.q <- function(x=NA, xcen=NA, param=c(quant, delta, xi, lambda), P, debug=FALSE){
    ## calculate log likelihhod for distribution and given parameters
    quant  <- param[[1]]  # replaced gamma with quant as a parameter
    delta  <- param[[2]]
    xi     <- param[[3]]
    lambda <- param[[4]]
    ## write gamma as a function of quant, delta, xi and lambda
    gamma <- qnorm(P) - delta * asinh( (quant-xi)/lambda )
    loglik.johnsonsu(x, xcen, param=c(gamma, delta, xi, lambda), debug=FALSE)
}

loglik.johnsonsu.q.set <- function(x=NA, xcen=NA, param=c(delta, xi, lambda), quant, P, debug=FALSE){
    ## calculate log likelihhod for distribution and given parameters
    delta  <- param[[1]]
    xi     <- param[[2]]
    lambda <- param[[3]]
    ## write gamma as a function of quant, delta, xi and lambda
    gamma <- qnorm(P) - delta * asinh( (quant-xi)/lambda )
    loglik.johnsonsu(x, xcen, param=c(gamma, delta, xi, lambda), debug=FALSE)
}

mle.johnsonsu.tol <- function(x, xcen=NA, param='auto',
                              side.which='upper', sided=1, P=0.99, conf=0.99, alpha.eff=NULL,
                              plots=FALSE, plots.nr=FALSE, debug=FALSE, main=NULL) {
    
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
    
    if (is.data.frame(x)) x <- x[1] # convert to vector

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
    ## Johnson SU parameters
    out <- mle.johnsonsu(x, xcen, param=param, plots=FALSE)
    params <- out$parms
    gamma  <- params$gamma
    delta  <- params$delta
    xi     <- params$xi
    lambda <- params$lambda
    params.save <- params
    params$type <- NULL
    out.all  <- as.data.frame(params)
    out.all$quant <- NA
    out.all$P     <- NA
    out.all$loglik <- out$loglik
    out.all$convergence <- out$convergence
    out.all$optimizer    <- 'maxLik'
    out.all$max.function <- 'loglik.johnsonsu'
    out.all$fit.params   <- 'gamma, delta, xi, lambda'
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
        ## following is identical to
        ## quant.P.orig <- ExtDist::qJohnsonSU(P, params = params)
        quant.P.orig  <- xi + lambda * sinh( (qnorm(P)-gamma)/delta )
        quant.P.save[k] <- quant.P.orig
        quant.param <- list(quant=quant.P.orig, delta=params$delta, xi=params$xi, lambda=params$lambda)
        
        ##----------------------
        ## calculate confidence limits using LR (Likelihood Ratio)
        ## confidene limit is defined at likelihood that is lower than max by chi-squared
        loglik.max <- loglik.johnsonsu.q(x, xcen, quant.param, P)
        loglik.tol <- loglik.max - qchisq(1 - alpha.eff, 1)/2   # qchisq(1-0.02, 1) = 5.411894
        cat('P =', P, 'MLE =', loglik.max, '; tolerance limit at MLE =', loglik.tol, '\n')
        temp <- data.frame(gamma=gamma, delta=delta, xi=xi, lambda=lambda,
                           quant=quant.P.orig, P=P, loglik=loglik.max, convergence=NA,
                           optimizer    = NA,
                           max.function = NA,
                           fit.params   = 'quant calculated for P')
        out.all <- rbind(out.all, temp)
        cat('\n')

        ##----------------------
        ## refit on alternate parameters to determine standard error for fit on quantile
        cat('Attempting MLE fit on alternate parameters for P=', P, '\n')
        ## constraints for quant, delta, xi, and/or lambda
        ## A %*% param + B > 0
        ## row 1: delta  > 0
        ## row 2: lambda > 0
        A <- matrix(c(0,1,0,0,  0,0,0,1), 2, 4, byrow=TRUE)
        B <- matrix(c(0,0),               2, 1)
        constraints <- list(ineqA=A, ineqB=B)
        out.qdxl <- NA
        out.qdxl <- maxLik::maxLik(loglik.johnsonsu.q,
                                        start = unlist(quant.param),  # quant, delta, xi, lambda
                                        x     = x,
                                        xcen  = xcen,
                                        P     = P,
                                        debug = debug,
                                        constraints = constraints,
                                        iterlim = 2000)
        print(summary(out.qdxl))
        convergence.qdxl <- if (out.qdxl$message == 'successful convergence ') {'successful'}
                            else {out.qdxl$message}
        if (convergence.qdxl != 'successful') {
            cat('####################################################################################\n')
            cat('WARNING: CONVERGENCE FAILURE IN mle.johnsonsu.tol when maximizing loglik.johnsonsu.q\n')
            cat('####################################################################################\n')
        }
        loglik.max.qdxl <- out.qdxl$maximum
        params.qdxl <- as.list(out.qdxl$estimate)
        quant.P  <- params.qdxl[[1]]
        delta.P  <- params.qdxl[[2]]
        xi.P     <- params.qdxl[[3]]
        lambda.P <- params.qdxl[[4]]
        gamma.P <- qnorm(P) - delta.P * asinh( (quant.P-xi.P)/lambda.P )
        params.qdxl$gamma <- gamma.P
        ## params.q.save[k] <- list(params.qdxl)
        print(as.data.frame(params.qdxl))
        ## following does not work for maxLik because hessian often has negative diagonals
        ## standard.error <- as.numeric( sqrt(diag(solve(out.qdxl$hessian))) )
        standard.error <- summary(out.qdxl)$estimate['quant', 'Std. error']
        temp <- data.frame(gamma=gamma.P, delta=delta.P, xi=xi.P, lambda=lambda.P,
                           quant=quant.P, P=P, loglik=loglik.max.qdxl, convergence=convergence.qdxl,
                           optimizer    = 'maxLik',
                           max.function = 'loglik.johnsonsu.q',
                           fit.params   = 'quant, delta, xi, lambda')
        out.all <- rbind(out.all, temp)
        cat('convergence for quant, delta, xi, lambda:', convergence.qdxl, '\n')
        cat('\n')
        
        ##----------------------
        ## estimate confidence limit using standard error to serve as starting point for search
        ## estimates      <- as.numeric( out.qdxl$par )
        ## coef           <- data.frame(estimates, standard.error)
        ## rownames(coef) <- names(out.qdxl$par)
        dof    <- length(x) - 4    # 4 independent fitting parameters in Johnson SU
        student.t <- qt(1 - (1-conf)/sided, dof)  # 2.3326 for dof=598
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
        ## function to fit on delta, xi and lambda
        loglik.fixedq <- function(quant, data, xcen, P, delta, xi, lambda, debug=FALSE) {
            ## first determine best fit delta, xi, and lambda for given quant (and P)
            x <- data  # using data rather than x as a parameter is needed for newton.raphson()
            ## constraints for delta, xi, and/or lambda
            ## A %*% param + B > 0
            ## row 1: delta  > 0
            ## row 2: lambda > 0
            A <- matrix(c(1,0,0,  0,0,1), 2, 3, byrow=TRUE)
            B <- matrix(c(0,0),              2, 1)
            constraints <- list(ineqA=A, ineqB=B)
            out.dxl <- NA
            out.dxl <- maxLik::maxLik(loglik.johnsonsu.q.set,
                                      start = c(delta, xi, lambda),  # delta, xi, lambda
                                      x     = x,
                                      xcen  = xcen,
                                      quant = quant,
                                      P     = P,
                                      debug = debug,
                                      constraints = constraints,
                                      iterlim = 2000)
            ## print(summary(out.dxl))
            if (out.dxl$message == 'successful convergence ') {
                params.dxl      <- as.list(out.dxl$estimate)
                loglik.dxl      <- out.dxl$maximum
                convergence.dxl <- 'successful'
            } else {
                cat('####################################################################################\n')
                cat('WARNING: CONVERGENCE FAILURE IN mle.johnsonsu when maximizing loglik.johnsonsu.q.set\n')
                cat('         for quant =', quant, 'and P =', P,                                        '\n')
                cat('####################################################################################\n')
                params.dxl      <- list(delta=NA, xi=NA, lambda=NA)
                loglik.dxl      <- NA
                convergence.dxl <- out.dxl$message
            }
            if (isFALSE(debug)) {
                return(loglik.dxl)
            } else {
                return(list(out    = out.dxl,
                            loglik = loglik.dxl,
                            params = params.dxl,
                            convergence = convergence.dxl))
            }
        }


        ## test fit at quant.P (could comment out this block of code)
        cat('Attempting MLE fit on alternate parameters for P=', P, 'but with fixed quant=quant.P\n')
        out.dxl <- loglik.fixedq(quant.P.orig, x, xcen, P, delta.P, xi.P, lambda.P, debug=TRUE)
        summary(out.dxl$out)
        temp <- data.frame(gamma=NA, delta=out.dxl$params[[1]], xi=out.dxl$params[[2]],
                           lambda=out.dxl$params[[3]], quant=quant.P.orig, P=P,
                           loglik=out.dxl$loglik, convergence=out.dxl$convergence,
                           optimizer    = 'maxLik',
                           max.function = 'loglik.johnsonsu.q.set',
                           fit.params   = 'delta, xi, lambda for given quant + P')
        out.all <- rbind(out.all, temp)
        cat('convergence for delta, xi, lambda at quant.P:', out.dxl$convergence, '\n\n')
        
        
        ##----------------------
        if (isTRUE(plots)) {
            ## initial guess for 1st point
            delta.plot <- delta
            xi.plot <- xi
            lambda.plot <- lambda
            ## plot fit above quant.P
            xplot <- seq(quant.P, xmax, length.out=51)
            yplot <- NA
            for (ploti in 1:length(xplot)) {
                ## if (ploti == 3) browser()
                out <- loglik.fixedq(xplot[ploti], x, xcen, P, delta.plot, xi.plot, lambda.plot, debug=TRUE)
                yplot[ploti] <- out$loglik
                points(xplot[ploti], yplot[ploti], col='black')
                if (yplot[ploti] < loglik.tol) break  # exit for loop
                ## better initial guess for next point
                ## delta.plot   <- out$params[[1]]
                ## xi.plot      <- out$params[[2]]
                ## lambda.plot  <- out$params[[3]]
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
                out <- loglik.fixedq(xplot[ploti], x, xcen, P, delta.plot, xi.plot, lambda.plot, debug=TRUE)
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
                yplot <- loglik.fixedq(xplot, x, xcen, P, delta.P, xi.P, lambda.P, debug=FALSE)
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
                                      delta  = delta,  # initial guess for loglik.fixedq
                                      xi     = xi,     # initial guess for loglik.fixedq
                                      lambda = lambda, # initial guess for loglik.fixedq
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
                                      delta  = delta,  # initial guess for loglik.fixedq
                                      xi     = xi,     # initial guess for loglik.fixedq
                                      lambda = lambda, # initial guess for loglik.fixedq
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
    tolerance <- data.frame(sided, alpha.eff=alpha.eff, P, conf, tol.lower, tol.upper)

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


mle.johnsonsu.tol.test <- function() {
    source('setup.r')

    ##-------------------------------------------------------------------
    ## create Johnson SU dataset
    set.seed(1)
    jparms <- list(gamma=-3.3, delta=5.3, xi=1.8, lambda=1.9, type='SU')
    x <- ExtDist::rJohnsonSU(1000, param=jparms)

    ## determine johnson su parameters
    plotspace(2,2)
    out.fit <- mle.johnsonsu(x, plots=TRUE)
    jparms.mle <- out.fit$parms
    print(out.fit$parms.compare)

    ## johnson su tolerance limits
    out.tol <- mle.johnsonsu.tol(x, plots=TRUE)

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
    ## wikipedia order: "100×p%/100×(1−α) tolerance interval"
    out.lower <- mle.johnsonsu.tol(data=x, param='auto',
                                   side.which='lower', sided=1, P=0.01,  conf=0.99,
                                   plots=TRUE, plots.nr=FALSE, debug=FALSE, main='lower, 1-sided 1/99')
    ## using default parameters except for 'plots'
    out.upper <- mle.johnsonsu.tol(data=x, param='auto',
                                   side.which='upper', sided=1, P=0.99,  conf=0.99,
                                   plots=TRUE, plots.nr=FALSE, debug=FALSE, main='upper, 1-sided 99/99')
    ## lower and upper 2-sided tolerance limits
    ## test that 1-sided 99/99 is the same as a 2-sided 98/98
    out.twosided <- mle.johnsonsu.tol(data=x, param='auto',
                                      side.which='both', sided=2, P=0.98,  conf=0.98,
                                      plots=TRUE, plots.nr=FALSE, debug=FALSE, main='2-sided 98/98')

    
    ## test that upper and lower 1-sided 75/99 are the same as a 2-sided 50/98  
    plotspace(2,2)
    ## lower and upper 1-sided tolerance limits
    out.lower <- mle.johnsonsu.tol(data=x, param='auto',
                                   side.which='lower', sided=1, conf=0.99, P=0.25, 
                                   plots=TRUE, plots.nr=FALSE, debug=FALSE, main='lower, 1-sided 99/25')
    out.upper <- mle.johnsonsu.tol(data=x, param='auto',
                                   side.which='upper', sided=1, conf=0.99, P=0.75, 
                                   plots=TRUE, plots.nr=FALSE, debug=FALSE, main='upper, 1-sided 99/75')
    ## lower and upper 2-sided tolerance limits
    out.twosided <- mle.johnsonsu.tol(data=x, param='auto',
                                      side.which='both', sided=2, conf=0.98, P=0.5, 
                                      plots=TRUE, plots.nr=FALSE, debug=FALSE, main='2-sided 98/50')


    ##-------------------------------------------------------------------
    ## TEST OUT CENSORED DATA FUNCTIONALITY
    ##-------------------------------------------------------------------
    fit.compare.cen <- function(x, xcen, main=NULL) {
        ## plot histogram with no censored data and fit
        hist(x, freq=FALSE, border='black', main=main, xlim=c(2,5), ylim=c(0,1))
        out.fit0 <- mle.johnsonsu.tol(x, xcen=NA, plots=FALSE)
        jparms0 <- out.fit0$params
        curve(ExtDist::dJohnsonSU(x, params = jparms0), min(x), max(x), col='black', add=TRUE)
        abline(v=out.fit0$tolerance$tol.upper, col='black', lty=2)
        ## plot histogram with all data treated as known fit at numeric value
        xcen.avg <- rowMeans(xcen, na.rm=TRUE) # use the average for interval data
        x.all <- c(x, xcen.avg)
        hist(x.all, freq=FALSE, border='red', add=TRUE)
        out.fit1 <- mle.johnsonsu.tol(x.all, xcen=NA, plots=FALSE)
        jparms1 <- out.fit1$params
        curve(ExtDist::dJohnsonSU(x, params = jparms1), min(x), max(x), col='red', add=TRUE)
        abline(v=out.fit1$tolerance$tol.upper, col='red', lty=2)
        ## plot fit if treat xcen as censored
        out.fit2 <- mle.johnsonsu.tol(x, xcen=xcen, plots=FALSE)
        jparms2 <- out.fit2$params
        curve(ExtDist::dJohnsonSU(x, params = jparms2), min(x), max(x), col='blue', type='p', add=TRUE)
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

    ## create data set
    set.seed(1)
    jparms <- list(gamma=-3.3, delta=5.3, xi=1.8, lambda=1.9, type='SU')
    x <- ExtDist::rJohnsonSU(1000, param=jparms)
    plotspace(1,1)
    
    ## left censored data should reduce the upper 99/99
    ## my expectation did not happen
    xnum <- 50
    x.low  <- NA
    x.high <- 2.0
    xcen1 <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit1a <- fit.compare.cen(x, xcen1, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    print(fit1a)

    plotspace(1,1)
    ## data that are high should push out the upper 99/99
    ## following shows that is true (99/99 for fit with xcen is >> that for fit to x)
    x.low  <- 3.7
    x.high <- NA
    xcen2 <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit2 <- fit.compare.cen(x, xcen2, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    print(fit2)
    
    ## data that are unknown over the middle range of data should strengthen confidence in 99/99
    ## following shows that is true (99/99 for fit with xcen is < that for fit to x)
    x.low  <- 2.6
    x.high <- 3.5
    xcen3 <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit3 <- fit.compare.cen(x, xcen3, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    print(fit3)
    
    ## data that are unknown between almost the lowest value and Inf should have no real impact
    ## following shows that is true (99/99 for fit with xcen is about the same as that for fit to x)
    x.low  <- 2.2
    x.high <- NA
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit4 <- fit.compare.cen(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    print(fit4)
   
    ## many low censored points
    xnum <-length(x)
    x.low  <- NA
    x.high <- 2.2
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit5 <- fit.compare.cen(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    print(fit5)
    
    ## multiple censored data
    xnum <- 2
    x.low  <- NA
    x.high <- 2.0
    xcen1 <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    x.low  <- 3.7
    x.high <- NA
    xcen2 <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    x.low  <- 2.6
    x.high <- 3.5
    xcen3 <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    xcen <- rbind(xcen1, xcen2, xcen3)
    
    out.tol <- mle.johnsonsu.tol(x, xcen, plots=TRUE)
    out.qq <- qqplot_censored(x, xcen)
    
}

