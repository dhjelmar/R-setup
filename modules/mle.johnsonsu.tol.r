mle.johnsonsu.tol <- function(data, param='auto', lambda.control=2,
                              side.which='upper', sided=1, alpha=0.01, P=0.99,
                              plots=FALSE, plots.nr=FALSE, debug=FALSE, main.adder=NULL) {
    
    ## johnsonsu distribution
    ## MLE (Maximum Likelihood Estimate) fit to determine parameters
    ## LR (Likelihood Ratio) appraoch to find tolerance limit

    ## input: data  = vector of data
    ##        param = initial guess for fit parameters for: gamma, delta, xi, and lambda
    ##                if type is also provided, it will not be used
    ##              = 'auto' (default) uses ExtDist::eJohnsonSU() for initial guess of parameters,
    ##                will switch to SuppDists::JohnsonFit(x) if that fails, and will switch to
    ##                try using list(gamma=xx, delta=xx, xi=xx, lambda=xx) if that fails
    ##        side.which = 'upper', 'lower', 'both' where 'both' calculates upper and lower limits
    ##                     not used if sided = 2
    ##        sided = 1 (default) means 1-sided tolerance limit
    ##              = 2 means 2-sided tolerance limit
    ##        alpha = 1 - confidence
    ##                For 1-sided, 99/99 tolerance limit, specify alpha = 0.01
    ##                For 2-sided, 99/99 tolerance limit, also specify alpha = 0.01
    ##        P     = proportion (or coverage)
    ##                e.g., if 1-sided with P = 0.99, then will base lower and upper tolerance
    ##                       limits on the quantiles associated with 1% and 99% coverage
    ##                e.g., if 2-sided with P = 0.99, then will base lower and upper tolerance limits
    ##                       limits on the quantiles associated with 0.5% and 99.5% coverage
    
    x <- data   
    ## if (isTRUE(plots) & sided == 2) par(mfrow=c(1,2))
    out <- NULL

    ##-----------------------------------------------------------------------------
    ## Johnson SU parameters
    out <- mle.johnsonsu(x)
    params <- out$jparms
    gamma  <- params$gamma
    delta  <- params$delta
    xi     <- params$xi
    lambda <- params$lambda
    params.save <- params
    params$type <- NULL
    
    ##-----------------------------------------------------------------------------
    ## redefine nnl function to fit on desired quantile to find standard error
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
    nll.q <- function(data, param, P, lambda.control, debug=FALSE){
        ## calculate nll (negative log likelihhod) for distribution
        x       <- data
        quant  <- param[[1]]  # replaced gamma with quant as a parameter
        delta  <- param[[2]]
        xi     <- param[[3]]
        lambda <- param[[4]]
        lambda <- lambda.fix(lambda, lambda.control)
        ## write gamma as a function of quant, delta, xi and lambda
        gamma <- qnorm(P) - delta * asinh( (quant-xi)/lambda )
        ## PDF for Johnson SU
        pdf <- delta /( lambda * sqrt(2 * pi)   ) *
            1 / sqrt(1 +                    ( (x-xi)/lambda   )^2)  *
            exp( -0.5*(gamma + delta * asinh( (x-xi)/lambda ) )^2 )
        nll     <- -sum(log(pdf))
        if (isTRUE(debug)) cat('quant=', signif(quant,11),
                               'delta=', signif(delta,11),
                               'xi   =', signif(xi,11),
                               'lambda=', signif(lambda,11), "\n")
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
        if (side.which == 'lower') {
            P.values <- 1-P
        } else if (side.which == 'upper') {
            P.values <- P
        } else {
            ## side.which == 'both'
            P.values <- c(1-P, P)
        }
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
        ## quant.P  <- xi + lambda * sinh( (qnorm(P)-gamma)/delta )
        quant.P <- ExtDist::qJohnsonSU(P, params = params)
        quant.P.save[k] <- quant.P
        quant.param <- c(quant=quant.P, params[2:4])

        ##----------------------
        ## calculate confidence limits using LR (Likelihood Ratio)
        ## confidene limit is defined at likelihood that is lower than max by chi-squared
        ll.max <- -nll.q(x, quant.param, P, lambda.control=lambda.control)
        ## Factor of 2 in the following is because alpha is generally specified for 1-sided
        ## but use for confidence part of tolerance limit is always 2-sided. The difference
        ## in a 1 or 2-sided tolerance limit comes from P, not alpha.
        ll.tol <-  ll.max - qchisq(1 - 2*alpha, 1)/2   # qchisq(1-2*0.01/1, 1) = 5.411894
        cat('MLE=', ll.max, '; tolerance limit at MLE=', ll.tol, '\n')

        ##----------------------
        ## refit on alternate parameters to determine standard error for fit on quantile
        cat('Attempting MLE fit on alternate parameters for P=', P, '\n')
        out.bestfit.q <- optim(par     = quant.param, 
                               fn      = nll.q, 
                               data    = x,
                               P       = P,
                               lambda.control = lambda.control,
                               debug   = debug,
                               control = list(trace=TRUE),
                               hessian = TRUE,
                               method  = "BFGS")
        nll.max.bestfit.q <- out.bestfit.q$value
        quant.P  <- out.bestfit.q$par[[1]]
        delta.P  <- out.bestfit.q$par[[2]]
        xi.P     <- out.bestfit.q$par[[3]]
        lambda.P <- out.bestfit.q$par[[4]]
        lambda.P <- lambda.fix(lambda.P, lambda.control)
        params.q <- as.list(out.bestfit.q$par)
        params.q$lambda <- lambda.P
        params.q$gamma <- qnorm(P) - delta.P * asinh( (quant.P-xi.P)/lambda.P )
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
        dof    <- length(x) - 4    # 4 independent fitting parameters in Johnson SU
        student.t <- qt(1 - alpha/sided, dof)  # 2.3326 for dof=598
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
            plot(quant.P, ll.max, col='blue', pch=16, cex=2,
                 xlab='quantile', ylab='log likelihood',
                 xlim=range(xmin, xmax),
                 ylim=range(ll.max, ll.tol))
            abline(h=ll.tol)
            abline(v=c(quant.P.alpha.l.guess, quant.P.alpha.u.guess), col='red', lty=2)
        }
        
        ##----------------------
        ## function to fit on delta, xi and lambda
        nll.fixedq <- function(data, param, quant, P,
                               lambda.control=lambda.control, debug=FALSE) {
            ## calculate nll (negative log likelihhod) for distribution
            ## for specified quant (i.e., only fit delta, xi, and lambda)
            x       <- data
            delta  <- param[[1]]
            xi     <- param[[2]]
            lambda <- param[[3]]
            lambda <- lambda.fix(lambda, lambda.control)
            ## write gamma as a function of quant, delta, xi and lambda
            gamma <- qnorm(P) - delta * asinh( (quant-xi)/lambda )
            ## PDF for Johnson SU
            pdf <- delta /( lambda * sqrt(2 * pi)   ) *
                1 / sqrt(1 +                    ( (x-xi)/lambda   )^2)  *
                exp( -0.5*(gamma + delta * asinh( (x-xi)/lambda ) )^2 )
            nll     <- -sum(log(pdf))
            if (isTRUE(debug)) cat('quant=', signif(quant,11),
                                   'delta=', signif(delta,11),
                                   'xi   =', signif(xi,11),
                                   'lambda=', signif(lambda,11),
                                   'nll   =', signif(nll,11), "\n")
            return(nll)
        }

        ll.fixedq <- function(x0, data, P, delta=delta.P, xi=xi.P, lambda=lambda.P,
                              lambda.control=lambda.control, debug=FALSE) {
            ## first determine best fit delta, xi, and lambda for given x0=quant (and P)
            fit <- NA
            tryCatch({
                fit <- optim(par     = c(delta, xi, lambda), 
                             fn      = nll.fixedq, 
                             data    = data,
                             quant   = x0,
                             P       = P,
                             lambda.control = lambda.control,
                             debug   = debug,
                             control = list(trace=TRUE,
                                            maxit=1e4),
                             hessian = TRUE,
                             method  = "BFGS")
                ll <- -fit$value
            }, error = function(e) {
                ## what to do if error
                cat('WARNING: CONVERGENCE FAILURE IN mle.johnsonsu.tol.r MODULE\n')
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
        mini.test <- function() {
            x <- iris$Sepal.Width
            P <- 0.90             # proportion or coverage
            out <- mle.johnsonsu(x)$jparms
            delta_test <- out$delta
            ## this one converges nicely
            ll.fixedq(4.4, x, P, delta=out$delta, xi=out$xi, lambda=out$lambda, debug=FALSE)
            ## this one does not converge
            ll.fixedq(4.5, x, P, delta=out$delta, xi=out$xi, lambda=out$lambda, debug=FALSE)
            ## this one converges nicely
            ll.fixedq(4.6, x, P, delta=out$delta, xi=out$xi, lambda=out$lambda, debug=FALSE)
        }

        ##----------------------
        if (isTRUE(plots)) {
            xplot <- seq(xmin, xmax, length.out=301)
            yplot <- NA
            for (ploti in 1:length(xplot)) {
                ## if (ploti == 3) browser()
                yplot[ploti] <- ll.fixedq(xplot[ploti], x, P, delta, xi, lambda,
                                          lambda.control=lambda.control, debug=FALSE)
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
                yplot <- ll.fixedq(x0     = xplot,
                                   data   = x,
                                   P      = P,
                                   delta  = delta.P,
                                   xi     = xi.P,
                                   lambda = lambda.P,
                                   lambda.control = lambda.control,
                                   debug  = FALSE)
                points(xplot, yplot, col=col, pch=16, cex=2)
                return(yplot)
            }
        }

        ## determine confidence bound (P alread determined whether this was a lower or upper bound)
        ## as point where the likelihood ratio equals 11.tol
        quant.P.alpha.l <- NA
        quant.P.alpha.u <- NA
        if (P <= 0.5 | side.which == 'lower') {
            out.nrl <- newton.raphson(f = ll.fixedq,
                                      xguess = quant.P.alpha.l.guess,
                                      ytarget = ll.tol,
                                      data   = x,
                                      P      = P,
                                      delta  = delta,
                                      xi     = xi,
                                      lambda = lambda,
                                      lambda.control = lambda.control,
                                      tol = 1e-5, n = 1000,
                                      plots=plots.nr)
            quant.P.alpha.l <- out.nrl$root
            if (is.null(quant.P.alpha.l)) { quant.P.alpha.l <- NA }
        } else if (P >= 0.5 | side.which == 'upper') {
            ## if (P > 0.5) browser()
            out.nru <- newton.raphson(f = ll.fixedq,
                                      xguess = quant.P.alpha.u.guess,
                                      ytarget = ll.tol,
                                      data   = x,
                                      P      = P,
                                      delta  = delta,
                                      xi     = xi,
                                      lambda = lambda,
                                      lambda.control = lambda.control,
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


        new.idea.not.working <- function() {
            ##----------------------
            ## alternate use of optim to find tolerance limit
            quant <- mean(x) + (max(x)-mean(x))/2
            tol.zero <- function(data, param, P, debug, ll.tol) {
                nll  <- nll.q(data, param, P, debug=debug)
                zero <- -abs(nll + ll.tol)
                return(zero)
            }
            out <- ll.fixedq(quant.P.alpha.u, x, P,
                             delta=delta, xi=xi, lambda=lambda, debug=TRUE)
            out.nll <- out$fit$value
            out.par <- out$fit$par
            tol.zero(x, c(quant.P.alpha.u-0.1, out.par), P, debug=FALSE, ll.tol)
            tol.zero(x, c(quant.P.alpha.u    , out.par), P, debug=FALSE, ll.tol)
            tol.zero(x, c(quant.P.alpha.u+0.1, out.par), P, debug=FALSE, ll.tol)

            fit <- optim(par     = c(quant, delta, xi, lambda), 
                         fn      = tol.zero,
                         data    = x,
                         P       = P,
                         debug   = FALSE,
                         ll.tol  = ll.tol,
                         control = list(trace=TRUE,
                                        maxit=1e4),   # a bit better without maxit
                         hessian = TRUE,
                         method  = "BFGS")

            
            mle.tol <- function(data, P, quant=quant, delta=delta, xi=xi, lambda=lambda, ll.tol) {
                ## first determine best fit delta, xi, and lambda for given x0=quant (and P)
                fit <- NA
                tryCatch({
                    fit <- optim(par     = c(quant, delta, xi, lambda), 
                                 fn      = tol.zero,
                                 data    = data,
                                 P       = P,
                                 debug   = FALSE,
                                 ll.tol  = ll.tol,
                                 control = list(trace=TRUE,
                                                maxit=1e4),
                                 hessian = TRUE,
                                 method  = "BFGS")
                    ll <- -fit$value
                    quant <- fit$par$quant
                    delta <- fit$par$delta
                    xi    <- fit$par$xi
                    lambda <- fit$par$lambda
                }, error = function(e) {
                    ## what to do if error
                    cat('WARNING: CONVERGENCE FAILURE IN mle.johnsonsu.tol.r MODULE\n')
                    cat('         IN mle.tol OPTIMIZATION ROUTINE.\n\n')
                    fit <- NA
                    ll  <- NA
                })
                return(list(fit=fit,
                            ll=ll,
                            parma=list(quant=quant, delta=delta, xi=xi, lambda=lambda)))
            }
            out.mle <- mle.tol(x, P, quant=quant, delta=delta, xi=xi, lambda=lambda, ll.tol)
        }

    }

    ## collect lower and upper tolerance limits
    tol.approx <- range(tol.approx, na.rm=TRUE)
    tol.limits <- range(tol.limits, na.rm=TRUE)
    if (side.which == 'both') {
        tol.lower <- tol.limits[1]
        tol.upper <- tol.limits[2]
    } else if (side.which == 'lower') {
        tol.lower <- tol.limits[1]
        tol.upper <- NA
    } else {
        # side.which == 'upper'
        tol.lower <- NA
        tol.upper <- tol.limits[2]
    }
    ## reset P to input value
    P <- P.in

    ## collect tolerance values in dataframe similar to extol.int for weibull
    tolerance <- data.frame(alpha, P, sided, tol.lower, tol.upper)

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


mle.johnsonsu.tol.test <- function() {
    source('setup.r')

    ## create Johnson SU dataset
    set.seed(1)
    jparms <- list(gamma=-3.3, delta=5.3, xi=1.8, lambda=1.9, type='SU')
    x <- ExtDist::rJohnsonSU(1000, param=jparms)

    ## determine johnson su parameters
    plotspace(2,2)
    out.fit <- mle.johnsonsu(x, plots=TRUE)
    jparms.mle <- out.fit$jparms
    print(out.fit$jparms.compare)

    ## johnson su tolerance limits
    out.tol <- mle.johnsonsu.tol(x, plots=TRUE)

    ## incorporated into other functions
    plotspace(2,2)
    out.hist <- hist_nwj(x)
    out <- qqplot_nwj(x, type='n')
    out <- qqplot_nwj(x, type='w')
    out <- qqplot_nwj(x, type='j')

    ## consider the following dataset
    x <- iris$Sepal.Width
    plotspace(3,2)
    ## lower tolerance limit
    out.lower <- mle.johnsonsu.tol(data=x, param='auto', lambda.control=2,
                                   side.which='lower', sided=1, alpha=0.01, P=0.99,
                                   plots=TRUE, plots.nr=FALSE, debug=FALSE, main.adder='lower, 1-sided')
    ## using default parameters except for 'plots'
    out.upper <- mle.johnsonsu.tol(data=x, param='auto', lambda.control=2,
                                   side.which='upper', sided=1, alpha=0.01, P=0.99,
                                   plots=TRUE, plots.nr=FALSE, debug=FALSE, main.adder='upper, 1-sided')
    ## lower and upper 1-sided tolerance limits
    out.both <- mle.johnsonsu.tol(data=x, param='auto', lambda.control=2,
                                   side.which='both', sided=1, alpha=0.01, P=0.99,
                                   plots=TRUE, plots.nr=FALSE, debug=FALSE, main.adder='both, 1-sided')
    ## lower and upper 2-sided tolerance limits
    out.twosided <- mle.johnsonsu.tol(data=x, param='auto', lambda.control=2,
                                      side.which='both', sided=2, alpha=0.01, P=0.99,
                                      plots=TRUE, plots.nr=FALSE, debug=FALSE, main.adder='2-sided')
}
