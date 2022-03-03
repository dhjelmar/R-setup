mle.johnsonsu.tol <- function(data, param='auto', alpha=0.01, P=0.99, sided=1, plots=FALSE, debug=FALSE) {
    
    ## johnsonsu distribution
    ## MLE (Maximum Likelihood Estimate) fit to determine parameters
    ## LR (Likelihood Ratio) appraoch to find tolerance limit

    ## input: data  = vector of data
    ##        param = initial guess for fit parameters for: gamma, delta, xi, and lambda
    ##                if type is also provided, it will not be used
    ##              = 'auto' (default) uses ExtDist::eJohnsonSU() for initial guess of parameters,
    ##                will switch to SuppDists::JohnsonFit(x) if that fails, and will switch to
    ##                try using list(gamma=xx, delta=xx, xi=xx, lambda=xx) if that fails
    
    x <- data   
    if (isTRUE(plots)) par(mfrow=c(1,2))   ## comment out if want to restrict to P and not also 1-P
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
    nll.q <- function(data, param, P, debug=FALSE){
        ## calculate nll (negative log likelihhod) for distribution
        x       <- data
        quant  <- param[[1]]  # replaced gamma with quant as a parameter
        delta  <- param[[2]]
        xi     <- param[[3]]
        lambda <- param[[4]]
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
    params.q.save <- NA
    k <- 0
    P.lower.upper <- c(1-P, P)   ## comment out if want to restrict to P and not also 1-P
    for (P in P.lower.upper) {   ## comment out if want to restrict to P and not also 1-P

        k <- k+1  # counter

        ##----------------------
        ## fist set estimated parameters
        ## following is identical to
        ## quant.P  <- xi + lambda * sinh( (qnorm(P)-gamma)/delta )
        quant.P <- ExtDist::qJohnsonSU(P, params = params)
        quant.param <- c(quant=quant.P, params[2:4])

        ##----------------------
        ## calculate confidence limits using LR (Likelihood Ratio)
        ## confidene limit is defined at likelihood that is lower than max by chi-squared
        ll.max <- -nll.q(x, quant.param, P)
        ll.tol <-  ll.max - qchisq(1 - alpha/sided, 1)/2   # qchisq(1-0.01/1, 1) = 6.634897 
        cat('MLE=', ll.max, '; tolerance limit at MLE=', ll.tol, '\n')

        ##----------------------
        ## refit on alternate parameters to determine standard error for fit on quantile
        cat('Attempting MLE fit on alternate parameters for P=', P, '\n')
        out.bestfit.q <- optim(par     = quant.param, 
                               fn      = nll.q, 
                               data    = x,
                               P       = P,
                               debug   = debug,
                               control = list(trace=TRUE),
                               hessian = TRUE,
                               method  = "BFGS")
        nll.max.bestfit.q <- out.bestfit.q$value
        quant.P  <- out.bestfit.q$par[[1]]
        delta.P  <- out.bestfit.q$par[[2]]
        xi.P     <- out.bestfit.q$par[[3]]
        lambda.P <- out.bestfit.q$par[[4]]
        params.q <- as.list(out.bestfit.q$par)
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


        ##----------------------
        if (isTRUE(plots)) {
            quant.dif <- quant.P.alpha.u.guess - quant.P
            xmin <- quant.P.alpha.l.guess - 2*quant.dif
            xmax <- quant.P.alpha.u.guess + 2*quant.dif
            xplot <- seq(xmin, xmax, length.out=101)
            yplot <- NA
            for (ploti in 1:101) {
                qparms <- list(quant  = xplot[ploti],
                               delta  = quant.param$delta,
                               xi     = quant.param$xi,
                               lambda = quant.param$lambda)
                yplot[ploti] <- -nll.q(x, qparms, P)
            }
            xyplot <- data.frame(quantile       = xplot,
                                 log.likelihood = yplot)
            plot(xyplot$quantile, xyplot$log.likelihood, col='black',
                 xlab='quantile', ylab='log likelihood',
                 ylim=range(ll.max, ll.tol))
            points(quant.P, ll.max, col='blue', pch=16, cex=2)
            abline(h=ll.tol)
        }
        
        ##----------------------
        ## function to fit on delta, xi and lambda
        nll.fixedq <- function(data, param, quant, P, debug=FALSE) {
            ## calculate nll (negative log likelihhod) for distribution
            ## for specified quant (i.e., only fit delta, xi, and lambda)
            x       <- data
            delta  <- param[[1]]
            xi     <- param[[2]]
            lambda <- param[[3]]
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

        ll.fixedq <- function(x0, data, P, delta=delta.P, xi=xi.P, lambda=lambda.P, debug=FALSE) {
            ## first determine best fit delta, xi, and lambda for given x0=quant (and P)
            fit <- NA
            tryCatch({
                fit <- optim(par     = c(delta, xi, lambda), 
                             fn      = nll.fixedq, 
                             data    = data,
                             quant   = x0,
                             P       = P,
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
            x <-     x <- iris$Sepal.Width
            ## P     <- 0.99  # proportion or coverage
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
            xplot <- seq(xmin, xmax*1.2, length.out=301)
            yplot <- NA
            for (ploti in 1:301) {
                yplot[ploti] <- ll.fixedq(xplot[ploti], x, P, delta, xi, lambda, debug=FALSE)
                cat(xplot[ploti], yplot[ploti], '\n')
            }
            xyplot <- data.frame(quantile       = xplot,
                                 log.likelihood = yplot,
                                 ll.tol         = ll.tol)
            points(xyplot$quantile, xyplot$log.likelihood, col='blue')
            ## converged <- which(!is.na(yplot))
            ## points(xplot[converged], yplot[converged],
            ##        xlab='quantile', ylab='log likelihood', col='red')
        }

        ## determine confidence bound (P alread determined whether this was a lower or upper bound)
        ## as point where the likelihood ratio equals 11.tol
        out.nrl <- newton.raphson(f = ll.fixedq,
                                  xguess = quant.P.alpha.l.guess,
                                  ytarget = ll.tol,
                                  data   = x,
                                  P      = P,
                                  delta  = delta,
                                  xi     = xi,
                                  lambda = lambda,
                                  tol = 1e-5, n = 1000, plot='no')
        quant.P.alpha.l <- out.nrl$root
        ## if (P > 0.5) browser()
        out.nru <- newton.raphson(f = ll.fixedq,
                                  xguess = quant.P.alpha.u.guess,
                                  ytarget = ll.tol,
                                  data   = x,
                                  P      = P,
                                  delta  = delta,
                                  xi     = xi,
                                  lambda = lambda,
                                  tol = 1e-5, n = 1000, plot='no')
        quant.P.alpha.u <- out.nru$root
        cat('Final confidence interval for P=', P, '\n')
        cat(quant.P.alpha.l, quant.P.alpha.u, '\n\n')

        ## collect tolerance limit calculations (i.e., for 1-P and P)
        tol.limits <- c(tol.limits, quant.P.alpha.l, quant.P.alpha.u)
        
        if (isTRUE(plots)) {
            ## plot intersection with log likelihood curve
            points(quant.P.alpha.l, ll.tol, col='red', pch=16, cex=2)
            points(quant.P.alpha.u, ll.tol, col='red', pch=16, cex=2)
        }


    }   ## comment out if want to restrict to P and not also 1-P

    ## collect lower and upper tolerance limits
    tol.limits <- range(tol.limits, na.rm=TRUE)
    tol.lower <- tol.limits[1]
    tol.upper <- tol.limits[2]
    tolerance <- data.frame(alpha, P, sided, tol.lower, tol.upper)

    ## add type to params list for use in other modules
    params$type <- 'SU'

    ## print final parameter comparison and tolerance limits
    print(params.save)
    cat('\n')
    print(tolerance)
    cat('\n')
    
    return(list(params     = params.save,
                tolerance  = tolerance))
}


mle.johnsonsu.test <- function() {
    source('setup.r')

    ## consider the following dataset
    x <- iris$Sepal.Width
    P     <- 0.99  # proportion or coverage
    P     <- 0.90  # proportion or coverage
    alpha <- 0.01
    sided <- 1

    out <- mle.johnsonsu.tol(x, 'auto', alpha=alpha, P=P, sided=sided, plots=TRUE, debug=FALSE)
    print(out$params)
    print(out$tolerance)
    
    ## test inside other modules
    plotspace(2,2)
    out.h <- hist_nwj(x, type = 'nwj', mle=TRUE, jfit='auto')
    out.n <- qqplot_nwj(x, type='n', mle=TRUE)
    out.w <- qqplot_nwj(x, type='w', mle=TRUE)
    out.j <- qqplot_nwj(x, type='j', mle=TRUE)
    
    out <- mle.johnsonsu(x, 'auto', plots=TRUE)
}
