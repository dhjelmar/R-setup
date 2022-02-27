mle.johnsonsu <- function(data, param='auto', fit.only=FALSE, alpha=0.01, P=0.99, sided=1, plots=FALSE, debug=FALSE) {
    
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
    ## https://personal.psu.edu/abs12/stat504/Lecture/lec3_4up.pdf
    
    x <- data   
    if (isTRUE(plots)) par(mfrow=c(1,2))
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
    if (param[1] == 'auto') {
        ## force the Johnson SU distribution
        jparms.out <- ExtDist::eJohnsonSU(x)
        jparms <- list(gamma   = jparms.out$gamma,
                       delta   = jparms.out$delta,
                       xi      = jparms.out$xi,
                       lambda  = jparms.out$lambda,
                       type    = 'SU')
        if (jparms$gamma == -0.5 & jparms$delta == 2.0 & jparms$xi == -0.5 & jparms$lambda == 2.0) {
            cat('\n')
            cat('ExtDist::eJohnsonSU(x) fit failed\n')
            print(as.data.frame(jparms[1:4]))
            cat('\n')
            ## try another packages
            jparms <- NA
            jparms <- SuppDists::JohnsonFit(x)
            if (jparms$type != 'SU') {
                cat('SuppDists::JohnsonFit(x) failed to converge or did not return type=SU\n')
                print(as.data.frame(jparms))
                cat('\nSetting jparms to:\n')
                jparms <- list(gamma   = 1,
                               delta   = 1,
                               xi      = 1,
                               lambda  = 1,
                               type    = 'SU')
            } else {
                cat('SuppDists::JohnsonFit(x) used instead for initial parameter guesses\n')
                print(as.data.frame(jparms))
                cat('\n')
            }
        }

    } else {
        ## use Johnson parameters specified in param
        ## needs to be in same list format as created by SuppDists::JohnsonFit
        jparms <- param
    }
    ## strip off type if that is provided in param list
    param <- jparms[1:4]   # list(gamma, delta, xi, lambda)
    params.guess <- param
    
    ##-----------------------------------------------------------------------------
    ## determine best fit using nnl
    nll <- function(data, param, debug=FALSE){
        ## calculate nll (negative log likelihhod) for distribution
        x      <- data
        gamma  <- param[[1]]  
        delta  <- param[[2]]
        xi     <- param[[3]]
        lambda <- param[[4]]
        ## PDF for Johnson SU
        pdf <- delta /( lambda * sqrt(2 * pi)   ) *
            1 / sqrt(1 +            ( (x-xi)/lambda )^2)  *
            exp( -0.5*(gamma + delta * asinh( (x-xi)/lambda ))^2 )
        ## the above is equivalent
        ## pdf <- ExtDist::dJohnsonSU(x, parms=c(gamma, delta, xi, lambda))
        nll     <- -sum(log(pdf))
        if (isTRUE(debug)) cat('gamma=', signif(gama,11),
                               'delta=', signif(delta,11),
                               'xi   =', signif(xi, 11),
                               'lambda=', signif(lambda,11),
                               'nll=', signif(nll,11), "\n")
        return(nll)
    }        
    print('Attempting MLE fit on regular parameters')
    out.bestfit <- stats::optim(par     = param, 
                                fn      = nll, 
                                data    = x,
                                debug   = debug,
                                control = list(trace=TRUE),
                                method  = "BFGS")
    nll.max.bestfit <- out.bestfit$value
    gamma  <- out.bestfit$par[[1]]
    delta  <- out.bestfit$par[[2]]
    xi     <- out.bestfit$par[[3]]
    lambda <- out.bestfit$par[[4]]
    params <- as.list(out.bestfit$par)
    print(as.data.frame(params))
    cat('\n')

    if (isTRUE(fit.only)) {
        jparms <- params
        jparms$type <- 'SU'
        return(list(params=jparms))
    }
    
    ##-----------------------------------------------------------------------------
    ## redefine nnl function to fit on desired quantile
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
    P.lower.upper <- c(1-P, P)
    for (P in P.lower.upper) {

        k <- k+1  # counter

        ##----------------------
        ## fist set estimated parameters
        quant.P <- ExtDist::qJohnsonSU(P, params = params)
        quant.param <- c(quant=quant.P, params[2:4])

        ##----------------------
        ## refit to find quantile, quant.P, associated with the fit
        print('Attempting MLE fit on regular parameters for P=', P)
        out.bestfit.q <- stats::optim(par     = quant.param, 
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
        cat('\n')
        
        ##----------------------
        ## estimate confidence limit using standard error to serve as starting point for search
        ## estimates      <- as.numeric( out.bestfit.q$par )
        standard.error <- as.numeric( sqrt(diag(solve(out.bestfit.q$hessian))) )
        ## coef           <- data.frame(estimates, standard.error)
        ## rownames(coef) <- names(out.bestfit.q$par)
        dof    <- length(x) - 4    # 4 independent fitting parameters in Johnson SU
        student.t <- qt(1 - alpha/sided, dof)  # 2.3326 for dof=598
        quant.P.alpha.l.guess <- quant.P - student.t * standard.error[1]
        quant.P.alpha.u.guess <- quant.P + student.t * standard.error[1]
        cat('Initial guesses for confidence interval for P=', P, '\n')
        cat(quant.P.alpha.l.guess, quant.P.alpha.u.guess, '\n\n')
        
        ##----------------------
        ## calculate confidence limits using LR (Likelihood Ratio)
        ## confidene limit is defined at likelihood that is lower than max by chi-squared
        ll.max.P <- -nll.max.bestfit
        ll.tol <-  ll.max.P - qchisq(1 - alpha/sided, 1)   # qchisq(1-0.01/1, 1) = 6.634897 
        cat('MLE=', ll.max.P, '; tolerance limit at MLE=', ll.tol, 'for P=', P, '\n')
        
        ## function for newton.raphson() iterates on x0
        ll.q <- function(x0, data, P, delta, xi, lambda) {
            ll <- -nll.q(data,
                         param=list(quant  = x0,
                                    delta  = delta,
                                    xi     = xi,
                                    lambda = lambda),
                         P,
                         debug=FALSE)
        }
        
        ## determine confidence bound (P alread determined whether this was a lower or upper bound)
        ## as point where the likelihood ratio equals 11.tol
        out.nrl <- newton.raphson(f = ll.q,
                                  xguess = quant.P.alpha.l.guess,
                                  ytarget = ll.tol,
                                  data   = x,
                                  P      = P,
                                  delta  = delta,
                                  xi     = xi,
                                  lambda = lambda,
                                  tol = 1e-5, n = 1000, plot='no')
        quant.P.alpha.l <- out.nrl$root
        out.nru <- newton.raphson(f = ll.q,
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
            ## plot the likelihood as a function of the quantile

            ## set range of quantile values to plot
            quant     <- seq(quant.P.alpha.l, quant.P.alpha.u, length.out=100)

            ## calculate and plot log likelihood
            ll <- NA
            for (i in 1:100) {
                ## vary quant parameter to see impact on likelihood
                param.vary.q <- list(quant  = quant[i],
                                     delta  = delta,
                                     xi     = xi,
                                     lambda = lambda)
                ll[i] <- -nll.q(x, param.vary.q, P=P, debug=FALSE)
            }
            plot(quant, ll, xlab='quantile', ylab='Log Likelihood')

            ## plot log likelihood at tolerance limit
            abline(h=ll.tol)

            ## plot intersection with log likelihood curve
            points(quant.P.alpha.l, ll.tol, col='red', pch=16, cex=2)
            points(quant.P.alpha.u, ll.tol, col='red', pch=16, cex=2)

        }


    }

    ## put parameters from the standard fit and fit on quantile P into dataframe for comparison
    a0 <- as.data.frame(t(unlist(params.guess)))
    a1 <- as.data.frame(t(unlist(params)))
    a2 <- as.data.frame(t(unlist(params.q.save[1])))
    a3 <- as.data.frame(t(unlist(params.q.save[2])))
    params.compare <- fastmerge(a0, a1)
    params.compare <- fastmerge(params.compare, a2)
    params.compare <- fastmerge(params.compare, a3)
    rownames(params.compare) <- c('initial guess parameters', 'standard fit',
                                  'fit on quantile at 1-P', 'fit on quantile at P')

    ## collect lower and upper tolerance limits
    tol.limits <- range(tol.limits, na.rm=TRUE)
    tol.lower <- tol.limits[1]
    tol.upper <- tol.limits[2]
    tolerance <- data.frame(alpha, P, sided, tol.lower, tol.upper)

    ## add type to params list for use in other modules
    params$type <- 'SU'

    ## print final parameter comparison and tolerance limits
    print(params.compare)
    cat('\n')
    print(tolerance)
    cat('\n')
    
    return(list(params     = params,
                params.q   = params.q,
                params.compare = params.compare,
                out.bestit = out.bestfit,
                tolerance  = tolerance))
}


mle.johnsonsu.test <- function() {
    source('modules/newton.raphson.r')
    source('modules/fastmerge.r')

    ## consider the following dataset
    x <- iris$Sepal.Width
    mean(x)  # 3.057333
    sd(x)    # 0.4358663
    P     <- 0.99  # proportion or coverage
    alpha <- 0.01
    sided <- 1

    out <- mle.johnsonsu(x, 'auto', alpha=alpha, P=P, sided=sided, plots=TRUE, debug=FALSE)
    print(out$params.compare)
    print(out$tolerance)

    jparms <- ExtDist::eJohnsonSU(x)
    jparms <- list(gamma=jparms[[1]], delta=jparms[[2]], xi=jparms[[3]], lambda=jparms[[4]])
    out <- mle.johnsonsu(x, jparms, alpha=alpha, P=P, sided=sided, plots=TRUE, debug=FALSE)
    print(out$params.compare)
    print(out$tolerance)

    jparms2 <- list(gamma=jparms[[1]], delta=jparms[[2]], xi=jparms[[3]], lambda=jparms[[4]], type='SU')
    out2 <- johnson_tol(x, jparms2)
    
    ## test inside other modules
    plotspace(2,2)
    out.h <- hist_nwj(x, type = 'nwj', mle=TRUE, jfit='auto')
    out.n <- qqplot_nwj(x, type='n', mle=TRUE)
    out.w <- qqplot_nwj(x, type='w', mle=TRUE)
    out.j <- qqplot_nwj(x, type='j', mle=TRUE)
    
    out <- mle.johnsonsu(x, 'auto', plots=TRUE)
}
