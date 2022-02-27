mle.normal <- function(data, param='auto', fit.only=FALSE, alpha=0.01, P=0.99, sided=1, plots=FALSE, debug=FALSE) {
    
    ## normal distribution
    ## MLE (Maximum Likelihood Estimate) fit to determine parameters
    ## LR (Likelihood Ratio) appraoch to find tolerance limit

    ## input: data  = vector of data
    ##        param = initial guess for fit parameters for: xbar and sdev
    ##                if type is also provided, it will not be used
    ##              = 'auto' uses mean(x) and sd(x) to determine initial guess

    ## based on approach found here:
    ## https://www.r-bloggers.com/2019/08/maximum-likelihood-estimation-from-scratch/
    ## https://personal.psu.edu/abs12/stat504/Lecture/lec3_4up.pdf
    
    x <- data    
    if (isTRUE(plots)) par(mfrow=c(1,2))
    out <- NULL

    ##-----------------------------------------------------------------------------
    ## initial guesses for Johnson parameters
    if (param[1] == 'auto') {
        ## use R to guess fit parameters
        xbar <- mean(x)
        sdev <- sd(x)
        param <- list(xbar=xbar, sdev=sdev)
    }
    params.guess <- param
    
    ##-----------------------------------------------------------------------------
    ## determine best fit using nnl
    nll <- function(data, param, debug=FALSE){
        ## calculate nll (negative log likelihhod) for distribution
        x       <- data
        xbar    <- param[[1]]
        sdev    <- param[[2]]
        z       <- (x - xbar) / sdev
        pdf     <- 1/(sqrt(2 * pi) * sdev) * exp(-z^2 / 2) 
        ## the following is equivalent
        ## pdf     <- dnorm(x, mean = xbar, sd = sdev, log = FALSE)
        nll     <- -sum(log(pdf))
        if (isTRUE(debug)) cat('xbar=', signif(xbar,11), 'sdev=', signif(sdev,11), 'nll=', signif(nll,11), "\n")
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
    xbar <- out.bestfit$par[[1]]
    sdev <- out.bestfit$par[[2]]
    params <- list(xbar=xbar, sdev=sdev)
    print(as.data.frame(params))
    cat('\n')

    if (isTRUE(fit.only)) return(list(params=params))
    
    ##-----------------------------------------------------------------------------
    ## redefine nnl function to fit on desired quantile
    nll.q <- function(data, param, P, debug=FALSE){
        ## calculate nll (negative log likelihhod) for distribution
        x       <- data
        quant   <- param[[1]]  # substituted for xbar
        sdev    <- param[[2]]
        ## the following does not work well
        ## xbar    <- quant - sqrt( 2*sdev^2 * (-log(qnorm(P)) - log(sqrt(2*pi)*sdev)) )      
        ## repalced it using an alternate form of the normal PDF found here:
        ## F(x) = 0.5 *( 1 + erf( (x-xbar)/(sdev * sqrt(2))))
        xbar <- quant - sqrt(2) * sdev * pracma::erfinv(2*P-1)
        z       <- (x - xbar) / sdev
        pdf     <- 1/(sqrt(2 * pi) * sdev) * exp(-z^2 / 2) 
        nll     <- -sum(log(pdf))
        if (isTRUE(debug)) cat('quant=', signif(quant,11),
                               'sdev=', signif(sdev,11),
                               'nll=', signif(nll,11), "\n")
        return(nll)
    }

    ##---------------------- aside
    ## f <- function(x, xbar, sdev) {
    ##   1/(sqrt(2*pi)*sdev) * exp(-(x-xbar)^2/(2*sdev^2))
    ## }
    ## f(xbar, xbar, sdev)
    ## par(mfrow=c(2,2))
    ## curve(f(x, xbar, sdev), xbar-3*sdev, xbar+3*sdev, xlab='x=observation', ylab='PDF')
    ## curve(dnorm(x, xbar, sdev), xbar-3*sdev, xbar+3*sdev, xlab='x=observation') # dnorm(x, mean, sd)
    ## curve(pnorm(x, xbar, sdev), 1, 5, xlab='x=quantile')                # pnorm(q, mean, sd)
    ## curve(qnorm(x, xbar, sdev), 0, 1, xlab='x=probability')                     # qnorm(p, mean, sd)
    ##---------------------- aside
    
    
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
        quant.P <- qnorm(P, xbar, sdev)   # for p=0.5, this should be xbar
        quant.param <- c(quant=quant.P, sdev=sdev)

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
        sdev.P   <- out.bestfit.q$par[[2]]
        params.q <- as.list(out.bestfit.q$par)
        params.q$xbar <- quant.P - sqrt(2) * sdev.P * pracma::erfinv(2*P-1)

        params.q.save[k] <- list(params.q)
        print(as.data.frame(params.q))
        cat('\n')
        
        ##----------------------
        ## estimate confidence limit using standard error to serve as starting point for search
        ## estimates      <- as.numeric( out.bestfit.q$par )
        standard.error <- as.numeric( sqrt(diag(solve(out.bestfit.q$hessian))) )
        ## coef           <- data.frame(estimates, standard.error)
        ## rownames(coef) <- names(out.bestfit.q$par)
        dof    <- length(x) - 2    # 2 independent fitting parameters in normal fit
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
        ll.q <- function(x0, data, P, sdev) {
            ll <- -nll.q(data,
                         param=list(quant  = x0,
                                    sdev  = sdev),
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
                                  sdev  = sdev,
                                  tol = 1e-5, n = 1000, plot='no')
        quant.P.alpha.l <- out.nrl$root
        out.nru <- newton.raphson(f = ll.q,
                                  xguess = quant.P.alpha.u.guess,
                                  ytarget = ll.tol,
                                  data   = x,
                                  P      = P,
                                  sdev  = sdev,
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
                                     sdev   = sdev)
                ll[i] <- -nll.q(x, param.vary.q, P=P, debug=FALSE)
            }
            plot(quant, ll, xlab='quantile', ylab='Log Likelihood')

            ## plot mean
            abline(v=xbar)

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


mle.normal.test <- function() {
    source('modules/newton.raphson.r')
    source('modules/fastmerge.r')

    ## consider the following dataset
    x <- iris$Sepal.Width
    mean(x)  # 3.057333
    sd(x)    # 0.4358663
    alpha <- 0.01
    P     <- 0.99  # proportion or coverage
    sided <- 1
    tol_out_norm <- tolerance::normtol.int(x, alpha = alpha, P=P, side=sided)
    upper_tol <- tol_out_norm$'1-sided.upper'
    upper_tol

    ## out <- mle(x, c(mean=0, sd=1), alpha=0.01, P=0.5 , sided=1, plots=TRUE, debug=FALSE)
    P     <- 0.99
    alpha <- 0.01
    sided <- 1
    out <- mle.normal(x, c(mean=0, sd=1), alpha=alpha, P=P, sided=sided, plots=TRUE, debug=FALSE)
    print(out$params.compare)
    print(out$tolerance)
    print(tolerance::normtol.int(x, alpha = alpha, P=P, side=sided))
}
