mle.weibull <- function(data, param, fit='n', alpha=0.01, P=0.99, sided=1, plots=FALSE, debug=FALSE) {
    
    ## weibull distribution
    ## MLE (Maximum Likelihood Estimate) fit to determine parameters
    ## LR (Likelihood Ratio) appraoch to find tolerance limit

    ## input: data  = vector of data
    ##        param = initial guess for fit parameters for: gamma, delta, xi, and lambda
    ##                if type is also provided, it will not be used
    
    x <- data    
    if (isTRUE(plots)) par(mfrow=c(1,2))
    out <- NULL

    ##-----------------------------------------------------------------------------
    ## initial guesses for Johnson parameters
    if (param[1] == 'exttol.int') {
        ## use R to guess fit parameters
        tol_out_weib <-  tolerance::exttol.int(x, alpha=alpha, P=P, side=1, dist="Weibull")
        shape   <- tol_out_weib$'shape.1'
        scale   <- tol_out_weib$'shape.2'
        param <- list(shape=shape, scale=scale)
    }
    
    ##-----------------------------------------------------------------------------
    ## determine best fit using nnl
    nll <- function(data, param, debug=FALSE){
        ## calculate nll (negative log likelihhod) for weibull distribution
        x      <- data
        shape  <- param[[1]]  
        scale  <- param[[2]]
        pdf <- shape / scale^shape * x^(shape-1) * exp(-(x/scale)^shape)
        ## the above is equivalent
        ## pdf <- dWeibull(x, shape, scale)
        nll     <- -sum(log(pdf))
        if (isTRUE(debug)) cat('shape=', signif(xbar,11), 'scale=', signif(sdev,11), 'nll=', signif(nll,11), "\n")
        return(nll)
    }        
    out.bestfit <- optim(par     = param, 
                         fn      = nll, 
                         data    = x,
                         debug   = debug,
                         control = list(trace=TRUE),
                         method  = "BFGS")
    nll.max.bestfit <- out.bestfit$value
    shape  <- out.bestfit$par[[1]]
    scale  <- out.bestfit$par[[2]]
    params <- as.list(out.bestfit$par)

    ##-----------------------------------------------------------------------------
    ## redefine nnl function to fit on desired quantile
    nll.q <- function(data, param, P, debug=FALSE){
        ## calculate nll (negative log likelihhod) for weibull distribution
        x       <- data
        quant  <- param[[1]]  # replaced gamma with quant as a parameter
        shape  <- param[[2]]
        scale  <- quant/(-log((1-P)))^(1/shape)
        pdf <- shape / scale^shape * x^(shape-1) * exp(-(x/scale)^shape)
        ## the above is equivalent
        ## pdf <- dWeibull(x, shape, scale)
        nll     <- -sum(log(pdf))
        if (isTRUE(debug)) cat('quant=', signif(quant,11),
                               'shape=', signif(shape,11),
                               'scale=', signif(scale,11),"\n")
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
        quant.P <- qweibull(P, shape=shape, scale=scale)
        quant.param <- list(quant=quant.P, shape=shape)

        ##----------------------
        ## refit to find quantile, quant.P, associated with the fit
        out.bestfit.q <- optim(par     = quant.param, 
                               fn      = nll.q, 
                               data    = x,
                               P       = P,
                               debug   = debug,
                               control = list(trace=TRUE),
                               method  = "BFGS")
        nll.max.bestfit.q <- out.bestfit.q$value
        quant.P  <- out.bestfit.q$par[[1]]
        shape.P  <- out.bestfit.q$par[[2]]
        params.q <- as.list(out.bestfit.q$par)
        params.q$scale  <- quant.P/(-log((1-P)))^(1/shape.P)

        params.q.save[k] <- list(params.q)
        
        ##----------------------
        ## calculate confidence limits using LR (Likelihood Ratio)
        ## confidene limit is defined at likelihood that is lower than max by chi-squared
        ll.max.P <- -nll.max.bestfit
        ll.tol <-  ll.max.P - qchisq(1 - alpha/sided, 1)   # qchisq(1-0.01/1, 1) = 6.634897 

        ## function for newton.raphson() iterates on x0
        ll.q <- function(x0, data, P, shape) {
            ll <- -nll.q(data,
                         param=list(quant  = x0,
                                    shape  = shape),
                         P,
                         debug=FALSE)
        }
        
        ## determine confidence bound (P alread determined whether this was a lower or upper bound)
        ## as point where the likelihood ratio equals 11.tol
        out.nrl <- newton.raphson(f = ll.q,
                                  xguess = quant.P - quant.P/100, # to move to the left of the max likelihood
                                  ytarget = ll.tol,
                                  data   = x,
                                  P      = P,
                                  shape  = shape,
                                  tol = 1e-5, n = 1000, plot='no')
        quant.P.alpha.l <- out.nrl$root
        out.nru <- newton.raphson(f = ll.q,
                                  xguess = quant.P + quant.P/100, # to move to the right of the max likelihood
                                  ytarget = ll.tol,
                                  data   = x,
                                  P      = P,
                                  shape  = shape,
                                  tol = 1e-5, n = 1000, plot='no')
        quant.P.alpha.u <- out.nru$root

        ## 
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
                                     shape  = shape)
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
    a1 <- as.data.frame(t(unlist(params)))
    a2 <- as.data.frame(t(unlist(params.q.save[1])))
    a3 <- as.data.frame(t(unlist(params.q.save[2])))
    params.compare <- fastmerge(a1, a2)
    params.compare <- fastmerge(params.compare, a3)
    rownames(params.compare) <- c('standard fit', 'fit on quantile at 1-P', 'fit on quantile at P')

    ## collect lower and upper tolerance limits
    tol.limits <- range(tol.limits, na.rm=TRUE)
    tol.lower <- tol.limits[1]
    tol.upper <- tol.limits[2]
    tolerance <- data.frame(alpha, P, sided, tol.lower, tol.upper)
    
    return(list(params     = params,
                params.q   = params.q,
                params.compare = params.compare,
                out.bestit = out.bestfit,
                tolerance  = tolerance))
}


mle.weibull.test <- function() {
    source('modules/newton.raphson.r')
    source('modules/fastmerge.r')

    ## consider the following dataset
    x <- iris$Sepal.Width
    mean(x)  # 3.057333
    sd(x)    # 0.4358663
    P     <- 0.99  # proportion or coverage
    alpha <- 0.01
    sided <- 1

    ## mle lr approach
    out <- mle.weibull(x, 'exttol.int', alpha=alpha, P=P, sided=sided, plots=TRUE, debug=FALSE)
    params <- out$params
    print(out$params.compare)
    print(out$tolerance)

    shape <- 2
    scale <- 2.5
    out <- mle.weibull(x, list(shape=2, scale=1), alpha=alpha, P=P, sided=sided, plots=TRUE, debug=FALSE)
    print(out$params.compare)
    print(out$tolerance)
    
    ## open source tolerance limit calculation
    tol_out_weib <-  tolerance::exttol.int(x, alpha=alpha, P=P, side=1, dist="Weibull")
    print(tol_out_weib)
    
    plotspace(1,1)
    xmin <- min(0, x, out$tolerance$tol.lower, tol_out_weib$`1-sided.lower`)
    xmax <- max(0, x, out$tolerance$tol.upper, tol_out_weib$`1-sided.upper`)
    hist(x, xlim=c(xmin, xmax), freq=FALSE)
    ## add initial guess pdf
    curve(dweibull(x, shape, scale), xmin, xmax, add=TRUE,                        col='blue', lty=3)
    ## add pdf and tolerance limits from stats package
    curve(dweibull(x, shape=tol_out_weib$shape.1, scale=tol_out_weib$shape.2), xmin, xmax, add=TRUE, cex=9)
    abline(v=c(tol_out_weib$`1-sided.lower`, tol_out_weib$`1-sided.upper`),                   lty=2)
    ## add pdf and tolerance limis from MLE
    curve(dweibull(x, shape=out$params$shape, scale=out$params$scale), xmin, xmax, col='red', add=TRUE)
    abline(v=c(out$tolerance$tol.lower, out$tolerance$tol.upper),                  col='red', lty=2)
    ## add legend
    legend('topleft', 
           legend = c('initial guess', 'stats package pdf', 'stats package tol', 'MLE PDF', 'MLE tol'),
           col    = c('blue',          'black',             'black',             'red',     'red'), 
           lty    = c( 3,               1,                   2,                   1,         1))
    
}

