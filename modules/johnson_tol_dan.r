johnson_tol <- function(x, jfit='all', alpha=0.01, P=0.99, side=1, plots=FALSE, plot.ll=FALSE, breaks=NULL) {

    ## input:   x     = data set
    ##          alpha = 1 - confidence
    ##          P     = proportion = coverage (tolerance interval only)
    ##          side  = 1 = 1-sided tolerance limit
    ##                = 2 = 2-sided tolerance limit
    ##          jfit  = 'all' = uses SuppDists::JohnsonFit to determine parameters
    ##                = 'SU'  = uses ExtDist::eJohnsonSU to determine parameters
    ##                = list of user specified parameters
    ##                  e.g., jparms <- list(gamma   = -1.039,
    ##                                       delta   = 1.66,
    ##                                       xi      = 14.46,
    ##                                       lambda  = 6.95,
    ##                                       type    = 'SU',
    ##                       johnson_tol(mtcars$mpg, jfit=jparms)
    ##
    ## output:  jparms = Johnson fit parameters used in tolerance limit calculation
    ##          xtol_lower = lower bound tolerance limit
    ##          xtol_upper = upper bound tolerance limit
    
    coverage <- P

    ##--------------------------------------------------------
    ## initialize tolerance limits
    xntol_lower  = NA
    xntol_upper  = NA
    xtol_lower   = NA
    xtol_upper   = NA

    
    ##--------------------------------------------------------
    ## determine Johnson parameters
    
    if (jfit[1] == 'all') {
        ## let R figure out which Johnson distribution fits best
        jparms  <- SuppDists::JohnsonFit(x)
    } else if (jfit[1] == 'SU') {
        ## force the Johnson SU distribution
        jparms.out <- ExtDist::eJohnsonSU(x)
        jparms <- list(gamma   = jparms.out$gamma,
                       delta   = jparms.out$delta,
                       xi      = jparms.out$xi,
                       lambda  = jparms.out$lambda,
                       type    = 'SU')
    } else {
        ## use Johnson parameters specified in jfit
        ## needs to be in same list format as created by SuppDists::JohnsonFit
        jparms <- jfit
    }
    gamma  <- jparms$gamma
    delta  <- jparms$delta
    xi     <- jparms$xi
    lambda <- jparms$lambda
    type   <- jparms$type

    
    ##--------------------------------------------------------
    ## Calculate tolerance limits if JohnsonSU fit
    jparms.optim <- NA
    jparms.q     <- NA
    if (type == 'SU') {

        ##----------------------
        ## first check to see what parameters an MLE fit gives

        ## define function to return the nnl (negative log likelihood) for
        ## a given data set and JohnsonSU parameters
        nll.johnsonsu <- function(x, jparms, debug=FALSE) {
            ## x is the name of the variable containing the data to be fitted
            ## if access nnl.johnsonsu with optim, cannot access parameters within
            ## jparms (e.g., jparms$gamma) with $ because jparms is passed in as atomic vector
            gamma <- jparms[[1]]  
            delta <- jparms[[2]]
            xi    <- jparms[[3]]
            lambda <- jparms[[4]]
            ## PDF for Johnson SU
            pdf <- delta /( lambda * sqrt(2 * pi)   ) *
                1 / sqrt(1 +            ( (x-xi)/lambda )^2)  *
                exp( -0.5*(gamma + delta * asinh( (x-xi)/lambda ))^2 )
            -sum(log(pdf))
            ## the following is equivalent to the above as long as it picks Johnson SU
            ## type <- 'SU'
            ## -sum(SuppDists::dJohnson(x, parms=c(gamma, delta, xi, lambda, type), log=TRUE))
            if (isFALSE(debug)) {
                return( -sum(log(pdf)) )
            } else {
                return(list(gamma=gamma, delta=delta, xi=xi, lambda=lambda, pdf=pdf))
            }
        }

        ## find the JohnsonSU parameters that minimize the nnl
        out.optim <- optim(par = jparms[1:4], # parameters to be optimized and initial guess
                           fn  = nll.johnsonsu,     # NNL function
                           x   = x,                 # data
                           control = list(parscale = jparms[1:4]) # not clear this was needed
                           )
        jparms.optim <- as.list(out.optim$par)
        jparms.optim$type <- 'SU'

        ##----------------------
        ## define function to return the nnl (negative log likelihood) for
        ## an alternate set and JohnsonSU parameters where replace gamma with quantile, i.e.:
        ##      quantile, delta, xi, lambda
        ## instead of:
        ##      gamma,    delta, xi, lambda
        nll.johnsonsu.quant <- function(x, coverage, jparms, debug=FALSE) {
            ## x is the name of the variable containing the data to be fitted
            ## coverage is the percent coverage
            ## if use optim, cannot access as jparms$gamma because
            ## jparms is passed in as atomic vector
            quant  <- jparms[[1]]  # replaced gamma with quant as a parameter
            delta  <- jparms[[2]]
            xi     <- jparms[[3]]
            lambda <- jparms[[4]]
            ## write gamma as a function of quant, delta, xi and lambda
            gamma <- qnorm(coverage) - delta * asinh( (quant-xi)/lambda )
            ## PDF for Johnson SU
            pdf <- delta /( lambda * sqrt(2 * pi)   ) *
                1 / sqrt(1 +            ( (x-xi)/lambda )^2)  *
                exp( -0.5*(gamma + delta * asinh( (x-xi)/lambda ))^2 )
            if (isFALSE(debug)) {
                return( -sum(log(pdf)) )
            } else {
                return(list(quant=quant, gamma=gamma, delta=delta, xi=xi, lambda=lambda, pdf=pdf))
            }
        }

        ## find the alternate JohnsonSU parameters
        quant.coverage <- ExtDist::qJohnsonSU(coverage, params = jparms.optim)
        quant.parms <- c(quant=quant.coverage, jparms.optim[2:4])
        out.optim.quant <- optim(par = quant.parms, # parameters to be optimized and initial guess
                                 fn  = nll.johnsonsu.quant,     # NLL function
                                 x   = x,                       # data
                                 coverage   = 0.99,             # percent coverage
                                 control = list(parscale = quant.parms), # not clear this was needed
                                 hessian = TRUE
                                 )
        estimates      <- as.numeric(out.optim.quant$par)
        standard.error <- as.numeric( sqrt(diag(solve(out.optim.quant$hessian))) )
        coef           <- data.frame(estimates, standard.error)
        rownames(coef) <- names(out.optim.quant$par)

        jparms.optim.quant <- as.list(out.optim.quant$par)

        ## write equivalent Johnson SU parameters in terms of gamma
        gamma <- qnorm(coverage) - jparms.optim.quant$delta *
            asinh( (jparms.optim.quant$quant-jparms.optim.quant$xi)/jparms.optim.quant$lambda )
        jparms.q <- list(gamma  = gamma,
                         delta  = jparms.optim.quant$delta,
                         xi     = jparms.optim.quant$xi,
                         lambda = jparms.optim.quant$lambda,
                         type   = 'SU',
                         quant  = jparms.optim.quant$quant)
        
        ##----------------------
        ## calculate the tolerance limit
        ## factor <- qnorm(1 - alpha/side)   # alpha = 0.01, sided=1, z = 2.326
        ## dof    <- length(x) - 4           # 4 independent fitting parameters in Johnson SU
        ## student.t <- qt(1 - alpha/side, dof)  # 2.3326 for dof=598
        ## factor <- student.t * sqrt(1 + 1/dof) # 2.3345
        ## factor <- ExtDist::qJohnsonSU(1 - alpha/side, params = jparms.q[1:5])
        ## 
        ##         jparms.dan <- list(gamma  = 3.0202,
        ##                            delta  = 1.4506,
        ##                            xi     = 1.0148,
        ##                            lambda = 00531,
        ##                            type   = 'SU')
        ##         factor <- ExtDist::qJohnsonSU(1 - alpha/side, params = jparms.dan) # 38
        ## 
        ## factor <- 2.6645    # not sure where 2.6645 comes from
        ## xtol_lower <- jparms.q$quant - factor * standard.error[1]
        ## xtol_upper <- jparms.q$quant + factor * standard.error[1]

        ## following approach here:
        ## https://personal.psu.edu/abs12/stat504/Lecture/lec3_4up.pdf

        if (isTRUE(plot.ll)) {
            ## plot the likelihood as a function of the quantile
            quant.min <- estimates[1] - 2 * standard.error[1]
            quant.max <- estimates[1] + 2 * standard.error[1]
            quant     <- seq(quant.min, quant.max, length.out=100)
            ll <- NA
            for (i in 1:100) {
                ## vary quant parameter to see impact on likelihood
                jparms.vary.q <- list(quant  = quant[i],
                                      delta  = jparms.optim.quant$delta,
                                      xi     = jparms.optim.quant$xi,
                                      lambda = jparms.optim.quant$lambda)
                ll[i] <- -nll.johnsonsu.quant(x, coverage, jparms=jparms.vary.q, debug=FALSE)
            }
            plot(quant, ll, xlab='quantile', ylab='Log Likelihood')
            abline(v=estimates[1])
        }

        ## find the likelihood associated with the best estimate
        ll.max <- -nll.johnsonsu.quant(x,
                                       coverage,
                                       jparms=list(quant = estimates[1],
                                                   delta  = jparms.optim.quant$delta,
                                                   xi     = jparms.optim.quant$xi,
                                                   lambda = jparms.optim.quant$lambda),
                                       debug=FALSE)

        ## confidene limit is defined at likelihood that is lower than max by chi-squared
        ll.tol <-  ll.max - qchisq(1 - alpha/side, 1)   # qchisq(1-0.01/1, 1) = 6.634897 
        if (isTRUE(plot.ll)) {
            abline(h=ll.tol)
        }

        ## find the quantiles associaed with the likelihood tolerance
        zero <- function(x0, data, coverage, delta, xi, lambda, ll.tol) {
            ll <- -nll.johnsonsu.quant(data,
                                       coverage,
                                       jparms=list(quant  = x0,
                                                   delta  = delta,
                                                   xi     = xi, 
                                                   lambda = lambda),
                                       debug=FALSE)
            zero <- ll - ll.tol
        } 
        out.nrl <- newton.raphson(f = zero,
                                  data = x,
                                  coverage = coverage,
                                  delta  = jparms.optim.quant$delta,
                                  xi     = jparms.optim.quant$xi,
                                  lambda = jparms.optim.quant$lambda,
                                  ll.tol = ll.tol,
                                  xguess = estimates[1] - standard.error[1],  # xguess is used for 1st param in f
                                  tol = 1e-5, n = 1000, plot='no')
        quant.tol.lower <- out.nrl$root
        out.nru <- newton.raphson(f = zero,
                                  data = x,
                                  coverage = coverage,
                                  delta  = jparms.optim.quant$delta,
                                  xi     = jparms.optim.quant$xi,
                                  lambda = jparms.optim.quant$lambda,
                                  ll.tol = ll.tol,
                                  xguess = estimates[1] + standard.error[1],  # xguess is used for 1st param in f
                                  tol = 1e-5, n = 1000, plot='no')
        quant.tol.upper <- out.nru$root
        if (isTRUE(plot.ll)) {
            points(estimates[1]   , ll.max, col='blue', pch=16, cex=2)
            points(quant.tol.lower, ll.tol, col='red' , pch=16, cex=2)
            points(quant.tol.upper, ll.tol, col='red' , pch=16, cex=2)
        }
        
        if (isTRUE(plots)) {
            ## plot fit from R packages vs my MLE fit vs. my MLE fit on quantiles
            plotspace(2,2)
            hist(x, freq=FALSE, main='Johnson SU Fits', xlim=range(quant.tol.lower, x, quant.tol.upper))
            curve(SuppDists::dJohnson(x, jparms)       , min(x, quant.tol.lower), max(x, quant.tol.upper),
                  col='black', lty=1, add=TRUE)
            curve(SuppDists::dJohnson(x, jparms.optim) , min(x, quant.tol.lower), max(x, quant.tol.upper),
                  col='red', lty=1, add=TRUE, lwd=2)
            curve(SuppDists::dJohnson(x, jparms.q[1:5]), min(x, quant.tol.lower), max(x, quant.tol.upper),
                  col='blue', lty=2, add=TRUE)
            abline(v=quant.tol.upper, col='red', lwd=2)
            legend("topright", c("SuppDists or ExtDist", "MLE", "MLE fit on quantile", "upper tolerance limit"),
                   col=c("black", "red", "blue", 'red'), lty=c(1,1,2,2))
            qqplot_nwj(x, type='j', jfit=jparms       , main='ExtDist or SuppDists auto fit')
            qqplot_nwj(x, type='j', jfit=jparms.optim , main='MLE fit with optim')
            qqplot_nwj(x, type='j', jfit=jparms.q[1:5], main='MLE fit with optim on quantile')
            
        }

    }

    return(list(side         = side, 
                alpha        = alpha,
                coverage     = coverage,
                jfit         = jfit,
                jparms       = jparms,
                jparms.optim = jparms.optim,
                jparms.q     = jparms.q,
                x            = x,
                xtol_lower   = quant.tol.lower, 
                xtol_upper   = quant.tol.upper))
}

## ## test

## johnson_tol(mtcars$mpg)

## out <- johnson_tol(mtcars$mpg)
## upper_tolerance_limit_john <- out$xtol_upper
## upper_tolerance_limit_john
## jparms <- out$jparms

## johnson_tol(mtcars$mpg, jfit=jparms)

## johnson_tol(mtcars$mpg, plots=TRUE)

## x <- mtcars$mpg
## johnson_tol(x, jfit='all', alpha=0.01, P=0.99, side=1, plots=FALSE, breaks=NULL)
