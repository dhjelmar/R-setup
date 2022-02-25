qqplot_nwj <- function(x, type='nwj', mle=TRUE, wfit=NULL, jfit='auto', mainadder=NULL) {
    ## creates side by side, normal, Weibull and/or Johnson qq plots

    ## make room for 1, 2 or 3 plots depending on length of string 'type'
    nplots <- nchar(type)
    if (nplots != 1) par(mfrow=c(1, nplots))

    nparms <- NULL
    wparms <- NULL
    jparms <- NULL
    if (grepl('n', type)) {
        ## make normal QQ plot
        main <- paste('Normal QQ Plot', mainadder, sep=" ")
        nparms <- qualityTools::qqPlot(x, "normal",  col='black', main=main)
    }
    
    if (grepl('w', type) & min(x)>0) {        
        ## make Weibull QQ plot
        main <- paste('Weibull QQ Plot', mainadder, sep=" ")
        ## tolerance package seems more robust than qualityTools
        ##     wparms <- qualityTools::qqPlot(x, "Weibull", col='black', main=main)
        ## also, tolerance packages is used in hist_nwj so more consistent
        out     <- tolerance::exttol.int(x, alpha=0.05, P=0.95, side=1, dist='Weibull')
        shape   <- out$'shape.1'
        scale   <- out$'shape.2'
        if (!is.null(wfit[1])) {
            shape <- wfit[[1]]
            scale <- wfit[[2]]
        }
        if (isFALSE(mle)) {
            out <- mle.weibull(x, list(shape=shape, scale=scale), fit.only=TRUE,
                               alpha=alpha, P=P, sided=1,
                               plots=FALSE, debug=FALSE)
            shape <- out$params$shape
            scale <- out$params$scale
        }
        wparms  <- list(shape=shape, scale=scale)
        qualityTools::qqPlot(x, "Weibull", col='black', main=main, start=wparms)
    } else if (grepl('w', type)) {
        ## some values are negative or zero so Weibull is not appropriate
        cat('Weibull plot not made because not appropriate; some values not > 0\n')
    }

    if (grepl('j', type)) {        
        ## obtain Johnson parameter estimates
        x <- sort(x, na.last=NA)
##         if (jfit[1] == 'SuppDists') {
##             ## let R figure out which Johnson distribution fits best
##             jparms  <- SuppDists::JohnsonFit(x)
##             main <- paste('Johnson QQ Plot', mainadder, '; Type=', jparms$type, sep=" ")
##         } else if (jfit[1] == 'ExtDist') {
##             ## force the Johnson SU distribution
##             jparms.out <- ExtDist::eJohnsonSU(x)
##             jparms <- list(gamma   = jparms.out$gamma,
##                            delta   = jparms.out$delta,
##                            xi      = jparms.out$xi,
##                            lambda  = jparms.out$lambda,
##                            type <- 'SU')
##             main <- paste('JohnsonSU QQ Plot', mainadder, sep=" ")
##         } else {
##             ## use Johnson parameters specified in jfit
##             ## needs to be in same list format as created by SuppDists::JohnsonFit
##             jparms <- jfit
##             main <- paste('User Specified Johnson QQ Plot', mainadder, sep=" ")
##         }
        ##         ## refit using MLE if specifed
        ##         if (isTRUE(mle)) {
        out <- mle.johnsonsu(x, jfit, fit.only=TRUE,
                             alpha=alpha, P=P, sided=1,
                             plots=FALSE, debug=FALSE)
        jparms <- out$params
        main <- paste('MLE Johnson QQ Plot', mainadder, sep=" ")
        ##         }

        
        ## make Johnson QQ plot
        xtheoretical <- SuppDists::qJohnson(ppoints(length(x)), jparms)
        plot(x, xtheoretical, xlab='Observed value, x', ylab='Expected Value', main=main)
        abline(0,1, col='red')

        ## quantiles_john <- qJohnson(ppoints(length(x)), jparms)
        ## stats::qqplot(
        ##     x = x,
        ##     xlab = "Data",
        ##     y = quantiles_john,
        ##     ylab = 'Quantiles from "Johnson" Distribution',
        ##     main = expression('Q-Q plot for "Johnson" Distribution'),
        ##     col='black')
        ## stats::qqline(x, distribution = function(p) qJohnson(p, jparms), col='red')
        ## abline(0,1, col='blue', lty=2)
        ## legend("bottomright",col=c('red', 'blue'), lty=c(1,2), legend=c('qqline', '45 degree line'))
    }
    
    return(list(nparms=nparms, wparms=wparms, jparms=jparms))
    
}



##-----------------------------------------------------------------------------
qqplot_nwj_tests <- function() {

    qqplot_nwj(mtcars$mpg)

    x <- sort(mtcars$mpg)

    ## comparisons for normal qq plot using same method as for Johnson qqplot
    par(mfrow=c(1,2))
    qualityTools::qqPlot(x, "normal",  col='black')
    xtheoretical <- qnorm(ppoints(length(x)), mean = mean(x), sd = sd(x))
    plot(x, xtheoretical, xlab='Observed value, x', ylab='Expected Value')
    abline(0,1, col='red')
    
    ## comparisons for Weibull qq plot using same method as for Johnson qqplot
    par(mfrow=c(1,2))
    qualityTools::qqPlot(x, "Weibull",  col='black')
    tol_out <-  exttol.int(x, alpha =0.1, P=0.99, side=1)
    shape <- tol_out$'shape.1'
    scale <- tol_out$'shape.2'
    xtheoretical <- qweibull(ppoints(length(x)), shape = shape, scale = scale)
    plot(x, xtheoretical, xlab='Observed value, x', ylab='Expected Value')
    abline(0,1, col='red')

}
