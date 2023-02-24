qqplot_nwj_xonly <- function(x, type='nwj', wfit='mle', jfit='mle', mainadder=NULL) {
    ## creates normal, Weibull and/or Johnson SU qq plots

    ## input: x    = vector of data
    
    ## make room for 1, 2 or 3 plots depending on length of string 'type'
    nplots <- nchar(type)
    if (nplots != 1) par(mfrow=c(1, nplots))

    nparms <- NULL
    wparms <- NULL
    jparms <- NULL
    if (grepl('n', type)) {
        ## make normal QQ plot
        main <- paste('Normal QQ Plot', mainadder, sep=" ")
        ## qualityTools package not being maintained so replaced here
        ## nparms <- qualityTools::qqPlot(x, "normal",  col='black', main=main)
        nparms <- list(mean=mean(x), sd=sd(x))
        x <- sort(x, na.last=NA)
        xtheoretical <- qnorm(ppoints(length(x)), mean = mean(x), sd = sd(x))
        main <- paste('Normal QQ Plot', mainadder, sep=" ")
        plot(x, xtheoretical, xlab='Observed value, x', ylab='Expected Value', main=main)
        abline(0,1, col='red')
        ## car::qqPlot(x, envelope=0.95, xlab='Expected Quantiles', ylab='Observed Values')
        
    }
    
    if (grepl('w', type) & min(x)>0) {
        ## obtain Weibull parameters
        wparms <- wfit
        if (wfit[1] == 'mle') {
            out <- mle.weibull(x)
            wparms <- out$parms
        } else if (wfit[1] == 'exttol') {
            ## tolerance package seems more robust than qualityTools
            ##     wparms <- qualityTools::qqPlot(x, "Weibull", col='black', main=main)
            ## also, tolerance packages is used in hist_nwj so more consistent
            out     <- tolerance::exttol.int(x, alpha=0.05, P=0.95, side=1, dist='Weibull')
            shape   <- out$'shape.1'
            scale   <- out$'shape.2'
            wparms  <- list(shape=shape, scale=scale)
        }
        
        ## make QQ plot
        main <- paste('Weibull QQ Plot', mainadder, sep=" ")
        ## qualityTools has nice confidence bounds on QQ plot but no ability for censored data
        ## qualityTools::qqPlot(x, "Weibull", col='black', main=main, start=wparms)
        ## sort data and add censored data, if needed, to prepare for QQ plot
        ## na.last = NA removes missing values
        ##         = TRUE puts missing values last
        ##         = FALSE puts missing values first
        x <- sort(x, na.last=NA)
        xtheoretical <- stats::qweibull(ppoints(length(x)), shape=wparms$shape, scale=wparms$scale)
        plot(x, xtheoretical, xlab='Observed value, x', ylab='Expected Value', main=main)
        abline(0,1, col='red')

        
    } else if (grepl('w', type)) {
        ## some values are negative or zero so Weibull is not appropriate
        cat('Weibull plot not made because not appropriate; some values not > 0\n')
    }

    if (grepl('j', type)) {        
        ## obtain Johnson SU parameter estimates
        jparms <- jfit
        if (jfit[1] == 'mle') {
            out <- mle.johnsonsu(x)
            jparms <- out$parms
        }

        ## make QQ plot
        main <- paste('Johnson SU QQ Plot', mainadder, sep=" ")
        ## sort data and add censored data, if needed, to prepare for QQ plot
        ## na.last = NA removes missing values
        ##         = TRUE puts missing values last
        ##         = FALSE puts missing values first
        x <- sort(x, na.last=NA)
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

    x <- mtcars$mpg
    set.seed(15)
    x <- rweibull(1000, shape=55, scale = 2)
    x <- sort(x)
    qqplot_nwj(x)

    par(mfrow=c(2,2))

    ## comparisons for normal qq plot using same method as for Johnson SU qqplot
    ## qualityTools::qqPlot(x, "normal",  col='black')
    car::qqPlot(x, envelope=0.95, xlab='Expected Quantile', ylab='Observed Value', 
                main='Normal Distribution QQ Plot using car::qqPlot()')
    xtheoretical <- qnorm(ppoints(length(x)), mean = mean(x), sd = sd(x))
    plot(x, xtheoretical, xlab='Observed value, x', ylab='Expected Value',
         main='Normal Distribution QQ Plot using qnorm(ppoints())')
    abline(0,1, col='red')
    
    ## comparisons for Weibull qq plot using same method as for Johnson SU qqplot
    ## qualityTools::qqPlot(x, "Weibull",  col='black')
    out     <- tolerance::exttol.int(x, alpha=0.05, P=0.95, side=1, dist='Weibull')
    shape   <- out$'shape.1'
    scale   <- out$'shape.2'
    car::qqPlot(x, distribution='weibull', shape=shape, scale=scale,
                envelope=0.95, xlab='Expected Quantile', ylab='Observed Value', 
                main='Weibull Distribution QQ Plot using car::qqPlot()')
    tol_out <-  exttol.int(x, alpha =0.1, P=0.99, side=1)
    shape <- tol_out$'shape.1'
    scale <- tol_out$'shape.2'
    xtheoretical <- qweibull(ppoints(length(x)), shape = shape, scale = scale)
    plot(x, xtheoretical, xlab='Observed value, x', ylab='Expected Value',
         main='Weibull Distribution QQ Plot using qnorm(ppoints())')
    abline(0,1, col='red')

}
