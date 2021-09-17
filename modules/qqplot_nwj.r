qqplot_nwj <- function(x, type='nwj', jfit='all', mainadder=NULL) {
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
        wparms <- qualityTools::qqPlot(x, "Weibull", col='black', main=main)
    } else if (grepl('w', type)) {
        ## some values are negative or zero so Weibull is not appropriate
        cat('Weibull plot not made because not appropriate; some values not > 0\n')
    }

    if (grepl('j', type)) {        
        ## obtain Johnson parameter estimates
        x <- sort(x, na.last=NA)
        if (jfit == 'all') {
            ## let R figure out which Johnson distribution fits best
            jparms  <- SuppDists::JohnsonFit(x)
            main <- paste('Johnson QQ Plot', mainadder, '; Type=', jparms$type, sep=" ")
        } else if (jfit == 'SU') {
            ## force the Johnson SU distribution
            jparms.out <- ExtDist::eJohnsonSU(x)
            jparms <- list(gamma   = jparms.out$gamma,
                           delta   = jparms.out$delta,
                           xi      = jparms.out$xi,
                           lambda  = jparms.out$lambda,
                           type <- 'SU')
            main <- paste('JohnsonSU QQ Plot', mainadder, sep=" ")
        } else if (jfit == 'SB') {
            ## force the Johnson SB distribution
            jparms <- ExtDist::eJohnsonSU(x)
            jparms$type <- 'SB'
            main <- paste('JohnsonSB QQ Plot', mainadder, sep=" ")
        } else {
            ## use Johnson parameters specified in jfit
            ## needs to be in same list format as created by SuppDists::JohnsonFit
            jparms <- jfit
            main <- paste('User Specified Johnson QQ Plot', mainadder, sep=" ")
        }

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

    if (grepl('k', type)) {
        ## alternate approach for Johnson QQ plot to get confidence intervals

        ## transpose Johnson distributed data to normal
        out <- johnson_tol(x, jfit='all', alpha=0.01, P=0.99, side=1, plots='no')
        xn <- out$xn
        
        ## now make normal QQ plot using 
        main <- paste('Johnson Transposed to Normal QQ Plot', mainadder, sep=" ")
        nparms <- qualityTools::qqPlot(xn, "normal",  col='black', main=main)        
    }
    
    return(list(nparms=nparms, wparms=wparms, jparms=jparms))
    
}
## qqplot_nwj(mtcars$mpg)


## ## comparisons for normal qq plot
## par(mfrow=c(1,2))
## qualityTools::qqPlot(x, "normal",  col='black')
## quantiles_norm <- qnorm(ppoints(length(x)), mean = mean(x), sd = sd(x))
## stats::qqplot(
##            x = x,
##            xlab = "Data",
##            y = quantiles_norm,
##            ylab = 'Quantiles from Normal Distribution',
##            main = expression('Q-Q plot for Normal Distribution'),
##            col='black')
## qqline(x, distribution = function(p) qnorm(p, mean=mean(x), sd=sd(x)), col=2)
## 
## ## comparisons for Weibull qq plot
## par(mfrow=c(1,2))
## qualityTools::qqPlot(x, "Weibull",  col='black')
## tol_out <-  exttol.int(x, alpha =0.1, P=0.99, side=1)
## shape <- tol_out$'shape.1'
## scale <- tol_out$'shape.2'
## quantiles_weib <- qweibull(ppoints(length(x)), shape = shape, scale = scale)
## stats::qqplot(
##            x = x,
##            xlab = "Data",
##            y = quantiles_weib,
##            ylab = 'Quantiles from Weibull Distribution',
##            main = expression('Q-Q plot for Weibull Distribution'),
##            col='black')
## qqline(x, distribution = function(p) qweibull(p, shape = shape, scale = scale), col=2)

## ## test method of using qJohnson by using same method for qnorm
## plotspace(1,2)
## x <- sort(mtcars$mpg)
## ## by hand
## probabilities <- ppoints(length(x))
## xtheoretical  <- qnorm(probabilities, mean=mean(x), sd=sd(x))
## plot(x, xtheoretical, xlab='Observed value, x', ylab='Expected Value')
## abline(0,1,col='red')
## ## using Rfunction
## qualityTools::qqPlot(x)
