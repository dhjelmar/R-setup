qqplot_nwj <- function(x, type='nw') {
    ## creates side by side, normal and Weibull qq plots
    if (type == 'nw') {
        par(mfrow=c(1,2))
        qualityTools::qqPlot(x, "normal",  col='black')
        qualityTools::qqPlot(x, "Weibull", col='black')
    } else {
        par(mfrow=c(1,3))
        qualityTools::qqPlot(x, "normal",  col='black')
        qualityTools::qqPlot(x, "Weibull", col='black')
        library(SuppDists) # need for Johnson distribution
        jparms  <- JohnsonFit(x)
        quantiles_john <- qJohnson(ppoints(length(x)), jparms)
        stats::qqplot(
            x = x,
            xlab = "Data",
            y = quantiles_john,
            ylab = 'Quantiles from "Johnson" Distribution',
            main = expression('Q-Q plot for "Johnson" Distribution'),
            col='black')
        qqline(x, distribution = function(p) qJohnson(p, jparms), col=2)
    }
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
