qqplot_dlh <- function(x, dist, suppress='yes') {

    ## grab name of variable passed in as x
    xlabel <- deparse(substitute(x))

    ## first sort x so it increases
##    x <- sort(x)
    good <- !is.na(x)
    ord  <- order(x[good])
    x    <- x[good][ord]

    ## determine quantiles for data
    num  <- length(x)
    P    <- ppoints(num)
    if (dist == 'normal') {
        quantiles <- qnorm(P, mean = mean(x), sd = sd(x))
        
    } else if (dist == 'weibull') {
        ## mass::fitdistr
        ## Distributions "beta", "cauchy", "chi-squared", "exponential", "gamma", "geometric",
        ##               "log-normal", "lognormal", "logistic", "negative binomial", "normal",
        ##               "Poisson", "t" and "weibull" are recognised, case being ignored.
        out <- MASS::fitdistr(x, densfun = 'weibull')
        shape    <- out$estimate[[1]]
        scale    <- out$estimate[[2]]
        quantiles <- qweibull(P, shape = shape, scale = scale)

    } else if (dist == 'johnson') {
        library(SuppDists) # need for Johnson distribution
        jparms    <- SuppDists::JohnsonFit(x)
        quantiles <- SuppDists::qJohnson(P, jparms)

    }

    ylabel <- paste('Quantiles for', dist, 'distribution', sep=' ')
    ## plot data vs. quantiles
    ## df <- data.frame(x, quantiles)
    ## plotfit(df, x, quantiles,
    ##         main   = paste('QQ Plot for', xlabel, sep=' '),
    ##         xlabel = xlabel,
    ##         ylabel = ylabel,
    ##         interval='conf', alpha=0.05, sided=2,  # confidence intervals are wrong
    ##         suppress=suppress)
    plot(x, quantiles,
         main = paste('QQ Plot for', xlabel, sep=' '),
         xlab = xlabel,
         ylab = ylabel)
    abline(a=0, b=1, col='red')  # theoretical line fit

    ##    ## students-t for dataset
    ##    num <- length(x)
    ##    conf <- 0.05
    ##    t  <- qt(1-conf/2, num-1)
    ##    uncert    <- t * sd(x) / sqrt(num)

    ##    ## confidence limits for plot
    ##    ## following makes confidence limits but they do not match standard R and following is only for normal
    ##    ## https://stats.stackexchange.com/questions/111288/confidence-bands-for-qq-line
    ##    z    <- qnorm(P)
    ##    plot(z, x, type="p")
    ##    coef <- coef(rlm(x ~ z))
    ##    b    <- coef[1]
    ##    m    <- coef[2]
    ##    fit.value <- b + m*z
    ##    ## abline(b, m, col="red", lwd=1)
    ##    lines(z, fit.value, lty=1, lwd=1, col='red')
    ##    conf <- 0.95
    ##    zz   <- qnorm(1-(1-conf)/2)
    ##    SE   <- (m/dnorm(z))*sqrt(P*(1-P)/num)     #[WHY?]
    ##    upper <- fit.value+zz*SE
    ##    lower <- fit.value-zz*SE
    ##    lines(z,upper, lty=2, lwd=1, col="red")
    ##    lines(z,lower, lty=2, lwd=1, col="red")

    ## lines(x, lquant, col='black', lty='dashed')  # lower confidence line
    ## lines(x, uquant, col='black', lty='dashed')  # upper confidence line

}



qqplot_dlh_test <- function() {
    ## test qqplot_dlh against qualityTools::qqplot

    plotspace(2,2)

    ## normal distribution plots
    set.seed(1)
    x <- rnorm(n=1E2, mean=10, sd=1)
    qualityTools::qqPlot(x, "normal",  col='black')
    qqplot_dlh(x, 'normal')

    ## Weibull distribution plots
    set.seed(1)
    x <- rweibull(1E2, shape=10, scale=2)
    qualityTools::qqPlot(x, "Weibull", col='black')
    qqplot_dlh(x, 'weibull')
    
}
