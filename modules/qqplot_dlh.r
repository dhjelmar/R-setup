qqplot_dlh <- function(x, dist, suppress='yes') {

    ## grab name of variable passed in as x
    xlabel <- deparse(substitute(x))

    ## first sort x from smallest to largest and drop NA entries
    x <- sort(x, na.last=NA)
    ## following also works
    ##     good <- !is.na(x)
    ##     ord  <- order(x[good])
    ##     x    <- x[good][ord]

    ## create vector of probabilities
    ## this is simply a set of num equally spaced points between 0 and 1 but not including 0 and 1
    num  <- length(x)
    pp   <- ppoints(num)

    ## determine quantiles for data
    ## quantiles are cut points dividing the range of a probability distribution
    ## into continuous intervals with equal probabilities
    if (dist == 'normal') {
        xth <- qnorm(pp, mean = mean(x), sd = sd(x))
 
    } else if (dist == 'weibull') {
        ## mass::fitdistr
        ## Distributions "beta", "cauchy", "chi-squared", "exponential", "gamma", "geometric",
        ##               "log-normal", "lognormal", "logistic", "negative binomial", "normal",
        ##               "Poisson", "t" and "weibull" are recognised, case being ignored.
        out <- MASS::fitdistr(x, densfun = 'weibull')
        shape    <- out$estimate[[1]]
        scale    <- out$estimate[[2]]
        xth <- qweibull(pp, shape = shape, scale = scale)

    } else if (dist == 'johnson') {
        library(SuppDists) # need for Johnson distribution
        jparms    <- SuppDists::JohnsonFit(x)
        xth <- SuppDists::qJohnson(pp, jparms)

    }

    ylabel <- paste('Quantiles for', dist, 'distribution', sep=' ')
    ## plot data vs. quantiles
    ## df <- data.frame(x, xth)
    ## plotfit(df, x, xth,
    ##         main   = paste('QQ Plot for', xlabel, sep=' '),
    ##         xlabel = xlabel,
    ##         ylabel = ylabel,
    ##         interval='conf', alpha=0.05, sided=2,  # confidence intervals are wrong
    ##         suppress=suppress)
    plot(x, xth,
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
    ##    z    <- qnorm(pp)
    ##    plot(z, x, type="p")
    ##    coef <- coef(rlm(x ~ z))
    ##    b    <- coef[1]
    ##    m    <- coef[2]
    ##    fit.value <- b + m*z
    ##    ## abline(b, m, col="red", lwd=1)
    ##    lines(z, fit.value, lty=1, lwd=1, col='red')
    ##    conf <- 0.95
    ##    zz   <- qnorm(1-(1-conf)/2)
    ##    SE   <- (m/dnorm(z))*sqrt(pp*(1-pp)/num)     #[WHY?]
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
    qualityTools::qqPlot(x, "normal",  col='black', main='qualityTools::qqPlot for normal')
    qqplot_dlh(x, 'normal')

    ## Weibull distribution plots
    set.seed(1)
    x <- rweibull(1E2, shape=10, scale=2)
    qualityTools::qqPlot(x, "Weibull", col='black', main='qualityTools::qqPlot for Weibull')
    qqplot_dlh(x, 'weibull')
    
}

qqplot_dlh_experiment <- function(x=NULL) {
    if (is.null(x)) {
        set.seed(20200825)
        x <- sort(rnorm(20, 10, 3))
    }
    x <- sort(x, na.last=NA)
    z <- (x - mean(x))/sd(x)                 # data convereted to standard normal
    num  <- length(x)
    pp   <- ppoints(num)
    xth <- qnorm(pp, mean=mean(x), sd=sd(x)) # theoretical          normal quantiles for p values
    ##                                       # xth <- zth * sd(x) + mean(x)
    zth <- qnorm(pp, mean=0, sd=1)           # theoretical standard normal quantiles for p values
    ##                                       # zth <- (xth - mean(x))/sd(x)

    plotspace(1,6)
    ## option 1
    car::qqPlot(x, line='quartiles', main='car::qqPlot quartiles') # other options are robust or none
                                  # quartile option is default
                                  # robust option looks like 45 degree line but on data-quantile plot
    stats::qqline(x, distribution = function(p) qnorm(p, mean = mean(x), sd = sd(x)),
                  probs = c(0.25, 0.75), qtype=7, col = 'red')
    ## option 2
    car::qqPlot(x, line='robust', main='car::qqPlot robust') # other options are robust or none
                                  # quartile option is default
                                  # robust option looks like 45 degree line but on data-quantile plot
    stats::qqline(x, distribution = function(p) qnorm(p, mean = mean(x), sd = sd(x)),
                  probs = c(0.25, 0.75), qtype=7, col = 'red')
    ## option 3
    stats::qqnorm(x, main='stats::qqnorm'); qqline(x)
    ## option 4
    qualityTools::qqPlot(x, main='qualityTools::qqPlot')
    abline(0,1, col='black', lty=2)  # shows 45 degree line is used in qqPlot
    ## option 5
    stats::qqplot(x, xth, main='stats:qqplot')
    stats::qqline(x, distribution = function(p) qnorm(p, mean = mean(x), sd = sd(x)),
                  probs = c(0.25, 0.75), qtype=7, col = 'blue') # not the same as 40 degree line
    abline(0,1, col='red')
    ## option 6
    plot(x, xth, xlab='Observed value, x', ylab='Expected normal value', main='plot(x, xth)')
    abline(0,1, col='red')

    ##------------------------------------------------------
    ## now to figure out confidence bounds
    ##------------------------------------------------------

    plotspace(1,3)

    ## website qqplot of x vs zth with confidence bands
    ## https://stats.stackexchange.com/questions/111288/confidence-bands-for-qq-line
    ## confidence bands
    ## standard error, se
    coef<-coef(rlm(x~zth))    # best fit to data
    inter<-coef[1]
    slope<-coef[2]
    SE<-(slope/dnorm(zth))*sqrt(pp*(1-pp)/num)
    ## use 2-sided 95% confidence bands
    conf<-0.95
    sided <- 2
    zz<-qnorm(1-(1-conf)/sided)   # zz = 1.96
    ## calculate upper and lower confidence bands
    xfit   <- inter + slope * zth
    xlower <- xfit - zz*SE
    xupper <- xfit + zz*SE
    ## put everything in a dataframe
    df <- tibble(x, z, xth, zth, xfit, xlower, xupper)

    ## make plots as presented in website
    plot(zth,x)                             # plot points
    abline(inter,slope,col="red",lwd=2)     # plot line
    lines(zth,xupper,lty=2,lwd=2,col="red") # plot upper band
    lines(zth,xlower,lty=2,lwd=2,col="red") # plot lowwer band

    ## dlh
    plot(x, xth, xlab='Observed value, x', ylab='Expected Value', main='Normal Distribution plot(x, xth)')
    abline(0,1, col='blue')
    lines(xfit,   xth,        col='red')
    lines(xupper, xth, lty=2, col='red')
    lines(xlower, xth, lty=2, col='red')

    ## quality tools qqplot of x vs xth
    qualityTools::qqPlot(x, main='qualityTools::qqPlot')
    abline(0,1, col='black', lty=2)  # shows 45 degree line is used in qqPlot of x vs xth
    
    return(df)
}
##qqplot_dlh_experiment()
