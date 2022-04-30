qqplot_nwj <- function(x=NA, xcen=NA, type='nwj', wfit='mle', jfit='mle', mainadder=NULL) {
    ## creates normal, Weibull and/or Johnson qq plots
    
    ## input: x    = vector of data
    ##        xcen = dataframe of censored data

    if (!is.data.frame(xcen)) {
        ## no censored data were supplied, so use standard qqplot_nwj function
        out <- qqplot_nwj_xonly(x, type, wfit, jfit, mainadder)
        return(out)

    } else {
        ## redefine xcen with specific names and only 2 columns if more were included
        x.low  <- xcen[[1]]
        x.high <- xcen[[2]]
        xcen   <- data.frame(x.low, x.high)
        
        if (!is.na(x[1])) {
            ## known values also supplied, so add to censored values
            xcen2 <- data.frame(x.low=x, x.high=x)
            xcen <- rbind(xcen, xcen2)
        }
    }

    ## Kaplan and Meier method
    ## see https://pdixon.stat.iastate.edu/stat505/Chapter%2011.pdf, p. 4
    ## method seems to assume that known values provide boundaries to censored values
    ## That is not generally true unless everything is evaluated at every inspection time.
    ## As such, the algorithm in the paper does not directly apply.

    ## instead, try: order observations based on x.high from low to high
    ##               

    
    ## sort data from low to high
    xcen <- xcen[order(xcen$x.low),]
    xcen <- xcen[order(xcen$x.high),]
    
    ## consider: N = observations
    ##           J = number of unique known (uncensored) observations
    N  <- nrow(xcen)
    ## xu <- na.omit(unique(xcen[xcen3$x.low == xcen$x.high,]))
    ## J  <- nrow(xu)
    ## if (J < 1) {
    ##     cat('ERROR: No known values\n')
    ##     return()
    ## }

    ## create dataframe of the unique x values regardless of whether start or end of an interval
    xcen <- as_tibble(xcen)
    nobs <- nrow(xcen)
    x.all <- c(xcen$x.low, xcen$x.high)
    x <- unique(sort(x.all))
    xu <- as.data.frame(x=x)
    
    xu$fail <- NA
    xu$surv <- NA
    xu$fail.ratio <- NA
    xu$surv.ratio <- NA
    for (i in 1:nrow(xu)) {
      
        ## count values that are less than or equal to x
        xu$fail[i] <- sum(xcen$x.high <= x[i], na.rm=TRUE)
        
        ## count values that are greater than x
        xu$surv[i] <- sum(xcen$x.low > x[i], na.rm=TRUE)
        
        ## observed ratios
        xu$fail.ratio[i] <- xu$fail[i] / nobs
        xu$surv.ratio[i] <- xu$surv[i] / nobs
        
    }
    ## cumulative survival probability                         # useless?
    xu$fail.cumsurv <- cumprod(1 - xu$fail.ratio)
    xu$surv.cumsurv <- cumprod(1 - xu$surv.ratio)
    
    ## johnson su fit
    out.fit <- mle.johnsonsu(xcen=xcen)
    params <- out.fit$parms
    gamma <- params$gamma
    delta <- params$delta
    xi    <- params$xi
    lambda <- params$lambda
    
    ## determine CDF from predicted fit at each xu value
    xu$cdf <- stats::pnorm(gamma+delta*asinh((xu$x-xi)/lambda))
    
    ## plot cumulative distribution for data vs. fit
    plotspace(1,2)
    plot(xu$surv.cumsurv, xu$cdf)                               # useless?
    plot(xu$fail.ratio, xu$cdf)
    abline(0,1, col='red')
    
    return(xu)
}

qqplot_nwj_test <- function() {
    source('setup.r')
    x <- iris$Sepal.Width
    xnum <- 20
    x.low  <- NA
    x.high <- 2.0
    xcen1 <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    x.low  <- 3.7
    x.high <- NA
    xcen2 <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    x.low  <- 2.6
    x.high <- 3.5
    xcen3 <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    x.low  <- 3.1234
    x.high <- 3.1234
    xcen4 <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    xcen <- rbind(xcen1, xcen2, xcen3, xcen4)
    out.qq <- qqplot_nwj(x, xcen, type='j')
    return(out.qq)
}
