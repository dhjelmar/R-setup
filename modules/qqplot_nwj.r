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
        xall   <- data.frame(x.low, x.high)
        
        if (!is.na(x[1])) {
            ## known values also supplied, so add to censored values
            x2 <- data.frame(x.low=x, x.high=x)
            xall <- rbind(x2, xall)
        }

        ## now that everything is in xcen, separate again into known, x, and censored, xcen
        xcen.na  <- xall[ is.na(rowSums(xall)),]                  # censored rows with NA
        xall.val <- xall[!is.na(rowSums(xall)),]                  # rows w/o NA
        x <- xall.val[xall.val[[1]] == xall.val[[2]],][[1]]       # known values from xall
        x <- as.numeric(na.omit(x))                               # new set of known values
        xcen.lowhigh <- xall.val[xall.val[[1]] != xall.val[[2]],] # censored rows with no NA
        xcen <- rbind(xcen.lowhigh, xcen.na)                      # new set of censored rows
        
        ## create new xall with another column that identifies if known or censored
        xall.x    <- data.frame(x.low=x         , x.high=x          , type='known')
        if (nrow(xcen) > 0) {
             xall.xcen <- data.frame(x.low=xcen$x.low, x.high=xcen$x.high, type='censored')
             xall <- rbind(xall.x, xall.xcen)
        } else {
             xall <- xall.x
        }

    }

    ## Kaplan and Meier method
    ## see https://pdixon.stat.iastate.edu/stat505/Chapter%2011.pdf, p. 7

    ## sort data from low to high and renumber
    ## for the same x.high value, sort so censored observations are below known observations
    xall <- xall[order(xall$x.low, na.last=FALSE),]
    xall <- xall[order(xall$x.high),]
    rownames(xall) <- 1:nrow(xall)

    ## ##-----------------------------------------------------------------------------
    ## ## consider: N = observations
    ## ##           J = number of unique known (uncensored) observations
    ## N  <- nrow(xall)
    ## ## xu <- na.omit(unique(xall[xall$x.low == xall$x.high,]))
    ## ## J  <- nrow(xu)
    ## ## if (J < 1) {
    ## ##     cat('ERROR: No known values\n')
    ## ##     return()
    ## ## }
    ## 
    ## ## create dataframe of the unique x values regardless of whether start or end of an interval
    ## xall <- as_tibble(xall)
    ## nobs <- nrow(xall)
    ## x.all <- c(xall$x.low, xall$x.high)
    ## xu <- data.frame(x=unique(sort(x.all)))
    ## 
    ## xu$fail <- NA
    ## xu$surv <- NA
    ## xu$fail.ratio <- NA
    ## xu$surv.ratio <- NA
    ## for (i in 1:nrow(xu)) {
    ##   
    ##     ## count values that are less than or equal to x
    ##     xu$fail[i] <- sum(xall$x.high <= xu$x[i], na.rm=TRUE)
    ##     
    ##     ## count values that are greater than x
    ##     xu$surv[i] <- sum(xall$x.low > xu$x[i], na.rm=TRUE)
    ##     
    ##     ## observed ratios
    ##     xu$fail.ratio[i] <- xu$fail[i] / nobs
    ##     xu$surv.ratio[i] <- xu$surv[i] / nobs
    ##     
    ## }
    ## ## cumulative survival probability                         # useless?
    ## xu$fail.cumsurv <- cumprod(1 - xu$fail.ratio)
    ## xu$surv.cumsurv <- cumprod(1 - xu$surv.ratio)
    ## 
    ## ## johnson su fit
    ## out.fit <- mle.johnsonsu(xcen=xall)
    ## params <- out.fit$parms
    ## gamma <- params$gamma
    ## delta <- params$delta
    ## xi    <- params$xi
    ## lambda <- params$lambda
    ## 
    ## ## determine CDF from predicted fit at each xu value
    ## xu$cdf <- stats::pnorm(gamma+delta*asinh((xu$x-xi)/lambda))
    ## 
    ## ## plot cumulative distribution for data vs. fit
    ## ## plotspace(1,2)
    ## ## plot(xu$surv.cumsurv, xu$cdf)                               # useless?
    ## plot(xu$cdf, xu$fail.ratio, xlab='CDF of observed value', ylab='Failure ratio', main=mainadder)
    ## abline(0,1, col='red')
    ## 
    ## return(xu)
    ##
    ##-----------------------------------------------------------------------------

    type.list <- strsplit(type, "")[[1]]
    for (ichar in type.list) {
        
        if (grepl('n', ichar)) {
            cat('##################################################################\n')
            cat('functionality for censored, normal distribution not programmed yet\n')
            cat('##################################################################\n')

        } else if (grepl('w', ichar)) {
            if (length(type) > 1) {
                ## user provided parameters
                params <- wfit
            } else {
                ## determine fit parameters
                fit <- mle.weibull(x, xcen)
                params <- fit$parms
            }
            ## add theoretical quantiles
            xall$wtheoretical <- stats::qweibull(ppoints(nrow(xall)), shape=params$shape, scale=params$scale)
            main <- paste('Censored Weibull QQ Plot', mainadder, sep=" ")
            ## plot only known points
            xplot <- xall[xall$type == 'known',]
            plot(xplot$x.high, xplot$wtheoretical, xlab='Observed value, x', ylab='Expected value', main=main)
            abline(0,1,col='red')    
            
        } else if (grepl('j', ichar)) {
            if (length(type) > 1) {
                ## user provided parameters
                params <- jfit
            } else {
                ## determine fit parameters
                fit <- mle.johnsonsu(x, xcen)
                params <- fit$parms
            }
            ## add theoretical quantiles
            xall$jtheoretical <- SuppDists::qJohnson(ppoints(nrow(xall)), parms=params)
            main <- paste('Censored Johnson SU QQ Plot', mainadder, sep=" ")
                ## plot only known points
            xplot <- xall[xall$type == 'known',]
            plot(xplot$x.high, xplot$jtheoretical, xlab='Observed value, x', ylab='Expected value', main=main)
            abline(0,1,col='red')    
    
        }
        
    }
    
    return(list(xall=xall, xplot=xplot))
}

qqplot_nwj_test <- function() {
    source('setup.r')
  
    set.seed(1)
    x <- runif(1000, 1, 2)
    plotspace(1,4)
    fit <- mle.johnsonsu(x)
    out.hist <- hist_nwj(x, type='j', tolerance=FALSE)
    xcen.known <- data.frame(x.low=x, x.high=x)
    out.qq <- qqplot_nwj(x, type='j')
    out.qq <- qqplot_nwj(x=NA, xcen.known, type='j', mainadder='x as censored')
    xcen <- data.frame(x.low=c(NA, NA), x.high=c(1.3, 1.4))
    out.qq <- qqplot_nwj(x, xcen, type='j', mainadder = 'x + xcen')
    
    x <- iris$Sepal.Width    # range 2 to 4.4
    x <- c(x, 6)             # much larger than largest x
    xcen.known <- data.frame(x.low=x, x.high=x)
    xnum <- 1
    x.low  <- NA
    x.high <- 3.0
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
    plotspace(1,3)
    out.qq <- qqplot_nwj(x, type='j')
    out.qq <- qqplot_nwj(x=NA, xcen.known, type='j')   # not useful
    out.qq <- qqplot_nwj(x, xcen, type='j')
    
    
    ## subset of above
    x <- head(iris$Sepal.Width, 10)    # range 2 to 4.4
    x <- c(x, 6)             # much larger than largest x
    xcen.known <- data.frame(x.low=x, x.high=x)
    plotspace(1,3)
    out.qq <- qqplot_nwj(x, type='j')
    out.qq <- qqplot_nwj(x=NA, xcen.known, type='j', mainadder='x=NA')   # not useful
    out.qq <- qqplot_nwj(x, xcen, type='j', mainadder='x and xcen')
    
    
    set.seed(1)
    x <- rnorm(1000)
    x2 <- 4
    x <- c(x, x2)
    xcen <- data.frame(x.low=x, x.high=x)
    plotspace(1,3)
    out.hist <- hist_nwj(x, type='j', tolerance=FALSE)
    out.qq <- qqplot_nwj(x, xcen=NA, type='j')
    out.qq <- qqplot_nwj(x, xcen, type='j')

    
    ## two distributions thrown together
    x1 <- rnorm(200, 0.3, 0.1)
    x2 <- rnorm(200, 0.7, 0.1)
    x <- c(x1, x2)
    out.hist <- hist_nwj(x, type='j', tolerance=FALSE)
    out.qq <- qqplot_nwj(x, xcen=NA, type='j')
    xcen <- data.frame(x.low=x, x.high=x)
    out.qq <- qqplot_nwj(x=NA, xcen, type='j')
    
    
    
    ## example in paper
    x.low  <- c(NA, NA, 5, NA, 15)
    x.high <- c( 4,  4, 5, 14, 15)
    xcen <- data.frame(x.low=x.low, x.high=x.high)
    out.qq <- qqplot_nwj(x=NA, xcen=xcen, type='j')
    ##    x fail surv fail.ratio surv.ratio fail.cumsurv surv.cumsurv       cdf
    ## 1  4    2    2        0.4        0.4        0.600        0.600 0.5212980
    ## 2  5    3    1        0.6        0.2        0.240        0.480 0.5748700 <- known
    ## 3 14    4    1        0.8        0.2        0.048        0.384 0.9160881
    ## 4 15    5    0        1.0        0.0        0.000        0.384 0.9342868 <- known
    
    return(out.qq)
}
