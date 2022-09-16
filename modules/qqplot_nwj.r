qqplot_nwj <- function(x=NA, xcen=NA, type='nwj', nfit='mle', wfit='mle', jfit='mle', mainadder=NULL) {
    ## creates normal, Weibull and/or Johnson SU qq plots
    
    ## input: x    = vector of data
    ##        xcen = dataframe of censored data

    if (!is.data.frame(xcen)) {
        ## no censored data were supplied, so:
	##    - use standard qqplot_nwj function
	##    - no reason to use mle fit for normal (and not currently supported by qqplot_nwj_xonly)
        out <- qqplot_nwj_xonly(x, type=type, wfit=wfit, jfit=jfit, mainadder=mainadder)
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

    ## create qqplot of only x (known) data
    out <- qqplot_nwj_xonly(x, type=type, wfit=wfit, jfit=jfit,
                            mainadder=paste(mainadder, '; censored data, if any, ignored', sep=''))
    return(list(out=out, x=x, xcen=xcen, xall=xall))


    ##-----------------------------------------------------------------------------
    ## Kaplan and Meier method for qqplot of censored data
    
    ## I found something here, but not clear it is what I am looking for and did not seem to work.
    ##     https://pdixon.stat.iastate.edu/stat505/Chapter%2011.pdf, p. 7
    ## The following is that attempt.
    
    ## ## sort data from low to high and renumber
    ## ## for the same x.high value, sort so censored observations are below known observations
    ## xall <- xall[order(xall$x.low, na.last=FALSE),]
    ## xall <- xall[order(xall$x.high),]
    ## rownames(xall) <- 1:nrow(xall)
    ## 
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
    ## ## ## johnson su fit
    ## ## out.fit <- mle.johnsonsu(xcen=xall)
    ## ## params <- out.fit$parms
    ## ## gamma <- params$gamma
    ## ## delta <- params$delta
    ## ## xi    <- params$xi
    ## ## lambda <- params$lambda
    ## ## ## determine CDF from predicted fit at each xu value
    ## ## xu$cdf <- stats::pnorm(gamma+delta*asinh((xu$x-xi)/lambda))
    ##
    ## ## weibull fit
    ## out.fit <- mle.weibull(xcen=xall)
    ## params <- out.fit$parms
    ## shape <- params$shape
    ## scale <- params$scale
    ## ## determine CDF from predicted fit at each xu value
    ## ## xu$cdf <- stats::pnorm(1 - exp(-(xu$x/scale)^shape))
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


    ##-----------------------------------------------------------------------------
    ## Following was another failed attempt using probability of known values for y-axis.
    ## Problem is probably knowing what to use for the x-axis for each point.

    ## type.list <- strsplit(type, "")[[1]]
    ## for (ichar in type.list) {
    ##     
    ##     if (grepl('n', ichar)) {
    ##         if (length(type) > 1) {
    ##             ## user provided parameters
    ##             params <- nfit
    ##         } else {
    ##             ## determine fit parameters
    ##             fit <- mle.normal(x, xcen)
    ##             params <- fit$parms
    ##         }
    ##         ## add theoretical quantiles
    ##         xall$wtheoretical <- stats::qnorm(ppoints(nrow(xall)), mean=params[[1]], sd=params[[2]])
    ##         main <- paste('Censored Weibull QQ Plot (not checked)', mainadder, sep=" ")
    ##         ## plot only known points
    ##         xplot <- xall[xall$type == 'known',]
    ##         plot(xplot$x.high, xplot$wtheoretical, xlab='Observed value, x', ylab='Expected value', main=main)
    ##         abline(0,1,col='red')    
    ## 
    ##     } else if (grepl('w', ichar)) {
    ##         if (length(type) > 1) {
    ##             ## user provided parameters
    ##             params <- wfit
    ##         } else {
    ##             ## determine fit parameters
    ##             fit <- mle.weibull(x, xcen)
    ##             params <- fit$parms
    ##         }
    ##         ## add theoretical quantiles
    ##         xall$wtheoretical <- stats::qweibull(ppoints(nrow(xall)), shape=params$shape, scale=params$scale)
    ##         main <- paste('Censored Weibull QQ Plot (not checked)', mainadder, sep=" ")
    ##         ## plot only known points
    ##         xplot <- xall[xall$type == 'known',]
    ##         plot(xplot$x.high, xplot$wtheoretical, xlab='Observed value, x', ylab='Expected value', main=main)
    ##         abline(0,1,col='red')    
    ##         
    ##     } else if (grepl('j', ichar)) {
    ##         if (length(type) > 1) {
    ##             ## user provided parameters
    ##             params <- jfit
    ##         } else {
    ##             ## determine fit parameters
    ##             fit <- mle.johnsonsu(x, xcen)
    ##             params <- fit$parms
    ##         }
    ##         ## add theoretical quantiles
    ##         xall$jtheoretical <- SuppDists::qJohnson(ppoints(nrow(xall)), parms=params)
    ##         main <- paste('Censored Johnson SU QQ Plot (not checked)', mainadder, sep=" ")
    ##             ## plot only known points
    ##         xplot <- xall[xall$type == 'known',]
    ##         plot(xplot$x.high, xplot$jtheoretical, xlab='Observed value, x', ylab='Expected value', main=main)
    ##         abline(0,1,col='red')    
    ## 
    ##     }
    ##     
    ## }
    ## 
    ## return(list(xall=xall, xplot=xplot))
}

qqplot_nwj_test <- function() {
    set.seed(1)
    x <- rnorm(1000, 1, 2)
    out <- qqplot_nwj(x)
}
