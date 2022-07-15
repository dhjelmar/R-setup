loglik.normal <- function(x=NA, xcen=NA, param=c(xbar, sdev), debug=FALSE){
    ## calculate log likelihhod for distribution
    if (is.data.frame(x)) x <- x[[1]] # convert to vector
    if (is.data.frame(xcen)) {
        ## if known data are inside xcen, move them to x and keep remainder in xcen
        xcen.na  <- xcen[ is.na(rowSums(xcen)),]                  # censored rows with NA, if any
        xcen.val <- xcen[!is.na(rowSums(xcen)),]                  # rows w/o NA, if any
        x.add   <- xcen.val[xcen.val[[1]] == xcen.val[[2]],][[1]] # known values from xcen, if any
        x <- as.numeric(na.omit(c(x, x.add)))                     # new set of known values
        xcen.lowhigh <- xcen.val[xcen.val[[1]] != xcen.val[[2]],] # censored rows with no NA
        xcen <- rbind(xcen.lowhigh, xcen.na)                      # new set of censored rows
        names(xcen) <- c('x.low', 'x.high')                       # rename
    }
    ## set parameters
    xbar  <- param[[1]]  
    sdev  <- param[[2]]
    ## likelihood contribution from x
    pdf <- 1
    if (!is.na(x[1])) {
        z   <- (x - xbar) / sdev
        pdf <- 1/(sqrt(2 * pi) * sdev) * exp(-z^2 / 2)
        ## the following is equivalent to the above
        ## pdf <- stats::dnorm(x, mean=xbar, sd=sdev, log=FALSE)
    }
    ## likelihood contribution from xcen
    probability <- 1
    if (is.data.frame(xcen)) {
        xcen$F.low  <- stats::pnorm(xcen$x.low , mean=xbar, sd=sdev)
        xcen$F.high <- stats::pnorm(xcen$x.high, mean=xbar, sd=sdev)
        ## if low CDF is NA, set to 0
        xcen$F.low[is.na(xcen$F.low)]   <- 0
        ## if high CDF is NA, set to 1
        xcen$F.high[is.na(xcen$F.high)] <- 1
        ## calculate probability for the censored interval
        probability <- xcen$F.high - xcen$F.low
    }
    loglik <- sum(log(pdf), log(probability))
    if (isTRUE(debug)) cat('mean  =', signif(xbar,11),
                           'sd    =', signif(sdev,11),
                           'loglik=', signif(loglik,11), "\n")
    return(loglik)
}        

mle.normal <- function(x=NA, xcen=NA, param='auto', plots=FALSE, plot3d=FALSE, debug=FALSE) {
    
    ## normal distribution
    ## MLE (Maximum Likelihood Estimate) fit to determine parameters
    ## LR (Likelihood Ratio) appraoch to find tolerance limit

    ## input: x     = vector of known data
    ##        xcen  = dataframe of censored data
    ##                (1st column = low value or NA; 2nd column = high value or NA)
    ##        param = initial guess for fit parameters for: mean and sd
    ##                if type is also provided, it will not be used
    ##              = 'auto' (default) uses calculated mean and standard deviation for initial guess of parameters
    
    ## options: plots  = create histogram with fit and qqplot (qqplot currently based only on x)
    ##          plot3d = create 3D plot of the fit; max likelihood shown with red sphere
    
    ## based on approach found here:
    ## https://www.r-bloggers.com/2019/08/maximum-likelihood-estimation-from-scratch/
    ## https://personal.psu.edu/abs12/stat504/Lecture/lec3_4up.pdf

    ## alternately, could probably use bbmle::mle2() to determine coefficients
    ## an example is here: https://stackoverflow.com/questions/45208176/the-weibull-distribution-in-r-extdist

    if (is.data.frame(x)) x <- x[[1]] # convert to vector
  
    if (is.data.frame(xcen)) {
        ## censored data also provided (only reason for following is if
        ## x.low and x.high were not the names of the two columns of data)
        xcen <- data.frame(x.low = xcen[[1]], x.high = xcen[[2]])

        ## calculate average value, ignoring NA, xcen to for use in estimating parameters
        ## from packages that do not have censor capability
        xcen.avg <- rowMeans(xcen, na.rm=TRUE)

    } else {
        xcen.avg <- NA
    }
    x.avg <- as.numeric( na.omit( c(x, xcen.avg) ) )
    
    ## if (isTRUE(plots)) par(mfrow=c(1,2))
    out <- NULL

    ##-----------------------------------------------------------------------------
    ## set initial guess for MLE fit
    if (param[1] == 'auto') {
        ## initial parameters not specified, so let R figure guess
        tol_out_norm <- NA
        xbar   <- mean(x)
        sdev   <- sd(x)
        param <- list(xbar=xbar, sdev=sdev)
        params.compare <- as.data.frame(param)
        params.compare$loglik <- NA
        params.compare$description <- 'calculated mean and sd for known points'
    } else {
        ## use provided initial parameters
        xbar <- param[[1]]
        sdev <- param[[2]]
        params.compare <- as.data.frame(param)
        params.compare$loglik <- NA
        params.compare$description <- 'user provided parameters'
    }
    
    temp <- as.data.frame(t(unlist(param)))
    temp$loglik <- loglik.normal(x, xcen, param=c(xbar, sdev), debug=FALSE)
    temp$description <- 'initial guess for MLE'
    params.compare <- fastmerge(params.compare, temp)
    
    ##-----------------------------------------------------------------------------
    ## determine best fit using nll
    ## print('Attempting MLE fit on regular parameters')
    out <- NA
    out <- maxLik::maxLik(loglik.normal,
                          start = unlist(param),
                          x     = x,
                          xcen  = xcen,
                          debug = debug,
                          iterlim = 2000)
    parms.mle <- as.list(out$estimate)
    xbar <- parms.mle$xbar
    sdev <- parms.mle$sdev
    loglik <- out$maximum
    convergence <- if (out$code == 0) {'successful'}
                   else               {out$message}
    if (convergence != 'successful') {
        cat('###############################################\n')
        cat('WARNING: CONVERGENCE FAILURE IN mle.normal()   \n')
        cat('###############################################\n')
    }
    
    ## add MLE parameters to params.compare dataframe
    temp <- as.data.frame(parms.mle)
    temp$loglik <- loglik.normal(x, xcen, param=c(xbar, sdev), debug=FALSE)
    temp$description <- 'MLE'
    params.compare <- fastmerge(params.compare, temp)

    ## remove factor levels from params.compare
    params.compare <- droplevels.all(params.compare)

    if (isTRUE(plots)) {
        out.hist <- hist(x, plot=FALSE)
        curve.points <- stats::dnorm(x, mean=xbar, sd=sdev)
        hist(x, freq=FALSE, ylim=range(out.hist$density, curve.points))
        curve(stats::dnorm(x, mean=xbar, sd=sdev), min(x), max(x), add=TRUE)
        ## xcen likely not correctly used in qqplot
        qqplot_nwj(x, xcen, type='n', jfit=parms.mle)
    }
    
    if (isTRUE(plot3d)) {
        ## plot the fit
        loglik.plot <- df.init(c('param1', 'param2', 'loglik'))
        param1 <- seq(parms.mle[[1]]/1.05, parms.mle[[1]]*1.05, length.out=100)
        param2 <- seq(parms.mle[[2]]/1.05, parms.mle[[2]]*1.05, length.out=100)
        j <- 0
        for (i1 in 1:length(param1)) {
            for (i2 in 1:length(param2)) {
                j <- j+1
                loglik.plot[j, 1] <- param1[i1]
                loglik.plot[j, 2] <- param1[i2]
                loglik.plot[j, 3] <- loglik.normal(x, xcen, param=c(loglik.plot$param1[j], loglik.plot$param2[j]))
            }
        }
        rgl::plot3d(loglik.plot, xlab='mean', ylab='sd', zlab='loglik')
        rgl::points3d(x = parms.mle[[1]], y = parms.mle[[2]], z=loglik,
                      col='red', size=10, add=TRUE)
    }

    return(list(parms=parms.mle, parms.compare=params.compare, loglik=loglik, convergence=convergence))
}



mle.normal.test <- function() {
    source("F:\\Documents\\01_Dave's Stuff\\Programs\\GitHub_home\\R-setup\\setup.r")
    source('setup.r')
    
    fit.compare <- function(x, xcen, main=NULL) {
        ## plot histogram with no censored data and fit
        hist(x, freq=FALSE, border='black', main=main)
        out.fit0 <- mle.normal(x, xcen=NA, plots=FALSE)
        parms0 <- out.fit0$parms
        curve(stats::dnorm(x, parms0$xbar, parms0$sdev), min(x), max(x), col='black', add=TRUE)
        ## plot histogram with all data treated as known fit at numeric value
        xcen.avg <- rowMeans(xcen, na.rm=TRUE) # use the average for interval data
        x.all <- c(x, xcen.avg)
        hist(x.all, freq=FALSE, border='red', add=TRUE)
        out.fit1 <- mle.normal(x.all, xcen=NA, plots=FALSE)
        parms1 <- out.fit1$parms
        curve(stats::dnorm(x, parms1$xbar, parms1$sdev), min(x), max(x), col='red', add=TRUE)
        ## plot fit if treat xcen as censored
        out.fit2 <- mle.normal(x, xcen=xcen, plots=FALSE)
        parms2 <- out.fit2$parms
        curve(stats::dnorm(x, parms2$xbar, parms2$sdev), min(x), max(x), col='blue', type='p', add=TRUE)
        ## add legend
        legend('topright', 
               legend=c('known', 'all',   'censored'),
               col   =c('black', 'red', 'blue'),
               lty   =c( 1     ,  1   ,  NA),
               pch   =c( NA    ,  NA  ,  1))
    }

    plotspace(2,2)
    x <- iris$Sepal.Width
    xnum <- 10
    x.low  <- NA
    x.high <- 2.2
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit.compare(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    
    x.low  <- 3.7
    x.high <- NA
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit.compare(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    
    x.low  <- 2.2
    x.high <- 3.7
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit.compare(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))
    
    x.low  <- 2.2
    x.high <- NA
    xcen <- data.frame(x.low=rep(x.low,xnum), x.high=rep(x.high,xnum))
    fit.compare(x, xcen, main=paste(xnum, 'censored points from', x.low, 'to', x.high, sep=' ' ))

    ## plot the fit
    out.fit <- mle.normal(x, xcen, plot3d=TRUE)

}

