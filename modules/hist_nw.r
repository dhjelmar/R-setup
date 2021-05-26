hist_nw <- function(x, alpha=0.01, pvalue=0.99, breaks=NULL, suppress='no') {
    ## plot histogram, normal distribution, and Weibull distribution
    ## add lines for upper tolerance limits for given alpha and pvalue
    ## alpha  = 1 - confidence
    ## pvalue = coverage (tolerance interval only)

    ## normal distribution calculations
    tol_out <- normtol.int(x, alpha = alpha, P=pvalue, side=1)
    upper_tolerance_limit_norm <- tol_out$'1-sided.upper'

    ## Weibull distribution calculations
    tol_out <-  exttol.int(x, alpha=alpha, P=pvalue, side=1, dist="Weibull")
    shape   <- tol_out$'shape.1'
    scale   <- tol_out$'shape.2'
    upper_tolerance_limit_weib <- tol_out$'1-sided.upper'

    ## create vectors with density distributions
    xmin    <- min(x)
    xmax    <- max(x,upper_tolerance_limit_norm,upper_tolerance_limit_weib)
    chuncks <- (xmax-xmin)/1000
    xrange  <- seq(xmin,xmax,by=chuncks)
    xmean   <- mean(x)
    xsd     <- sd(x)
    xdensity_norm <- dnorm(xrange,xmean,xsd)    
    xdensity_weib <- dweibull(xrange,shape=shape,scale=scale)
    maxdensity <- max(xdensity_norm, xdensity_weib)

    ## make histogram
    ## warning: xlim range can mess up x-axis
    ## obtain histogram parameters but suppress plot
    out <- hist(x, plot=FALSE)
    if (is.null(breaks)) breaks <- length(out$breaks)
    ymax <- max(out$density, maxdensity)
    ## create plot
    hist(x, breaks=breaks,
         xlim=c(xmin,xmax+(xmax-xmin)/breaks), 
         ylim=c(0,maxdensity),
         freq=FALSE)

    ## add distributions
    lines(x=xrange, y=xdensity_norm, col='red', lty=1)
    lines(x=xrange, y=xdensity_weib, col='blue', lty=1)

    ## add lines for mean and upper 1-sided 99/99 tolerance limits
    abline(v=xmean,col="red")
    abline(v=upper_tolerance_limit_norm,col="red",lty=2)
    abline(v=upper_tolerance_limit_weib,col="blue",lty=2)

    ## print to screen
    if (suppress == 'no') {
        cat("Tolerance limit input parameters:\n")
        cat("   upper, 1-sided,", pvalue*100,"/",(1-alpha)*100,"\n")
        cat("\n")
        cat("Normal distribution (red):\n")
        cat("   mean               =",xmean,"\n")
        cat("   standard deviation =",xsd,"\n")
        cat("   tolerance limit    =",upper_tolerance_limit_norm,"\n")
        cat("\n")
        cat("Weibull distribution (blue):\n")
        cat("   shape              =",shape,"\n")
        cat("   scale              =",scale,"\n")
        cat("   tolerance limit    =",upper_tolerance_limit_weib,"\n")
    }    
}
## hist_nw(mtcars$mpg)
