hist_nwj <- function(x, alpha=0.01, pvalue=0.99, breaks=NULL, jfit='all', suppress='no') {
    ## plot histogram and normal, Weibull, and Johnson distributions
    ## add lines for upper tolerance limits for given alpha and pvalue
    ## alpha  = 1 - confidence
    ## pvalue = coverage (tolerance interval only)

    ## normal distribution calculations
    tol_out_norm <- tolerance::normtol.int(x, alpha = alpha, P=pvalue, side=1)
    upper_tolerance_limit_norm <- tol_out_norm$'1-sided.upper'

    ## Weibull distribution calculations
    if (min(x) < 0) {
        cat('\n')
        cat('#####################################################################\n')
        cat('Skip Weibull: Negative values cannot be fit with Weibull distribution\n')
        cat('              Some day I should add logic to skip if that is the case\n')
        cat('#####################################################################\n')
    } else {
        tol_out_weib <-  tolerance::exttol.int(x, alpha=alpha, P=pvalue, side=1, dist="Weibull")
        shape   <- tol_out_weib$'shape.1'
        scale   <- tol_out_weib$'shape.2'
        upper_tolerance_limit_weib <- tol_out_weib$'1-sided.upper'
    }
        
    ## Johnson distribution calculations
    tol_out_john <- johnson_tol(x, alpha=alpha, P=pvalue, side=1, jfit=jfit)
    jparms   <- tol_out_john$jparms
    if (jparms$type == 'SB') {
        upper_tolerance_limit_john <- NA
    } else {
        upper_tolerance_limit_john <- tol_out_john$bounds$xtol[2]
    } 
       
    ## create vectors with density distributions
    xmin    <- min(x)
    xmax    <- max(x,
                   upper_tolerance_limit_norm,
                   upper_tolerance_limit_weib,
                   upper_tolerance_limit_john,
                   na.rm = TRUE)
    chuncks <- (xmax-xmin)/1000
    xrange  <- seq(xmin,xmax,by=chuncks)
    xmean   <- mean(x)
    xsd     <- sd(x)
    xdensity_norm <- stats::dnorm(xrange,xmean,xsd)    
    xdensity_weib <- stats::dweibull(xrange,shape=shape,scale=scale)
    xdensity_john <- SuppDists::dJohnson(xrange, jparms)

    ## make histogram
    ## warning: xlim range can mess up x-axis
    ## obtain histogram parameters but suppress plot
    if (is.null(breaks)) {
        ## breaks not specified so figure out what to use
        out <- hist(x, plot=FALSE)
        breaks <- length(out$breaks)
    }
    ## hist needed with number of breaks to be used to get density for max height of histogram
    out <- hist(x, breaks=breaks, plot=FALSE)
    maxdensity <- max(xdensity_norm, xdensity_weib, xdensity_john, out$density, na.rm=TRUE)
    ymax <- max(out$density, maxdensity)
    ## create plot
    hist(x, breaks=breaks,
         xlim=c(xmin,xmax+(xmax-xmin)/breaks), 
         ylim=c(0,maxdensity),
         freq=FALSE)

    ## add distributions
    lines(x=xrange, y=xdensity_norm, col='red',   lty=1)
    lines(x=xrange, y=xdensity_weib, col='blue',  lty=1)
    lines(x=xrange, y=xdensity_john, col='black', lty=1)

    ## add lines for mean and upper 1-sided 99/99 tolerance limits
    abline(v=xmean,col="red")
    abline(v=upper_tolerance_limit_norm,col="red",  lty=2)
    abline(v=upper_tolerance_limit_weib,col="blue", lty=2)
    abline(v=upper_tolerance_limit_john,col="black",lty=2)

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
        cat("\n")
        cat("Johnson distribution (black):\n")
        cat("   gamma              =",jparms$gamma ,"\n")
        cat("   delta              =",jparms$delta ,"\n")
        cat("   xi                 =",jparms$xi    ,"\n")
        cat("   lambda             =",jparms$lambda,"\n")
        cat("   type               =",jparms$type  ,"\n")
        cat("   tolerance limit    =",upper_tolerance_limit_john,"\n")
    }    

    return(list(tol_out_norm=tol_out_norm,
                tol_out_weib=tol_out_weib,
                tol_out_john=tol_out_john))

}
## out   <- hist_nwj(mtcars$mpg)
## gamma <- out$tol_out_john$jparms$gamma
