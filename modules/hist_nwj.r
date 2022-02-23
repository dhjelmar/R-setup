hist_nwj <- function(x, type='nwj', alpha=0.01, P=0.99, breaks=NULL, jfit='all',
                     upperbound=TRUE, main=NULL, subtitle='yes', suppress='no', plot=TRUE) {
    ## plot histogram and normal, Weibull, and Johnson distributions
    ## adds lines for upper tolerance limits for given alpha and proportion
    ## alpha  = 1 - confidence
    ## P = coverage proportion (tolerance interval only)
    proportion = P
    ## breaks = number of bins used in the histogram (default auto determines number)
    ## jfit  = 'all' = uses SuppDists::JohnsonFit to determine parameters
    ##       = 'SU'  = uses ExtDist::eJohnsonSU to determine parameters
    ##       = list of user specified parameters
    ##         e.g., jparms <- list(gamma   = -1.039,
    ##                              delta   = 1.66,
    ##                              xi      = 14.46,
    ##                              lambda  = 6.95,
    ##                              type    = 'SU')
    ##              johnson_tol(mtcars$mpg, jfit=jparms)
    ## main = title for histogram (default is "Histogram of ...")
    ## subtitle = 'yes' = puts description of lines in subtitle
    ##          = 'no'  = subtitle is blank
    ##          = user specified, single line subtitle
    ## suppress = 'yes' creates plot but does not return calculated values to the screen
    ##            (e.g., fit parameters and tolerance limits)
    
    ## name of variable passed in
    xname <- deparse(substitute(x))

    ## initialize parameters
    xmean <- NA
    xsd   <- NA
    shape <- NA
    scale <- NA
    jparms <- NA
    tol_out_norm <- NA
    tol_out_weib <- NA
    tol_out_john <- NA
    upper_tolerance_limit_norm <- NA
    upper_tolerance_limit_weib <- NA
    upper_tolerance_limit_john <- NA
    xdensity_norm <- NA
    xdensity_weib <- NA
    xdensity_john <- NA
    
    if (grepl('n', type)) {
        ## normal distribution calculations
        tol_out_norm <- tolerance::normtol.int(x, alpha = alpha, P=proportion, side=1)
        upper_tolerance_limit_norm <- tol_out_norm$'1-sided.upper'
    }
        
    if (grepl('w', type)) {
        ## Weibull distribution calculations
        if (min(x) < 0) {
            cat('\n')
            cat('#####################################################################\n')
            cat('Skip Weibull: Negative values cannot be fit with Weibull distribution\n')
            cat('#####################################################################\n\n')
            tol_out_weib <- NA
            shape <- NA
            scale <- NA
            upper_tolerance_limit_weib <- NA
        } else {
            tol_out_weib <-  tolerance::exttol.int(x, alpha=alpha, P=proportion, side=1, dist="Weibull")
            shape   <- tol_out_weib$'shape.1'
            scale   <- tol_out_weib$'shape.2'
            upper_tolerance_limit_weib <- tol_out_weib$'1-sided.upper'
        }
    }
        
    if (grepl('j', type)) {
        ## ## Johnson distribution calculations
        ## tol_out_john <- johnson_tol(x, alpha=alpha, P=proportion, side=1, jfit=jfit)
        ## jparms   <- tol_out_john$jparms
        ## upper_tolerance_limit_john <- tol_out_john$xtol_upper
        ## Johnson distribution calculations
        tol_out_john <- mle.johnsonsu(x, jfit, alpha=alpha, P=P, sided=1, plots=FALSE, debug=FALSE)
        jparms   <- tol_out_john$params
        upper_tolerance_limit_john <- tol_out_john$tolerance$tol.upper
    }
    
    ## create vectors with density distributions
    xmin    <- min(x)
    if (isTRUE(upperbound)) {
        xmax <- max(x,
                    upper_tolerance_limit_norm,
                    upper_tolerance_limit_weib,
                    upper_tolerance_limit_john,
                    na.rm = TRUE)
    } else {
        xmax <- max(x, na.rm = TRUE)
    }
    chuncks <- (xmax-xmin)/1000
    xrange  <- seq(xmin,xmax,by=chuncks)
    xmean   <- mean(x)
    xsd     <- sd(x)
    if (grepl('n', type)) xdensity_norm <- stats::dnorm(xrange,xmean,xsd)    
    if (grepl('w', type)) xdensity_weib <- stats::dweibull(xrange,shape=shape,scale=scale)
    if (grepl('j', type)) xdensity_john <- SuppDists::dJohnson(xrange, jparms)

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
    if (is.null(main)) main <- paste('Histogram of', xname, sep=" ")
    hist(x, breaks=breaks,
         xlab = xname,
         xlim=c(xmin,xmax+(xmax-xmin)/breaks), 
         ylim=c(0,maxdensity),
         freq=FALSE,
         main=main,
         plot=plot)
    if (subtitle == 'yes') {
        n <- w <- j <- NULL
        if (grepl('n', type)) n <- 'red = normal,'
        if (grepl('w', type)) w <- 'blue = Weibull,'
        if (grepl('j', type)) j <- 'black = Johnson'
        ## subtitle <- list('line color: red = normal, blue = Weibull, black = Johnson',
        ##                  'line type: solid = distribution or mean, dashed = 1-sided upper bound')
        subtitle <- list(paste('line color:', n, w, j, sep=' '),
                         'line type: solid = distribution or mean, dashed = 1-sided upper bound')
        mtext(subtitle, side=3, line=c(0.75, 0), cex=.75, col='black')
    } else if (subtitle == 'no') {
        ## skip subitile
    } else {
        ## use user specified subitile
        mtext(subtitle)
    }


    if (isTRUE(plot)) {
        
        ## add distributions
        if (grepl('n', type)) lines(x=xrange, y=xdensity_norm, col='red',   lty=1)
        if (grepl('w', type)) lines(x=xrange, y=xdensity_weib, col='blue',  lty=1)
        if (grepl('j', type)) lines(x=xrange, y=xdensity_john, col='black', lty=1)
        
        ## add lines for mean and upper 1-sided 99/99 tolerance limits
        abline(v=xmean,col="red")
        if (isTRUE(upperbound)) {
            if (grepl('n', type)) abline(v=upper_tolerance_limit_norm,col="red",  lty=2)
            if (grepl('w', type)) abline(v=upper_tolerance_limit_weib,col="blue", lty=2)
            if (grepl('j', type)) abline(v=upper_tolerance_limit_john,col="black",lty=2)
        }

    }

    ## print to screen
    if (suppress == 'no') {
        cat("Upper, 1-Sided Tolerance Limits\n")
        cat("confidence =", 1-alpha, ", coverage proportion", proportion, '\n')
        cat("-------------------------------")
        cat("\n")
        if (grepl('n', type)) {
            cat("Normal distribution (red):\n")
            cat("   mean                =",xmean,"\n")
            cat("   standard deviation  =",xsd,"\n")
            cat("   tolerance limit     =",upper_tolerance_limit_norm,"\n")
            cat("\n")
        }
        if (grepl('w', type)) {
            cat("Weibull distribution (blue):\n")
            cat("   shape               =",shape,"\n")
            cat("   scale               =",scale,"\n")
            cat("   tolerance limit     =",upper_tolerance_limit_weib,"\n")
            cat("\n")
        }
        if (grepl('j', type)) {
            cat("Johnson distribution (black):\n")
            cat("   gamma               =",jparms$gamma ,"\n")
            cat("   delta               =",jparms$delta ,"\n")
            cat("   xi                  =",jparms$xi    ,"\n")
            cat("   lambda              =",jparms$lambda,"\n")
            cat("   type                =",jparms$type  ,"\n")
            cat("   tolerance limit     =",upper_tolerance_limit_john,"\n")
            cat("   (tolerance limit is the distribution bound if type SB or SL)\n")
        }
    }    

    return(list(sided = '1-sided, upper limits',
                proportion = proportion,
                alpha = alpha,
                nparms = list(mean = xmean, sd=xsd),
                wparms = list(shape = shape, scale = scale),
                jparms = jparms,
                tol_out_norm=tol_out_norm,
                tol_out_weib=tol_out_weib,
                tol_out_john=tol_out_john,
                xrange = xrange,
                xdensity_norm = xdensity_norm,
                xdensity_weib = xdensity_weib,
                xdensity_john = xdensity_john))

}

## set.seed(1)
## jparms <- list(gamma = -1.039, delta = 1.66, xi = 14.46, lambda = 6.95, type    = 'SU')
## xjohn <- SuppDists::rJohnson(999, parms=jparms) + 2
## out   <- hist_nwj(xjohn, jfit=jparms)
## 
## out <- hist_nwj(mtcars$mpg)
## out <- hist_nwj(mtcars$mpg, type='n')
## out <- hist_nwj(mtcars$mpg, type='w')
## out <- hist_nwj(mtcars$mpg, type='j')
## out <- hist_nwj(mtcars$mpg, type='nw')
