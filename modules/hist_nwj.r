hist_nwj <- function(x, type='nwj', nfit='standard', wfit='mle', jfit='mle', breaks=NULL,
                     tolerance=TRUE, side='upper', sided=1, P=0.99, conf=0.99,
                     main=NULL, subtitle='yes', suppress='no', plot=TRUE) {
    ## plot histogram and normal, Weibull, and Johnson distributions
    ## adds lines for upper tolerance limits for given alpha and proportion
    ## P = coverage proportion (tolerance interval only)
    ## conf = confidence used for tolerance limit
    ## breaks = number of bins used in the histogram (default auto determines number)
    ## nfit  = 'standard' uses mean(x) and sd(x) for parameters
    ##       = 'mle' uses maximum likelihood estimate
    ##       = vector or list for user supplied parameters (order: mean, sd; 1st 2 parameters used; names not used)
    ## wfit  = 'mle' uses maximum likelihood estimate
    ##       = 'exttol.int' uses tolerance::exttol.int to determine shape and scale
    ##       = vector or list for user supplied parameters (order: shape, scale; 1st 2 parameters used; names not used)
    ## jfit  = 'mle' uses maximum likelihood estimate
    ##       = 'SuppDists' uses SuppDists::JohnsonFit
    ##       = 'ExtDist'   uses ExtDist::eJohnsonSU 
    ##       = vector or list for user supplied parameters (order: gamma, delta, xi, lambda; 1st 4 parameters used; names not used)
    ##         e.g., jparms <- list(gamma   = -1.039,
    ##                              delta   = 1.66,
    ##                              xi      = 14.46,
    ##                              lambda  = 6.95,
    ##                              type    = 'SU')
    ##               mle.johnsonsu(mtcars$mpg, jfit=jparms)
    ## mle  = TRUE uses maximum likelihood estimate for fit and likelihood ratio for tolerance limit
    ## main = title for histogram (default is "Histogram of ...")
    ## subtitle = 'yes' = puts description of lines in subtitle
    ##          = 'no'  = subtitle is blank
    ##          = user specified, single line subtitle
    ## suppress = 'yes' creates plot but does not return calculated values to the screen
    ##            (e.g., fit parameters and tolerance limits)
    
    proportion = P
    
    ## set effective alpha level for use in standard R packages based on desired confidence and whether 1 or 2 sided
    alpha <- (1-conf)/sided

    ## name of variable passed in
    xname <- deparse(substitute(x))

    ## initialize parameters
    xmean <- mean(x)         # may get overwritten later
    xsd   <- sd(x)           # may get overwritten later
    shape <- NA
    scale <- NA
    jparms <- NA
    tol_out_norm <- NA
    tol_out_weib <- NA
    tol_out_john <- NA
    tolerance_limit_norm.l <- NA
    tolerance_limit_weib.l <- NA
    tolerance_limit_john.l <- NA
    tolerance_limit_norm.u <- NA
    tolerance_limit_weib.u <- NA
    tolerance_limit_john.u <- NA
    xdensity_norm <- NA
    xdensity_weib <- NA
    xdensity_john <- NA
    
    if (grepl('n', type)) {
        ## normal distribution calculations
	if (nfit == 'mle') {
            out.fit <- mle.normal(x, xcen)
            xmean <- out.fit$xbar
            xsd   <- out.fit$sdev
            tol_out_norm <- NA
            tolerance_limit_norm.l <- NA
            tolerance_limit_norm.u <- NA

        } else if (is.numeric(nfit[[1]])) {
            ## user supplied parameters
            xmean <- nfit[[1]]
            xsd   <- nfit[[2]]
            tol_out_norm <- NA
            tolerance_limit_norm.l <- NA
            tolerance_limit_norm.u <- NA

        } else {
            xmean <- mean(x)
            xsd   <- sd(x)
            tol_out_norm <- tolerance::normtol.int(x, alpha = alpha, P=proportion, side=sided)
            tolerance_limit_norm.l <- tol_out_norm[[4]]
            tolerance_limit_norm.u <- tol_out_norm[[5]]
        }
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
        } else {
            if (wfit == 'mle') {
                tol_out_weib <- mle.weibull.tol(x, side.which=side, sided=sided, conf=conf, P=P)
                shape <- tol_out_weib$params$shape
                scale <- tol_out_weib$params$scale
                tolerance_limit_weib.l <- tol_out_weib$tolerance$tol.lower
                tolerance_limit_weib.u <- tol_out_weib$tolerance$tol.upper
            } else if (wfit == 'exttol.int') {
                tol_out_weib <-  tolerance::exttol.int(x, alpha=alpha, P=proportion,
                                                       side=sided, dist="Weibull")
                shape   <- tol_out_weib$'shape.1'
                scale   <- tol_out_weib$'shape.2'
                tolerance_limit_weib.l <- tol_out_weib[[5]]
                tolerance_limit_weib.u <- tol_out_weib[[6]]
            } else if (is.numeric(wfit[[1]])) {
                ## user supplied parameters
                shape <- wfit[[1]]
                scale <- wfit[[2]]
                tol_out_weib <- NA
            } else {
                ## no other method programmed
                shape <- NA
                scale <- NA
                tol_out_weib <- NA
            }
        }
    }
        
    if (grepl('j', type)) {
        ## ## Johnson distribution calculations
        fit.j <- mle.johnsonsu(x)
        jparms <- fit.j$parms.compare
        if (jfit == 'mle' | isTRUE(tolerance)) {
            jparms   <- fit.j$parms
            if (isTRUE(tolerance)) {
                tol_out_john <- mle.johnsonsu.tol(x, param=jparms,
                                                  side.which=side, sided=sided, P=P, conf=conf, 
                                                  plots=FALSE, debug=FALSE)
                tolerance_limit_john.l <- tol_out_john$tolerance$tol.lower
                tolerance_limit_john.u <- tol_out_john$tolerance$tol.upper
                tol_out_john <- tol_out_john$tolerance
            }                
        } else if (jfit == 'SuppDists') {
            jparms <- jparms[jparms$description == 'SuppDists::JohnsonFit(x)', 1:5]
        } else if (jfit == 'ExtDist') {
            jparms <- jparms[jparms$description == 'ExtDist', 1:5]
        } else {
            jparms <- jfit
        }            
    }
    
    ## create vectors with density distributions
    xmin <- min(x)
    xmax <- max(x)
    if (isTRUE(tolerance)) {
        if (side != 'lower') {
            if (grepl('n', type)) xmax <- max(xmax, tolerance_limit_norm.u, na.rm=TRUE)
            if (grepl('w', type)) xmax <- max(xmax, tolerance_limit_weib.u, na.rm=TRUE)
            if (grepl('j', type)) xmax <- max(xmax, tolerance_limit_john.u, na.rm=TRUE)
        }
        if (side != 'upper') {
            if (grepl('n', type)) xmin <- min(xmin, tolerance_limit_norm.l, na.rm=TRUE)
            if (grepl('w', type)) xmin <- min(xmin, tolerance_limit_weib.l, na.rm=TRUE)
            if (grepl('j', type)) xmin <- min(xmin, tolerance_limit_john.l, na.rm=TRUE)
        }            
    }
    chuncks <- (xmax-xmin)/1000
    xrange  <- seq(xmin,xmax,by=chuncks)
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
    if (length(breaks) == 1) {
        hist(x, breaks=breaks,
             xlab = xname,
             xlim=c(xmin,xmax+(xmax-xmin)/breaks), 
             ylim=c(0,maxdensity),
             freq=FALSE,
             main=main,
             plot=plot)
    } else {
        hist(x, breaks=breaks,
             xlab = xname,
             xlim=c(xmin,xmax+(xmax-xmin)/length(breaks)), 
             ylim=c(0,maxdensity),
             freq=FALSE,
             main=main,
             plot=plot,
             xaxt='n')                           # do not plot and label the x-axis
        axis(side=1, at=breaks, labels=breaks)   # use breaks vector for x-axis instead
    }
    if (subtitle == 'yes') {
        n <- w <- j <- NULL
        if (grepl('n', type)) n <- 'red = normal,'
        if (grepl('w', type)) w <- 'blue = Weibull,'
        if (grepl('j', type)) j <- 'black = Johnson SU'
        ## subtitle <- list('line color: red = normal, blue = Weibull, black = Johnson SU',
        ##                  'line type: solid = distribution or mean, dashed = bound')
        subtitle <- list(paste('line color:', n, w, j, sep=' '),
                         'line type: solid = distribution or mean; dashed = bound')
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
        if (isTRUE(tolerance)) {
            if (side != 'lower') {
                if (grepl('n', type)) abline(v=tolerance_limit_norm.u,col="red",  lty=2)
                if (grepl('w', type)) abline(v=tolerance_limit_weib.u,col="blue", lty=2)
                if (grepl('j', type)) abline(v=tolerance_limit_john.u,col="black",lty=2)
            }
            if (side != 'upper') {
                if (grepl('n', type)) abline(v=tolerance_limit_norm.l,col="red",  lty=2)
                if (grepl('w', type)) abline(v=tolerance_limit_weib.l,col="blue", lty=2)
                if (grepl('j', type)) abline(v=tolerance_limit_john.l,col="black",lty=2)
            }            
        }

    }

    ## print to screen
    if (suppress == 'no') {
        cat(side, ',', sided, "-Sided Tolerance Limits\n")
        cat("coverage proportion, P = ", P, " / confidence = ", conf, '\n')
        cat("------------------------------------------------")
        cat("\n")
        if (grepl('n', type)) {
            cat("Normal distribution (red):\n")
            cat("   mean                =",xmean,"\n")
            cat("   standard deviation  =",xsd,"\n")
            cat("   tolerance limit     =",tolerance_limit_norm.u,"\n")
            cat("\n")
        }
        if (grepl('w', type)) {
            cat("Weibull distribution (blue):\n")
            cat("   shape               =",shape,"\n")
            cat("   scale               =",scale,"\n")
            cat("   tolerance limit     =",tolerance_limit_weib.u,"\n")
            cat("\n")
        }
        if (grepl('j', type)) {
            cat("Johnson distribution (black):\n")
            cat("   gamma               =",jparms$gamma ,"\n")
            cat("   delta               =",jparms$delta ,"\n")
            cat("   xi                  =",jparms$xi    ,"\n")
            cat("   lambda              =",jparms$lambda,"\n")
            cat("   type                =",jparms$type  ,"\n")
            cat("   tolerance limit     =",tolerance_limit_john.u,"\n")
        }
    }    

    return(list(side = side,
                sided = sided,
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
## x <- rnorm(1000, 100, 1)
## out <- hist_nwj(x)
## out <- hist_nwj(x, type='n')
## out <- hist_nwj(x, type='w')
## out <- hist_nwj(x, type='j')
## out <- hist_nwj(x, type='nw')
