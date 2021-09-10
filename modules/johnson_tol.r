johnson_tol <- function(xall, alpha=0.01, P=0.99, side=1, jfit='all', plots='no') {

    ## xall = data set
    ## alpha  = 1 - confidence
    ## P = proportion = coverage (tolerance interval only)
    proportion <- P

    ##--------------------------------------------------------
    ## initialize tolerance limits
    ztol_lower   = NA
    ztol_upper   = NA
    xntol_lower  = NA
    xntol_upper  = NA
    xtol_lower   = NA
    xtol_upper   = NA

    
    ##--------------------------------------------------------
    ## determine Johnson parameters
    
    if (jfit == 'all') {
        ## let R figure out which Johnson distribution fits best
        jparms  <- SuppDists::JohnsonFit(xall)
    } else if (jfit == 'SU') {
        ## force the Johnson SU distribution
        jparms.out <- ExtDist::eJohnsonSU(xall)
        jparms <- list(gamma   = jparms.out$gamma,
                       delta   = jparms.out$delta,
                       xi      = jparms.out$xi,
                       lambda  = jparms.out$lambda,
                       type    = 'SU')
    } else {
        ## use Johnson parameters specified in jfit
        ## needs to be in same list format as created by SuppDists::JohnsonFit
        jparms <- jfit
    }
    gamma  <- jparms$gamma
    delta  <- jparms$delta
    xi     <- jparms$xi
    lambda <- jparms$lambda
    type   <- jparms$type

    
    ##--------------------------------------------------------
    ## transform xall to normal distribution, xnall
    ## https://www.sigmamagic.com/blogs/how-do-i-transform-data-to-normal-distribution/

    if (type == 'SU') {
        ## xn[1] <- gamma + delta * asinh((x[1] - xi)/lambda)
        ## xnall[1] <- gamma + delta * asinh( (xall[1] - xi)/ lambda )
        ## xn[1]
        ## xnall[1]
        xnall <- gamma + delta * asinh( (xall - xi)/ lambda )
    } else if (type == 'SB') {
        xnall <- try( gamma + delta *   log( (xall - xi)/(lambda + xi - xall) ) )
    } else {
        ## type == 'SL'
        xnall <- try( gamma + delta *   log( (xall - xi)/ lambda ) )  
    }
        
    ##--------------------------------------------------------
    ## transform to standard normal, zall

    zall <- (xnall - mean(xnall)) / sd(xnall)
    df <- as_tibble( data.frame(x=xall, z=zall) )

    if (plots != 'no') {
        plotspace(2,2)
        hist(xall)
        hist(xnall)
        hist(zall)
    }

    
    ##--------------------------------------------------------
    ## error check then proceed
    na_rows <- df[is.na(df$z),]
    if (nrow(na_rows) != 0) {

        ## z values could not be calculated for all x values
        cat('\n')
        cat('####################################\n')
        cat('  ERROR IN JOHNSON_TOL  \n')
        cat('  z(x) = NA or NaN for some x values \n')
        cat('####################################\n')
        print(jparms)
        print(na_rows)
        ztol_out <- normtol.int(zall, alpha = alpha, P=proportion, side=side)
        
    } else {
        ## johnson fit returned z values for all x values so possibly a decent fit
        
        ##--------------------------------------------------------
        ## calculate standard normal tolerance limit, ztol (i.e., for standard normal)
    
        ztol_out <- normtol.int(zall, alpha = alpha, P=proportion, side=side)
        xmin <- min(xall)
        xmax <- max(xall)
        if (side == 1) {
            ztol_lower <- ztol_out$'1-sided.lower'
            ztol_upper <- ztol_out$'1-sided.upper'
        } else if (side == 2) {
            ztol_lower <- ztol_out$'2-sided.lower'
            ztol_upper <- ztol_out$'2-sided.upper'
        }


        ##-----------------------------------------------------------------------------        
        ## transform ztol back to normal tolerance limit, xntol
        xntol_lower <- mean(xnall) + sd(xnall) * ztol_lower # lower z is negative
        xntol_upper <- mean(xnall) + sd(xnall) * ztol_upper

        
        ##-----------------------------------------------------------------------------        
        ## transform xntol back to johnson tolerance limit, xtol
        
        if (type == 'SB') {
            ## Not sure tolerance limits have meaning for a bounded fit
            xtol_lower <- NA
            xtol_upper <- NA
            
        } else {
            ## Tolerance limits should have meaning for SU and SN
            ## Upper tolerance limit should have meaning for SL

            if (type == 'SU') {
                ## xnall <- gamma + delta * asinh( (xall - xi)/ lambda )
                xtol_lower = xi + lambda * sinh( (xntol_lower - gamma)/delta )
                xtol_upper = xi + lambda * sinh( (xntol_upper - gamma)/delta )

            } else if (type == 'SB') {
                ## xnall <- gamma + delta *   log( (xall - xi)/(lambda + xi - xall) )
                zero <- function(x, gamma, delta, xi, xntol) {
                    ## find where xntol_calc == xntol
                    xntol_calc <- gamma + delta *   log( (xntol - xi)/(lambda + xi - xntol) )
                    zero <- xntol_calc - xntol
                }
                xtol_lower <- newton.raphson(zero, parms=jparms, xntol=xntol_lower, xguess=mean(xall))
                xtol_upper <- NA
                                                                                           
            } else {
                ## type == 'SL'
                ## xnall <- gamma + delta *   log( (xall - xi)/ lambda )
                xtol_lower <- NA
                xtol_upper <- xi + lambda * exp( (xntol_upper - gamma)/delta )
                
            }

        }
    }
    bounds <- data.frame(sided = c(side, side),
                         ztol  = c(ztol_lower, ztol_upper),
                         xntol = c(xntol_lower, xntol_upper),
                         xtol  = c(xtol_lower, xtol_upper))
    return(list(xz=df, jparms=jparms, ztol_out, bounds=bounds))
}

## ## test
##
## set.seed(1)
## x                          <- rnorm(n=1E5, mean=10, sd=1)
## tol_out_john               <- johnson_tol(x, alpha=0.1, P=0.99, side=1)
## jparms                     <- tol_out_john$jparms
## ztol_upper                <- tol_out_john$ztol_upper
## ztol_upper
## upper_tolerance_limit_john <- tol_out_john$xtol_upper$`root approximation`
## upper_tolerance_limit_john
## 
## xz <- tol_out_john$xz
## plotspace(1,3)
## hist(xz$x)
## hist(xz$z)
## plot(xz$x, xz$z)
