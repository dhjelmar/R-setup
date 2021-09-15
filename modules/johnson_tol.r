johnson_tol <- function(x, alpha=0.01, P=0.99, side=1, jfit='all', plots='no') {

    ## x = data set
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
        jparms  <- SuppDists::JohnsonFit(x)
    } else if (jfit == 'SU') {
        ## force the Johnson SU distribution
        jparms.out <- ExtDist::eJohnsonSU(x)
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
    ## transform x to normal distribution, xn
    ## https://www.sigmamagic.com/blogs/how-do-i-transform-data-to-normal-distribution/

    if (type == 'SU') {
        xn <- gamma + delta * asinh( (x - xi)/ lambda )
    } else if (type == 'SB') {
        xn <- try( gamma + delta *   log( (x - xi)/(lambda + xi - x) ) )
    } else {
        ## type == 'SL'
        xn <- try( gamma + delta *   log( (x - xi)/ lambda ) )  
    }
    df <- data.frame(x, xn)
    
    if (plots != 'no') {
        plotspace(1,2)
        hist(x)
        hist(xn)
    }
    
    ##--------------------------------------------------------
    ## error check then proceed
    na_rows <- df[is.na(df$xn),]
    if (nrow(na_rows) != 0) {

        ## z values could not be calculated for all x values
        cat('\n')
        cat('#####################################\n')
        cat('  ERROR IN JOHNSON_TOL  \n')
        cat('  xn(x) = NA or NaN for some x values \n')
        cat('#####################################\n')
        print(jparms)
        print(na_rows)
        xntol_out <- tolerance::normtol.int(xn, alpha = alpha, P=proportion, side=side)

    } else {
        ## johnson fit returned xn values for all x values so possibly a decent fit
        
        ##--------------------------------------------------------
        ## calculate normal tolerance limit, ztol (i.e., for standard normal)
    
        xntol_out <- tolerance::normtol.int(xn, alpha = alpha, P=proportion, side=side)
        if (side == 1) {
            xntol_lower <- xntol_out$'1-sided.lower'
            xntol_upper <- xntol_out$'1-sided.upper'
        } else if (side == 2) {
            xntol_lower <- xntol_out$'2-sided.lower'
            xntol_upper <- xntol_out$'2-sided.upper'
        }

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
                ## xn <- gamma + delta * asinh( (x - xi)/ lambda )
                xtol_lower = xi + lambda * sinh( (xntol_lower - gamma)/delta )
                xtol_upper = xi + lambda * sinh( (xntol_upper - gamma)/delta )

            } else if (type == 'SB') {
                ## xn <- gamma + delta *   log( (x - xi)/(lambda + xi - x) )
                zero <- function(x, gamma, delta, xi, xntol) {
                    ## find where xntol_calc == xntol
                    xntol_calc <- gamma + delta *   log( (xntol - xi)/(lambda + xi - xntol) )
                    zero <- xntol_calc - xntol
                }
                xtol_lower <- newton.raphson(zero, parms=jparms, xntol=xntol_lower, xguess=mean(x))
                xtol_upper <- NA
                                                                                           
            } else {
                ## type == 'SL'
                ## xn <- gamma + delta *   log( (x - xi)/ lambda )
                xtol_lower <- NA
                xtol_upper <- xi + lambda * exp( (xntol_upper - gamma)/delta )
                
            }

        }
    }
    return(list(sided=side, 
                alpha=alpha,
                P=P,
                jfit=jfit,
                jparms=jparms, 
                xtol_lower=xtol_lower, 
                xtol_upper=xtol_upper))
}

## ## test

## johnson_tol(mtcars$mpg)

## out <- johnson_tol(mtcars$mpg)
## upper_tolerance_limit_john <- out$xtol_upper
## upper_tolerance_limit_john

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
