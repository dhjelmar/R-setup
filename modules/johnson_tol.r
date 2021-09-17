johnson_tol <- function(x, jfit='all', alpha=0.01, P=0.99, side=1, plots='no') {

    ## input:   x     = data set
    ##          alpha = 1 - confidence
    ##          P     = proportion = coverage (tolerance interval only)
    ##          side  = 1 = 1-sided tolerance limit
    ##                = 2 = 2-sided tolerance limit
    ##          jfit  = 'all' = uses SuppDists::JohnsonFit to determine parameters
    ##                = 'SU'  = uses ExtDist::eJohnsonSU to determine parameters
    ##                = list of user specified parameters
    ##                  e.g., jparms <- list(gamma   = -1.039,
    ##                                       delta   = 1.66,
    ##                                       xi      = 14.46,
    ##                                       lambda  = 6.95,
    ##                                       type    = 'SU',
    ##                       johnson_tol(mtcars$mpg, jfit=jparms)
    ##
    ## output:  jparms = Johnson fit parameters used in tolerance limit calculation
    ##          xtol_lower = lower bound tolerance limit
    ##          xtol_upper = upper bound tolerance limit
    
    proportion <- P

    ##--------------------------------------------------------
    ## initialize tolerance limits
    xntol_lower  = NA
    xntol_upper  = NA
    xtol_lower   = NA
    xtol_upper   = NA

    
    ##--------------------------------------------------------
    ## determine Johnson parameters
    
    if (jfit[1] == 'all') {
        ## let R figure out which Johnson distribution fits best
        jparms  <- SuppDists::JohnsonFit(x)
    } else if (jfit[1] == 'SU') {
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
        ## calculate normal tolerance limit, xntol
    
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
                x=x,
                xn=xn,
                xtol_lower=xtol_lower, 
                xtol_upper=xtol_upper))
}

## ## test

## johnson_tol(mtcars$mpg)

## out <- johnson_tol(mtcars$mpg)
## upper_tolerance_limit_john <- out$xtol_upper
## upper_tolerance_limit_john
## jparms <- out$jparms

## johnson_tol(mtcars$mpg, jfit=jparms)
