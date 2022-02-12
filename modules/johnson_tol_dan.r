johnson_tol <- function(x, jfit='all', alpha=0.01, P=0.99, side=1, plots=FALSE, breaks=NULL) {

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
    ## Calculate tolerance limits if JohnsonSU fit
    if (type == 'SU') {

        ## transform x to normal distribution, xn
        ## https://www.sigmamagic.com/blogs/how-do-i-transform-data-to-normal-distribution/
        xn <- gamma + delta * asinh( (x - xi)/ lambda )
        df <- data.frame(x, xn)
        
    }
    
    return(list(sided=side, 
                alpha=alpha,
                P=P,
                jfit=jfit,
                jparms=jparms,
                xn=xn,
                xntol_lower=xntol_lower,
                xntol_upper=xntol_upper,
                x=x,
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

## johnson_tol(mtcars$mpg, plots=TRUE)

x <- mtcars$mpg

## find best fit JohnsonSU paramters
jparms.out <- ExtDist::eJohnsonSU(x)  # this failed
jparms <- list(gamma   = jparms.out$gamma,
               delta   = jparms.out$delta,
               xi      = jparms.out$xi,
               lambda  = jparms.out$lambda,
               type    = 'SU')
jparms  <- SuppDists::JohnsonFit(x) # returns SU
gamma <- jparms$gamma
delta <- jparms$delta
xi    <- jparms$xi
lambda <- jparms$lambda

## find Pth quantile for the best fit distribution 
## P = coverage (e.g., 99%)
P <- 0.99
quantile.P <- ExtDist::qJohnsonSU(P, params=jparms)

## find Z assocaited with some % confidence on standard normal
alpha <- 0.01
sided <- 2
z.conf <- qnorm(1 - alpha/sided)   # confidence = 1 - alpha; z.conf = 2.575829 for 1-sided 99

## replace gamma 
gamma.P.alpha <- z.conf - delta * asinh( (quantile.P - xi)/ lambda )

## refit the Johnson distributions
jparms  <- SuppDists::JohnsonFit(gamma.P.alpha) # returns SU

