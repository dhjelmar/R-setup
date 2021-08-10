newton.raphson <- function(f, ..., xguess=0, tol = 1e-5, n = 1000) {
    ## use newton raphson to find x where f(x, ...) = 0
    ## where "..." is a list of additional arguments needed by f, if any
    ## starts search from xguess

    args <- unlist( list(...) )
    
    ## following modified from https://rpubs.com/aaronsc32/newton-raphson-method
    library(numDeriv) # Package for computing f'(x)

    x0 <- xguess # Set start value to supplied guess
    k <- n       # Initialize for iteration results

    ## Check to see if xguess result in 0
    ## if (f(x0, ...) == 0.0) return(x0)
    if (length(args) == 1) 
    
    ## iterate to find where f(x, ...) = 0
    for (i in 1:n) {
        dx <- genD(func = f, ..., x = x0)$D[1] # First-order derivative f'(x0)
        x1 <- x0 - (f(x0, ...) / dx)           # Calculate next guess x1
        k[i] <- x1                        # Store x1
        ## Once the difference between x0 and x1 becomes sufficiently small, output the results.
        if (abs(x1 - x0) < tol) {
            root.approx <- tail(k, n=1)
            res <- list('root' = root.approx, 'iterations' = k)
            return(res)
        }
        ## If Newton-Raphson has not yet reached convergence set x1 as x0 and continue
        x0 <- x1
    }
    print('Too many iterations in method')
}

## example: this works
## x <- rnorm(1000)
## jparms  <- JohnsonFit(x)
## zero <- function(x, parms, z)  john_z(x, parms) - z
## x.out <- newton.raphson(
##     zero,
##     parms = jparms,
##     z=2.4,
##     xguess = mean(x),
##     tol=1E-10
## )
## x.out$root

## ## example: this works
## x <- rnorm(1000)
## zero <- function(x, mean, sd, z)  (x - mean)/sd - z
## x.out <- newton.raphson(
##     zero,
##     mean=mean(x),
##     sd=sd(x),
##     z=2.4,
##     xguess = mean(x),
##     tol=1E-10
## )
## x.out$root

## ## example: this fails because only 1 x is given to z at a time so sd(x) is NA
## ##          (mean(x) = x but that does not make the function abort
## x <- rnorm(1000)
## zero <- function(x, z)  (x - mean(x))/sd(x) - z
## x.out <- newton.raphson(
##     zero,
##     z=2.4,
##     xguess = mean(x),
##     tol=1E-10
## )
## x.out$root
