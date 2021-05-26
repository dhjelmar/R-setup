newton.raphson <- function(f, ..., xguess=0, tol = 1e-5, n = 1000) {
    ## use newton raphson to find x where f(x, ...) = 0
    ## where "..." is a list of additional arguments needed by f, if any
    ## starts search from xguess
    
    ## following modified from https://rpubs.com/aaronsc32/newton-raphson-method
    library(numDeriv) # Package for computing f'(x)

    x0 <- xguess # Set start value to supplied guess
    k <- n       # Initialize for iteration results

    ## Check to see if xguess result in 0
    if (f(x0, ...) == 0.0) return(x0)
    
    ## iterate to find where f(x, ...) = 0
    for (i in 1:n) {
        dx <- genD(func = f, ..., x = x0)$D[1] # First-order derivative f'(x0)
        x1 <- x0 - (f(x0, ...) / dx)           # Calculate next guess x1
        k[i] <- x1                        # Store x1
        ## Once the difference between x0 and x1 becomes sufficiently small, output the results.
        if (abs(x1 - x0) < tol) {
            root.approx <- tail(k, n=1)
            res <- list('root approximation' = root.approx, 'iterations' = k)
            return(res)
        }
        ## If Newton-Raphson has not yet reached convergence set x1 as x0 and continue
        x0 <- x1
    }
    print('Too many iterations in method')
}

##  xtol_upper <- newton.raphson(function(x, parms=jparms, z=ztol_upper) john_z(x, parms) - z,
##                               xguess = xmax, tol=1E-10)
