newton.raphson.target <- function(fun, ytarget=0, x0, maxiter=1000, tol=1E-5, ...) {
    ## wrapper to allow specification of a target value other than zero for function, fun

    ## input: ytarget = target value for function
    ##        see documentation for pracma::newtonRaphson() for other input descriptions

    zero <- function(...) {
        fun(...) - ytarget
    }

    out <- pracma::newtonRaphson(fun=zero, x0=x0, maxiter=maxiter, tol=tol, ...)
}


newton.raphson.target_test <- function() {
    ##----------------------
    ## FIRST TEST HOME GROWN newton.raphson()

    ## following correctly finds closest root (i.e., zero crossing in this case)
    plot(sin, 0, 5*pi)
    ytarget = 0
    abline(h=ytarget, col='red')
    xguess = 5
    points(xguess, ytarget, col='blue', pch=16)
    out <- newton.raphson(sin, xguess=xguess, ytarget=ytarget)
    root <- out$root
    abline(v=root, col='red')

    ## following misses closest root but does find a root
    plot(sin, 0, 5*pi)
    ytarget = 0.5
    abline(h=ytarget, col='red')
    xguess = 5
    points(xguess, ytarget, col='blue', pch=16)
    out <- newton.raphson(sin, xguess=xguess, ytarget=ytarget)
    ## plot intermediate attempts
    points(out$xiterations, out$yiterations, col='blue', pch=1)
    root <- out$root
    abline(v=root, col='red')


    ##----------------------
    ## NEXT TEST WRAPPER ON pracma::newtonRaphson()

    ## following similarly misses closest root but does find a root
    plot(sin, 0, 5*pi)
    ytarget = 0.5
    abline(h=ytarget, col='red')
    xguess = 5
    points(xguess, ytarget, col='blue', pch=16)
    out <- newton.raphson.target(sin, x0=xguess, ytarget=ytarget)
    root <- out$root
    abline(v=root, col='red')
}
