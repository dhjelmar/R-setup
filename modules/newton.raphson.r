newton.raphson <- function(f, ..., xguess=0, ytarget=0, tol = 1e-5, n = 1000,
                           relax=0.8, nrelax=30, plots=FALSE, plot.add=FALSE, main.adder=NULL) {
    ## use newton raphson to find x where f(x, ...) = ytarget
    ## where "..." is a list of additional arguments needed by f, if any
    ## starts search from xguess
    ## relax = relaxation factor to help convergence after 1st 100 iterations
    ##       = 1 is no relaxation
    ##       = 0 is no change (so relax must be > 0)

    args <- unlist( list(...) )
    
    ## following modified from https://rpubs.com/aaronsc32/newton-raphson-method
    library(numDeriv) # Package for computing f'(x)

    x0 <- xguess      # Set start value to supplied guess
    xvalue <- n       # Initialize for iteration results
    yvalue <- n       # Initialize for iteration results

    ## Check to see if xguess result in 0
    ## ## if (f(x0, ...) == 0.0) return(x0)
    ## if (length(args) == 1) {
    ##     if (f(x0) == ytarget) return(x0)
    ## } else {
    if (f(x0, ...) == ytarget) return(x0)
    ## }
    
    ## iterate to find where f(x, ...) = ytarget
    xvalue[1] <- x0          # store x values
    yvalue[1] <- f(x0, ...)  # store y values

    tried <- 'initial guess'
    for (i in 2:(n+1)) {

        ## function value at current guess for x
        y0 <- f(x0, ...)
        cat('Newton-Raphson iteration =', i-1, 'x =', x0, ', y =', y0, ', yerror =', y0-ytarget, 'tried =', tried, '\n')

        ## Use first order derivative to make next guess, x1
        dx <- genD(func = f, ..., x = x0)$D[1] # First-order derivative f'(x0)
        x1 <- x0 - ( (y0 - ytarget) / dx)
        
        ## save new values of x and y
        xvalue[i] <- x1
        y1      <- f(x1, ...)
        yvalue[i] <- y1
        if (plots == TRUE & plot.add == TRUE) points(x1, y1, col='red', pch=3)

        ## if (abs(x1 - x0) < tol & abs(y1 - y0) < tol) {
        ##     ## differences between x0 and x1 and between y0 and y1 are sufficiently small,
        ##     ## so output the results
        ## if (abs(x1 - x0) < tol & abs(y1 - ytarget) < tol) {
        ##     ## differences between x0 and x1 and between y1 and ytarget are sufficiently small,
        if (abs(y1 - ytarget) < tol) {
            ## difference between y1 and ytarget is sufficiently small,
            ## so output the results
            if (plots == TRUE & plot.add == TRUE) points(x1, y1, col='red', pch=2)
            cat('Newton-Raphson iteration =', i, 'x =', x1, ', y =', y1, ', yerror =', y1-ytarget, 'tried =', tried, '\n')
            ## root.approx <- tail(xvalue, n=1)
            ## res <- list('root' = root.approx, 'iterations' = xvalue)
            res <- list(xiterations=xvalue, yiterations=yvalue, ytarget=ytarget, root=x1, convergence="successful")
            if (isTRUE(plots)) {
                if (is.null(main.adder)) {
                    main <- 'Solution found within tolerance in newton.raphson()'
                } else {
                    main <- paste('Solution found within tolerance in newton.raphson(); ',
                                  main.adder, sep='')
                }
                plot(xvalue, yvalue, col='red', main=main)
                points(xvalue[1], yvalue[1], col='black', pch=19)  # solid black initial guess
                last <- length(xvalue)
                points(xvalue[last], yvalue[last], col='red', pch=19)  # solid red final result
                xplot <- seq(min(xvalue), max(xvalue), (max(xvalue)-min(xvalue))/100)
                ## yplot <- f(xplot, ...) # this was OK for most functions but not all
                yplot <- NA
                for (i in 1:length(xplot)) {
                    yplot[i] <- f(xplot[i], ...)
                }
                points(xplot,yplot,type='l')    # black line for the function
                abline(h=ytarget, col='blue')   # blue horizontal line for the target
                legend('topright',
                       legend=c('initial guess', 'intermediate steps', 'final result', 'function', 'ytarget'),
                       col   =c('black',         'red',                'red',          'black',    'blue'),
                       pch   =c(19,              1,                    19,              NA,         NA),
                       lty   =c(NA,              NA,                   NA,              1,          1))
            }
            return(res)
        }

        ## check whether solution is bracketed
        if ( (yvalue[i]-ytarget) / (yvalue[i-1]-ytarget) < 0) {
            ## bracketed solution since successive y values - ytarget have opposite sign
            ## interpolate to new x0 guess
            xnext <- (x1 - x0)/(yvalue[i] - yvalue[i-1])*(ytarget - yvalue[i-1]) + x0
            tried <- 'interpolation'
        } else {
            ## have not bracketed solution
            xnext <- x1
            tried <- 'slope used'
        }

        if (i > nrelax) {
            ## have already tried many iterations, so relax next guess
            ## added a little variablility to relaxation factor to help keep it from getting stuck
            relax.use <- relax*rnorm(1, 1, 0.01)
            xnext <- x0 + relax.use * (xnext-x0)
            tried <- paste(tried, ' with relaxation', sep='')
        }

        x0 <- xnext
        
    }

    ## below here will not be reached if the search is successful
    if (is.null(main.adder)) {
        main <- 'Too many iterations in newton.raphson()'
    } else {
        main <- paste('Too many iterations in newton.raphson(); ', main.adder, sep='')
    }
    plot(xvalue, yvalue, col='red', main=main)
    points(xvalue[1], yvalue[1], col='black', pch=19)  # solid black filled circle initial guess
    last <- length(xvalue)
    points(xvalue[last], yvalue[last], col='red', pch=19)  # solid red filled circle final result
    xplot <- seq(min(xvalue), max(xvalue), (max(xvalue)-min(xvalue))/100)
    ## yplot <- f(xplot, ...) # this was OK for most functions but not all
    yplot <- NA
    for (i in 1:length(xplot)) {
        yplot[i] <- f(xplot[i], ...)
    }
    points(xplot,yplot,type='l')    # black line for function
    abline(h=ytarget, col='blue')   # blue horizontal line for the target
    ## xplot1; points(xplot1, f(xplot1, ...), col='blue')
    legend('topright',
           legend=c('initial guess', 'intermediate steps', 'final result', 'function', 'ytarget'),
           col   =c('black',         'red',                'red',          'black',    'blue'),
           pch   =c(19,              1,                    19,              NA,         NA),
           lty   =c(NA,              NA,                   NA,              1,          1))
    res <- list(xiterations=xvalue, yiterations=yvalue, ytarget=ytarget, root=NA, convergence='failed')
    return(res)

}

test_newton.raphson <- function() {
    set.seed(1)
  
    ## example: this works
    x <- rnorm(1000)
    jparms  <- SuppDists::JohnsonFit(x)
    zero <- function(x, parms, z)  john_z(x, parms) - z
    x.out <- newton.raphson(
        zero,
        parms = jparms,
        z=2.4,
        xguess = mean(x),
        tol=1E-10,
        plots=TRUE
    )
    x.out$root
    
    ## example: this works
    set.seed(1)
    x <- rnorm(1000, mean=10, sd=20)
    zero <- function(x, parm1, parm2, z)  (x^2 + parm1 + parm2) - z
    plotspace(2,1)
    x.out <- newton.raphson(
      zero,
      parm1 = 11,
      parm2 = mean(x),
      z=400,
      xguess = 30,
      tol=1E-10,
      plots=TRUE
    )
    x.out$root
    
    ## example: same as above, but ytarget is not zero (this works)
    set.seed(1)
    x <- rnorm(1000, mean=10, sd=20)
    notzero <- function(x, parm1, parm2)  (x^2 + parm1 + parm2)
    plotspace(2,1)
    x.out <- newton.raphson(
      notzero,
      xguess = 30,
      ytarget=400,
      parm1 = 11,
      parm2 = mean(x),
      tol=1E-10,
      plots=TRUE
    )
    x.out$root
    ## plot(x, y = x^2 + 11 + mean(x), main='wider view of function')
    curve(x^2 + 11 + mean(x), -10, 32, main='wider view of function')
    x.iter <- x.out$xiterations
    y.iter <- x.out$yiterations
    num <- length(x.iter)
    points(x.iter     , y.iter     , col='red' , pch=1 , cex=2)
    points(x.iter[1]  , y.iter[1]  , col='blue', pch=16, cex=2)
    points(x.iter[num], y.iter[num], col='red' , pch=16, cex=2)
    
    ## example: this works
    x <- rnorm(1000)
    zero <- function(x, mean, sd, z)  (x - mean)/sd - z
    x.out <- newton.raphson(
        zero,
        mean=mean(x),
        sd=sd(x),
        z=2.4,
        xguess = mean(x),
        tol=1E-10,
        plots=TRUE
    )
    x.out$root

    ## example: This fails because only 1 xguess is passed into the function so sd(xguess) is NA
    ##          It aborts when the newton.raphson function tries using the zero function
    ##          e.g., zero(x=2, z=3) returns NA
    ##                zero(x=c(2,4), z=3) returns -3.71 and -2.29
    ##          Note that mean(x) = x also does not make much sense for a single x value, 
    ##          but that does not make the function abort
    x <- rnorm(1000)
    zero <- function(x, z)  (x - mean(x))/sd(x) - z
    x.out <- NA
    x.out <- newton.raphson(
        zero,
        z=2.4,
        xguess = mean(x),
        tol=1E-10,
        plots=TRUE
    )
    x.out$root
}
