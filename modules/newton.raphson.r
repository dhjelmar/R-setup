newton.raphson <- function(f, ..., xguess=0, ytarget=0, tol = 1e-5, n = 1000, plot='no') {
    ## use newton raphson to find x where f(x, ...) = ytarget
    ## where "..." is a list of additional arguments needed by f, if any
    ## starts search from xguess

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

    for (i in 2:(n+1)) {

        ## Use first order derivative to make next guess, x1
        dx <- genD(func = f, ..., x = x0)$D[1] # First-order derivative f'(x0)
        x1 <- x0 - ( (f(x0, ...)-ytarget) / dx)

        if (abs(x1 - x0) < tol) {
            ## difference between x0 and x1 is sufficiently small, so output the results
            root.approx <- tail(xvalue, n=1)
            res <- list('root' = root.approx, 'iterations' = xvalue)
            if (plot != 'no') {
                main <- 'Solution found within tolerance'
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
                points(xplot,yplot,type='l')
                legend('bottomright',
                       legend=c('initial guess', 'intermediate steps', 'final result'),
                       col   =c('black',         'red',                'red'),
                       pch   =c(19,              1,                    19))
            }
            return(res)
        }
        
        ## save new valeus of x and y
        xvalue[i] <- x1
        yvalue[i] <- f(x1, ...)
        
        ## check whether solution is bracketed
        if ( (yvalue[i]-ytarget) / (yvalue[i-1]-ytarget) < 0) {
            ## bracketed solution since successive y values have opposite sign
            ## use bisection to keep solution from diverging
            x0 <- (x0 + x1)/2
        } else {
            ## have not bracketed solution so use X1 as next guess
            x0 <- x1
        }
    }
    main <- 'Too many iterations in method'
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
    points(xplot,yplot,type='l')
    legend('bottomright',
           legend=c('initial guess', 'intermediate steps', 'final result'),
           col   =c('black',         'red',                'red'),
           pch   =c(19,              1,                    19))
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
        plot='yes'
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
      plot='yes'
    )
    x.out$root
    ## plot(x, y = x^2 + 11 + mean(x), main='wider view of function')
    curve(x^2 + 11 + mean(x), -10, 32, main='wider view of function')
    x.iter <- x.out$iterations
    y.iter <- x.iter^2 + 11 + mean(x)
    num <- length(x.iter)
    points(x.iter     , y.iter     , col='red' , pch=1 , cex=2)
    points(x.iter[1]  , y.iter[1]  , col='blue', pch=16, cex=2)
    points(x.iter[num], y.iter[num], col='red' , pch=16, cex=2)
    
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
      plot='yes'
    )
    x.out$root
    ## plot(x, y = x^2 + 11 + mean(x), main='wider view of function')
    curve(x^2 + 11 + mean(x), -10, 32, main='wider view of function')
    x.iter <- x.out$iterations
    y.iter <- x.iter^2 + 11 + mean(x)
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
        plot='yes'
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
    x.out <- newton.raphson(
        zero,
        z=2.4,
        xguess = mean(x),
        tol=1E-10,
        plot='yes'
    )
    x.out$root
}
