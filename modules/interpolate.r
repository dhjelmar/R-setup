interpolate <- function(x, xvec, yvec=seq(1, length(xvec)), plot=FALSE) {
    ## interpolate to find y that corresponds with x for given xvec and yvec vectors
    ## extrapolate if x is outside the bounds of xvec

    if (length(xvec) == 1) {
        ## no interpolation is needed so just replace x with yvec
        y <- yvec[1]

    } else if (x < min(xvec)) {
        ## extrapolate on low side to find y
        x1 <- xvec[1]
        x2 <- xvec[2]
        y1 <- yvec[1]
        y2 <- yvec[2]
        y <- (y2-y1)/(x2-x1)*(x-x1) + y1

    } else if (x > max(xvec)) {
        ## extrapolate on high side to find y
        num <- length(xvec)
        x1 <- xvec[num-1]
        x2 <- xvec[num]
        y1 <- yvec[num-1]
        y2 <- yvec[num]
        y <- (y2-y1)/(x2-x1)*(x-x1) + y1

    } else {
        ## interpolate to find y
        ## find index of x in xvec; return lower and higher indices if not found exactly
        xlow  <-  last(which(x >= xvec))
        xhigh <- first(which(x <= xvec))

        if (xlow == xhigh) {
            ## exact match so replace with scaled value
            y <- yvec[xlow]

        } else {
            ## x found between two values in xvec
            x1 <- xvec[xlow]
            x2 <- xvec[xhigh]
            y1 <- yvec[xlow]
            y2 <- yvec[xhigh]
            y <- (y2-y1)/(x2-x1)*(x-x1) + y1
        }

    }

    if (isTRUE(plot)) {
        plot(xvec, yvec, type='b', xlim=range(xvec,x), ylim=range(yvec,y))
        points(x,y,col='red')
    }

    return(y)
}

interpolate_test <- function() {
    plotspace(2,3)
    xvec <- c(5, 10, 20, 40, 80, 160)
    interpolate(3, xvec, plot=TRUE)
    interpolate(5, xvec, plot=TRUE)
    interpolate(30, xvec, plot=TRUE)
    interpolate(50, xvec, plot=TRUE)
    interpolate(160, xvec, plot=TRUE)
    interpolate(170, xvec, plot=TRUE)
    return()
}
