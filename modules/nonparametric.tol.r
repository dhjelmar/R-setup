nonparametric.tol <- function(x, conf=NA, P=NA, tol.index=NA, truncate=FALSE) {
    ## given x and any 2 of the following parameters: conf, P, and tol.index
    ## find the 3rd parameter

    ## parameters: x         = observations in a dataframe or vector (or number of obser ations)
    ##             conf      = desired confidence
    ##             P         = desired coverage
    ##             tol.index = index of observation for the desired tolerance bound
    ##
    ## options:    truncate  = FALSE (default) reports calculated parameter without truncation
    ##                       = FALSE truncates calcualted parameter to 2 decimal places

    ## convert x to a vector then sort from low to high
    if (is.data.frame(x)) x <- x[[1]]
    if (length(x) == 1)   x <- 1:x
    x <- sort(x)
    num <- length(x)

    if (is.na(tol.index)) {
        ## confidence and P specified, find tolerance bound

        alpha <- 1 - conf
        if (num < log(alpha)/log(P)) {
            message <- paste('ERROR: Insufficient number of points for requested confidence ',
                             'and coverage, P. Need at least ', ceiling(log(alpha)/log(P)),
                             ' points for ', conf*100, '% confidence on ', P*100, 'th percentile.',
                             sep='')
            return(list(message=message, df=NULL)) # exit function
        }

        lhs <- conf

        ## start with highest index and march down until find desired confidence
        for (i in num:1) {

            ## calculated confidence for index i
            rhs <- stats::pbinom(i-1, num, P)

            if (rhs < lhs) {
                ## index i confidence value dropped below the target confidence
                ## so base tolerance on the prior, higher index, i+1
                tol.index <- i+1
                tol <- x[tol.index]
                ## round tolerance down keeping 2 decimal places
                if (isTRUE(truncate)) tol <- floor(tol * 10000) / 10000
                break # exit for loop
            }

        }

    } else if (is.na(conf)) {
        ## P and index of tolearnce bound specified, find confidence
        tol <- x[tol.index]
        conf <- stats::pbinom(tol.index-1, num, P)
        if (isTRUE(truncate)) conf <- floor(conf * 10000) / 10000

    } else if (is.na(P)) {
        ## conf and index of tolearnce bound specified, find coverage (P)
        tol <- x[tol.index]
        alpha <- 1 - conf
        highpoints <- num - tol.index
        P <- stats::qbeta(alpha, num-highpoints, highpoints+1)
        if (isTRUE(truncate)) P<- floor(P * 10000) / 10000
        
    } else {
        message <- 'ERROR: One of the following must be NA: conf, P, tol.index'
        return(list(message=message, df=NULL)) # exit function
    }

    ## collect results
    df <- data.frame(observations=num, side='upper', sided=1, conf=conf, P=P,
                     tol.bound=tol, tol.index=tol.index, high.points=num-tol.index)
    message <- paste('For ', num, ' observations, the upper, 1-sided ',
                     conf*100, '% confidence bound on the ', P*100, 'th percentile is the ',
                     tol.index, 'th observation (', num-tol.index, ' high points) which is ', tol, '.',
                     sep='')
    return(list(message=message, df=df))
}

nonparametric.tol.test <- function() {
    test <- 'failure'
    if (nonparametric.tol(80, 0.9, 0.95)$df$tol.bound == 79) test <- 'success'
    cat('test 1 = ', test, '\n')
    
    test <- 'failure'
    if (nonparametric.tol(1:80, 0.9, 0.95)$df$tol.bound == 79) test <- 'success'
    cat('test 2 = ', test, '\n')
    
    test <- 'failure'
    if (nonparametric.tol(459, 0.99, 0.99)$df$tol.bound == 459) test <- 'success'
    cat('test 3 = ', test, '\n')
    
    test <- 'failure'
    if (nonparametric.tol(30000, 0.99, 0.99)$df$tol.bound == 29740) test <- 'success'
    cat('test 4 = ', test, '\n')
    
    test <- 'failure'
    if (is.null(nonparametric.tol(458, 0.99, 0.99)$df)) test <- 'success'
    cat('test 5 = ', test, '\n')
    
    test <- 'failure'
    if (is.null(nonparametric.tol(458, conf=0.99, P=0.99, tol.index=457)$df)) test <- 'success'
    cat('test 6 = ', test, '\n')    
    
    test <- 'failure'
    if (nonparametric.tol(1:80, conf=0.9, tol.index=79, truncate=TRUE)$df$P == 0.9522) test <- 'success'
    cat('test 7 = ', test, '\n')
    
    test <- 'failure'
    if (nonparametric.tol(1:80, P=0.95, tol.index=79, truncate=TRUE)$df$conf == 0.9139) test <- 'success'
    cat('test 8 = ', test, '\n')
}
