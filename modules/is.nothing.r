is.nothing <- function(x) {
    ## return true if x is NULL or NA
    any(is.null(x)) | any(is.na(x))
}
is.nothing.test <- function() {
    v1 <- 'failure'
    a <- NULL
    if (isTRUE(is.nothing(a))) v1 <- 'success'
    v2 <- 'failure'
    a <- NA
    if (isTRUE(is.nothing(a))) v2 <- 'success'
    v3 <- 'failure'
    a <- NaN
    if (isTRUE(is.nothing(a))) v3 <- 'success'
    v4 <- 'failure'
    a <- 1
    if (isFALSE(is.nothing(a))) v4 <- 'success'
    v5 <- 'failure'
    a <- c(1,2,3,4)
    if (isFALSE(is.nothing(a))) v5 <- 'success'
    v6 <- 'failure'
    a <- c('asdf', '1234')
    if (isFALSE(is.nothing(a))) v6 <- 'success'
    v7 <- 'failure'
    a <- list(one=1, two=2)
    if (isFALSE(is.nothing(a))) v7 <- 'success'
    v8 <- 'failure'
    a <- data.frame(one=c(1,2,3), two=c(10,20,30))
    if (isFALSE(is.nothing(a))) v8 <- 'success'
    cat('verif 1 =', v1, '\n')
    cat('verif 2 =', v2, '\n')
    cat('verif 3 =', v3, '\n')
    cat('verif 4 =', v4, '\n')
    cat('verif 5 =', v5, '\n')
    cat('verif 6 =', v6, '\n')
    cat('verif 7 =', v7, '\n')
    cat('verif 8 =', v8, '\n')
}
