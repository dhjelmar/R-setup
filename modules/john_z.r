john_z <- function(x, parms) {
    ## transform x values to z which is normally distributed if Johnson distribution
    ## x can be one value or a set of values
    ## parms is the set of parameters from SuppDists::JohnsonFit on the entire dataset
    u <- (x - parms$xi) / parms$lambda
    type <- parms$type
    if (type == 'SU') {
        ## unbounded distribution
        fu <- u + (1 + u^2)^0.5
    } else if (type == 'SL') {
        ## log normal
        fu <- u
    } else if (type == 'SB') {
        ## bounded distribution
        fu <- u / (1-u)
    } else if (type == 'SN') {
        ## normal distribution
        fu <- exp(u)
    } else {
        cat('\n')
        cat('#############################\n')
        cat('  ERROR                      \n')
        cat('  no type provided in parms  \n')
        cat('#############################\n')
        print(jparms)
        ## stop()
    }
    ## z is normally distributed
    z    <- parms$gamma + parms$delta * log(fu)
    return(z)
}
