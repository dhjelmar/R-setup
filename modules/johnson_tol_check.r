johnson_tol_check <- function(x=NULL, jfit=NULL, breaks=NULL) {

    ## default parameters if not supplied
    if (is.null(x))      x      <- mtcars$mpg
    if (is.null(jfit))   jfit   <- SuppDists::JohnsonFit(x)
    if (is.null(breaks)) breaks <- length( hist(x, plot=FALSE)$breaks )
    
    ## get transformation data
    jparms <- jfit
    out    <- johnson_tol(x, jfit=jparms, alpha=0.1, P=0.99, side=1, plots='no')
    xn     <- out$xn
    xntol_upper <- out$xntol_upper
    xtol_upper  <- out$xtol_upper

    ## plot normal tranformation histogram
    plotspace(2,2)
    xnmin <- min(xn)
    xnmax <- max(xn, xntol_upper)
    xnrange <- seq(xnmin, xnmax, (xnmax-xnmin)/100)
    hist(xn, breaks=breaks, freq=FALSE, xlim=range(xnrange))
    xdensity_norm <- stats::dnorm(xnrange, mean(xn), sd(xn))
    lines(xnrange, xdensity_norm)
    abline(v=xntol_upper)

    ## plot normal transformation qqplot
    qualityTools::qqPlot(xn)

    ## plot johnson histogram
    out <- hist_nwj(x, jfit=jparms, breaks=breaks, type='j')

    ## plot johnson qqplot
    out <- qqplot_nwj(x, type='j', jfit=jparms)

    return(list(x=x, jfit=jfit))
    
}
## both of the following return the same thing:
## johnson_tol(mtcars$mpg, plots=TRUE)
## johnson_tol_check(mtcars$mpg)

## need to document:
##    hqt
##    hist_nwj
##    qqplot_nwj
##    johnson_tol
##    johnson_tol_check
