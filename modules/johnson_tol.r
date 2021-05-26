johnson_tol <- function(xall, alpha=0.01, P=0.99, side=1) {

    ## xall = data set
    ## alpha  = 1 - confidence
    ## P = pvalue = coverage (tolerance interval only)
    pvalue <- P
    
    ## install.packages('tolerance')
    library(tolerance)
    ## install.packages('SuppDists')
    library(SuppDists)

    ## Johnson distribution calculations from SuppDists package
    jparms  <- JohnsonFit(xall)
    ## transform xall dataset to normally distributed z using Johnson distribution
    zall <- john_z(xall, jparms)
    df <- as_tibble( data.frame(x=xall, z=zall) )

    na_rows <- df[is.na(df$z),]
    if (nrow(na_rows) != 0) {

        ## z values could not be calculated for all x values
        cat('\n')
        cat('##############################\n')
        cat('  ERROR IN JOHNSON_TOL  \n')
        cat('  z(x) = NA for some x values \n')
        cat('##############################\n')
        print(jparms)
        print(na_rows)
        ztol_out <- normtol.int(zall, alpha = alpha, P=pvalue, side=side)
        ztol_upper <- NA
        xtol_upper <- NA
        
    } else {
        
        ## johnson fit returned z values for all x values so possibly a decent fit

        ## xz_plot <- plot(xall, zall)
        
        ## define function to find x for a given z
        ## zero <- function(x, parms=jparms, ztarget) john_z(x, parms) - ztarget
        zero <- function(x, parms, z) john_z(x, parms) - z
        
        ## check use of uniroot
        ## xcheck <- 13
        ## ztarget <- john_z(xcheck)
        ## ztox <- uniroot(zero, lower=min(xall), tol=1E-10, upper=2*max(xall), ztarget=ztarget)$root
        ## x_from_z <- uniroot(zero, ztarget=ztarget,
        ##                     lower=min(xall), upper=2*max(xall), tol=1E-10)$root
        
        ## find tolerance limits using normal distribution calculations on z
        ztol_out <- normtol.int(zall, alpha = alpha, P=pvalue, side=side)
        xmin <- min(xall)
        xmax <- max(xall)
        if (side == 1) {
            ztol_upper <- ztol_out$'1-sided.upper'
            ztol_lower <- ztol_out$'1-sided.lower'
        } else if (side == 2) {
            ztol_upper <- ztol_out$'2-sided.upper'
            ztol_lower <- ztol_out$'2-sided.lower'
        }

        ## transform z back to u
        ##  uniroot is not very robust so try newton raphson
        ##  xrange <- max(xall) - min(xall)
        ##  xtol_upper <- uniroot(zero, ztarget=ztol_upper,
        ##                        lower=xmax, upper=xmax + xrange, tol=1E-10)$root

        if (jparms$type == 'SB') {
            ## Not sure tolerance limits have meaning for a bounded fit
            xtol_lower <- NA
            xtol_upper <- NA
        } else {
            ## Tolerance limits should have meaning for SU and SN
            ## Upper tolerance limit should have meaning for SL
            xtol_upper <- newton.raphson(zero, parms=jparms, z=ztol_upper, xguess=mean(xall), tol=1E-10)
            if (jparms$type == 'SL') {
                xtol_lower <- NA
            } else {
                xtol_lower <- newton.raphson(zero, parms=jparms, z=ztol_lower, xguess=mean(xall), tol=1E-10)
            }
        }
    }
    return(list(xz=df, jparms=jparms, ztol_out=ztol_out, 
                ztol_upper=ztol_upper, xtol_upper=xtol_upper,
                ztol_lower=ztol_lower, xtol_lower=xtol_lower))
}

## ## test
##
## set.seed(1)
## x                          <- rnorm(n=1E5, mean=10, sd=1)
## tol_out_john               <- johnson_tol(x, alpha=0.1, P=0.99, side=1)
## jparms                     <- tol_out_john$jparms
## ztol_upper                <- tol_out_john$ztol_upper
## ztol_upper
## upper_tolerance_limit_john <- tol_out_john$xtol_upper$`root approximation`
## upper_tolerance_limit_john
## 
## xz <- tol_out_john$xz
## plotspace(1,3)
## hist(xz$x)
## hist(xz$z)
## plot(xz$x, xz$z)
