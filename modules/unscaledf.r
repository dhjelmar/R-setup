unscaledf <- function(dfs, meanvec, scaledata, scalefactor, zeromin) {
    ## reverses scaling
    if (zeromin == TRUE) dfs <- sweep(dfs, 2, (scalefactor/2), "-")
    scalerange <- (sapply(scaledata,max) - sapply(scaledata,min)) / scalefactor
    df <- sweep(dfs, 2, scalerange, "*")
    df <- sweep(df,  2, meanvec,    "+")
    df <- data.frame(df)
    ## names(df) <- gsub("_s", "", names(df))
    return(df)
}
## dfs.out <- scaled(dfs)
## dfs     <- dfs.out$dfs
## dfsu    <- unscaledf(dfs, meanvec=dfs.out$meanvec, scaledata=dfs.out$scaledata,
##                      scalefactor=dfs.out$scalefactor, zeromin=dfs.out$zeromin)

unscaled_model <- function(cs, y, meanvec, scaledata, scalefactor, zeromin, data=FALSE) {
    ## input: cs = list of scaled model coefficient values including intercept as 1st
    ##             (same order and number as when data was scaled)
    ##         y  = name of independent variable (e.g., 'mpg')
    ##        additional parameters are original scaledf() output
    ## note:  function assumes dependent and independent model parameters were scaled
    bs <- cs[1]
    cs <- cs[2:length(cs)]
    scalerange <- (sapply(scaledata,max) - sapply(scaledata,min)) / scalefactor
    ## identify vector location with independent variable
    yloc <- which(names(scalerange) == y)
    xloc <- -yloc
    ## unscale model intercept
    if (zeromin == FALSE) scalefactor <- 1
    a <- -meanvec[xloc]/scalerange[xloc] + scalefactor/2
    b <- (bs + cs %*% as.matrix(a) - scalefactor/2) * scalerange[yloc] + meanvec[yloc]
    b <- b[[1]]
    ## unscale model coefficients
    coef <- cs[xloc] * scalerange[yloc] / scalerange[xloc]
    ##
    fit <- 'fit <- b + as.matrix(xdata) %*% as.matrix(coef)'
    return(list(b=b, coef=coef, fit=fit))
}
