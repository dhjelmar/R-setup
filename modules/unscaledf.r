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

