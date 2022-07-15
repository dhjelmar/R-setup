scaledf <- function(df, meanvec = NULL, scaledata = NULL, scalefactor = 2, zeromin = FALSE, suffix=NA) {
    ## centers and scales a dataframe
    ## df          = dataframe of parameters to be scaled
    ## meanvec     = vector of mean values for centering
    ##               default uses means of df for centering
    ## scaledata   = dataframe used to provide range of each parameter for scaling
    ##               default uses range of df for scaling
    ## scalefactor = number used to set range of scaled data (multiplies data by scalefactor / range of scaledata)
    ##               (2 results in -1 to +1 if centered and ranges of scaledata are equal distant from mean)
    ##               (1 results in  0 to +1 if zeroed   and ranges of scaledata are equal distant from mean)
    ## zeromin     = FALSE (default) results in center of scaled df at 0 if ranges of scaledata are equal distant from mean)
    ##             = TRUE results in min of scaled df at 0 if ranges of scaledata are equal distant from mean)
    ## suffix      = NA (default) does nothing
    ##             = character (e.g., "_s") adds character to each variable to identify it is scaled
    
    ## center data
    if (missing(meanvec))   meanvec <- apply(df, 2, mean)  # use mean of df  to center
    # subtract meanvec from each row in df; output is a dataframe
    dfcenter <- sweep(df,       2, meanvec,    "-")  
    
    ## scale centered data
    if (missing(scaledata)) scaledata <- df                # use range of df to scale
    ## calculate vector (1 value for each parameter in scaledata) = range of parameter / scalefactor
    scalerange <- (sapply(scaledata,max) - sapply(scaledata,min)) / scalefactor
    ## divide every row in dfcenter by scalerange
    dfs        <- sweep(dfcenter, 2, scalerange, "/")

    ## adjust location center if needed
    if (zeromin == TRUE) dfs <- sweep(dfs, 2, (scalefactor/2), "+")
    
    ## convert to dataframe and a suffix to column names
    if (!is.na(suffix)) {
        dfs <- data.frame(dfs)
        names(dfs) <- paste(names(dfs), "_s", sep = "")
    }
    
    ## return parameters
    return( list(dfs=dfs, meanvec=meanvec, scaledata=scaledata, scalefactor=scalefactor, zeromin=zeromin) )
}
# dfs.out <- scaledf(df)$dfs
# dfs     <- dfs.out$dfs
# summary(dfs)

