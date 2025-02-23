standardize <- function(df, meanvec=NULL, sdvec=NULL) {
    ## standardizes data in dataframe (i.e., x --> (x-mean(x))/sd(x)
    ## df          = dataframe of parameters to be scaled
    ## meanvec     = vector of mean values for centering
    ##               default uses means of df for centering
    ## sdvec       = vector of standard deviations
    ##               default uses sd(x)
    
    ## determine mean and sd if not supplied
    if (missing(meanvec)) meanvec <- apply(df, 2, mean)  # use mean of df
    if (missing(sdvec))   sdvec   <- apply(df, 2, sd)  # use sd of df
    
                                        # standardize data
    dfcenter <- sweep(df,       2, meanvec,    "-")  
    dfs      <- sweep(dfcenter, 2, sdvec  ,    "/")
    
    ## return parameters
    return( list(dfs=dfs, meanvec=meanvec, sdvec=sdvec) )
}

standardize_rev <- function(dfs, meanvec, sdvec, coef_s=NULL, y=NULL) {
    ## reverses standardization of dataframe
    ##
    ## inputs:
    ##     dfs     = standardized dataframe
    ##     meanvec = vector of mean values used to standardize dataframe
    ##     sdvec   = vector of standard deviations used to standardize dataframe
    ##     coef_s  = vector of coefficients from regression on standardized dataframe
    ##     y       = name of indepenent variable in coef_s (required if coef_s != NULL)
    ##     fit     = how to calculate fit results sing non-standardized parameters
    ##     fit_s   = vector with fit results calculated using standardized parameters
    ## note:  function assumes dependent and independent model parameters were scaled
    df <- sweep(dfs, 2, sdvec  , "*")
    df <- sweep(df , 2, meanvec, "+")
    if (missing(coef_s)) {
        return(df)
    } else {
        ## convert standardized coefficients to non-standardized and return fit values
        if (missing(y)){return('ERROR: In standardize_rev(), coef_s supplied requires y to also be supplied')}
        cs <- coef_s
        bs <- cs[1]              # 1st coefficient is always the intercept
        cs <- cs[2:length(cs)]
        ## identify vector location with independent variable
        yloc <- which(names(meanvec) == y)
        xloc <- -yloc
        ## non-standardized model intercept
        b <- bs - cs %*% as.matrix(meanvec[xloc] / sdvec[xloc]) * sdvec[yloc] + meanvec[yloc]
        b <- b[[1]]
        ## unscale model coefficients
        coef <- cs[xloc] * sdvec[yloc] / sdvec[xloc]
        ## determine fit results
        fit_s <- (bs + as.matrix(dfs[xloc]) %*% as.matrix(cs)) * sdvec[yloc] + meanvec[yloc]
        fit <- 'fit <- b + as.matrix(non-standardized x-data) %*% as.matrix(coef)'
        return(list(b=b, coef=coef, fit=fit, fit_s=fit_s))
    }
}

standardize_example <- function() {
    ## define data
    df <- mtcars[,c('cyl', 'disp', 'hp', 'wt', 'mpg')]
    ## standardize data
    dfs.out <- standardize(df)
    dfs     <- dfs.out$dfs
    meanvec <- dfs.out$meanvec
    sdvec <- dfs.out$sdvec
    summary(dfs)
    ## fit data
    model <- lm(mpg ~ cyl + disp + hp + wt, data=dfs)
    coef_s <- coef(model)
    ## reverse standardization
    out <- standardize_rev(dfs, meanvec, sdvec, coef_s, y='mpg')
    cat('coef_s: ', coef_s, '\n')
    cat('b     : ', out$b, '\n')
    cat('coef  : ', out$coef, '\n')
    ## fit
    y <- 'mpg'
    yloc <- which(names(meanvec) == y)
    df$fit <- out$b + as.matrix(df[-yloc]) %*% as.matrix(out$coef)
    df$fit_s <- out$fit_s
    ## plot residuals
    plotspace(2,2)
    plotfit(df$cyl,  df$fit_s - df$mpg, main='standardized') 
    plotfit(df$cyl,  df$fit - df$mpg, main='non-standardized') 
    plotfit(df$disp, df$fit_s - df$mpg, main='standardized') 
    plotfit(df$disp, df$fit - df$mpg, main='non-standardized') 
    plotfit(df$hp,   df$fit_s - df$mpg, main='standardized') 
    plotfit(df$hp,   df$fit - df$mpg, main='non-standardized') 
    plotfit(df$wt,   df$fit_s - df$mpg, main='standardized') 
    plotfit(df$wt,   df$fit - df$mpg, main='non-standardized') 
    return(df)
}
#out <- standardize_example()
#head(out)
