plotfit2 <- function(df1, df2, xx, yy, byvar,
                     xrange=NULL, yrange=NULL,
                     vlines1=NULL, vlines2=NULL,
                     color=c('black', 'red'), main=NULL) {

    ##-----------------------------------------------------------------------------
    ## extract column locations and values from dataframe
    ## df1
    xxcol <- which(grepl(xx, names(df1)))  ## xx column
    xx1   <- df1[, xxcol]                  ## xx values
    yycol <- which(grepl(yy, names(df1)))  
    yy1   <- df1[, yycol]
    bycol <- which(grepl(byvar, names(df1)))  
    by1   <- df1[, bycol]
    ## df2
    xxcol <- which(grepl(xx, names(df2)))  
    xx2   <- df2[, xxcol]
    yycol <- which(grepl(yy, names(df2)))  
    yy2   <- df2[, yycol]
    bycol <- which(grepl(byvar, names(df2)))  
    by2   <- df2[, bycol]
    
    ##-----------------------------------------------------------------------------
    ## define plot parameters
    xmax   <- max(xrange, xx1, xx2, vlines1, vlines2)
    xmin   <- min(xrange, xx1, xx2, vlines1, vlines2)
    yrange <- range(yrange, yy1, yy2)
    ## make extra room for legend
    xmax_plot <- xmax + (xmax-xmin)*0.2
    xrange <- range(xmin, xmax_plot)
    
    ##-----------------------------------------------------------------------------
    ## create empty plot
    df <- rbind(df1, df2)
    plot(xx1, yy1, type='n',
         xlim=xrange, ylim=yrange,
         xlab=xx,     ylab=yy,
         main=main)

    ##-----------------------------------------------------------------------------
    ## add nofit points
    ## not programmed

    ##-----------------------------------------------------------------------------
    ## add points, fits and vlines
    ## fit is over range of data
    ## lines extended to max(vlines, data)
    eq1 <- addfit(xx1, yy1, col=color[1], vlines=vlines1)
    eq2 <- addfit(xx2, yy2, col=color[2], vlines=vlines2)

    ##-----------------------------------------------------------------------------
    ## add the 2 fit equations under main title
    subtitle <- list(eq1, eq2)
    mtext(subtitle, side=3, line=c(0.75, 0), cex=.75, col=color)

    ##-----------------------------------------------------------------------------
    ## add legend
    ## determine legend location by slope of all data
    fitall       <- lm(yy1~xx1, data=df)
    fitall_slope <- fitall$coefficients[[2]]
    if (fitall_slope > 0) {
        legendloc <- 'bottomright'
    } else {
        legendloc <- 'topright'
    }
    legend(legendloc, title=byvar, col=color, legend=unique(df[[bycol]]), pch=1)
    
    ##-----------------------------------------------------------------------------
    ## return fitted equations
    return(list(eq1=eq1, eq2=eq2))
}

## df1 <- mtcars
## df1$type <- 'type1'
## df2     <- df1
## df2$mpg <- df1$mpg * 1.1 + rnorm(1, 10, 1)
## df2$type <- 'type2'

## ## separate dataframes
## df1 <- df[df$type == 'type1',]
## df2 <- df[df$type == 'type2',]

## plotfit2(df1, df2, 'mpg', 'disp', 'type')
## plotfit2(df1, df2, 'mpg', 'disp', 'type',
##          xrange=c(0, 50), yrange=c(0, 500),
##          vlines1=c(10,35), vlines2=c(20,47),
##          main='my title')
