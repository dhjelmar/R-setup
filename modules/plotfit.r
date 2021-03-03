plotfit <- function(df, xx, yy, xlimspec=NULL, ylimspec=NULL, vlines=NULL, confbounds='yes', bg="grey90") {
    ## create scatter plot with linear regression line and 95% confidence lines
    ## xx and yy are variables within dataframe df
    ## xlimspec and ylimspec can be used to force the range of the plot
    ## vlines will plot exactly two vertical lines at specified x locations
    ## confbounds requests that 95% confidence bounds be applied to the plot
    ## bg controls the background color of the plot
    
    ## example usage: plotfit(df,hsub_in,power_ratio,ylimspec=c(0,1.1))
    ##                plotfit(mtcars, disp, mpg)
    ##                plotfit(mtcars, disp, mpg, vlines=c(100, 400))
    ##                plotfit(mtcars, disp, mpg, vlines=c(100, 400), confbounds='no')

    ## extract the names of xx and yy
    xlabel <- deparse(substitute(xx))
    ylabel <- deparse(substitute(yy))

    ## extract xx and yy from df and create dataframe
    xx <- eval(substitute(xx),df)   # need to recognize name passed into function as xx
    yy <- eval(substitute(yy),df)   # need to recognize name passed into function as yy
    newdf <- data.frame(xx,yy)

    ## set min and max for plot if not specified
    xmin <- min(newdf$xx,na.rm=TRUE)
    xmax <- max(newdf$xx,na.rm=TRUE)
    ## xmax_plot <- xmax + (xmax-xmin)*1.2  # did this to make room for a legend
    ymin <- min(newdf$yy,na.rm=TRUE)
    ymax <- max(newdf$yy,na.rm=TRUE)
    if (missing(xlimspec)) { xlimspec <- c(xmin,xmax) }
    if (missing(ylimspec)) { ylimspec <- c(ymin,ymax) }

    ## change background of plot to specified color
    par(bg=bg)  

    ## perform regrssion
    new.xx <- seq(min(newdf$xx,na.rm=TRUE),max(newdf$xx,na.rm=TRUE),len=100)
    fit <- lm(yy~xx,data=newdf)
    ## estbound(fit)
    pred   <- predict(fit, new=data.frame(xx=new.xx), interval="conf", level=1-0.05)
    intercept <- fit$coefficients[[1]]
    slope <- fit$coefficients[[2]]

    ## plot points and fit
    eq = paste0("y = ", signif(intercept,4), "* x + ", signif(slope,4))
    plot(newdf$xx,newdf$yy,xlim=xlimspec,ylim=ylimspec,
         xlab=xlabel,ylab=ylabel,
         main=eq)
    lines(new.xx,pred[,"fit"],lwd=2)

    ## if request confidence bounds, add to plot
    if (confbounds == 'yes') {
        lines(new.xx,pred[,"lwr"],lty=3)
        lines(new.xx,pred[,"upr"],lty=3)
    }

    ## if request upper/lower range bounds, add to plot
    ## return rise of fit over range bounds if supplied, otherwise return rise over range of data
    if ( !is.null(vlines) ) {
        abline(v=vlines[1], col='black', lty='dashed')
        abline(v=vlines[2], col='black', lty='dashed')
        rise <- (vlines[2] - vlines[1]) * slope
    } else {
        rise <- (xmax - xmin) * slope
    }        

    return( list(intercept = intercept, slope = slope, rise = rise) )
    
}
