plotfit <- function(df, xx, yy, xlimspec=NULL, ylimspec=NULL, vlines=NULL,
                    main = 'equation', xlabel = NULL, ylabel = NULL,
                    interval='conf', alpha=0.05, sided=2, bg="grey90",
                    suppress='no') {
    ## create scatter plot with linear regression line and 95% confidence lines
    ## xx and yy are variables within dataframe df
    ## xlimspec and ylimspec can be used to force the range of the plot
    ## vlines will plot exactly two vertical lines at specified x locations
    ## main = 'equation'    # default is to print the fitted equation as the title
    ##      = 'your title'  # can specify title instead
    ## confidence or prediction limits
    ##     interval = 'conf' # confidence limits
    ##              = 'pred' # prediction limits
    ##              = 'none' # no limits added to plot
    ##     alpha    = 1 - confidence level
    ##     sided    = 2      # 2-sided interval
    ##              = 1      # 1-sided interval (currently plots upper and lower)
    ## bg controls the background color of the plot
    
    ## example usage: plotfit(df,hsub_in,power_ratio,ylimspec=c(0,1.1))
    ##                plotfit(mtcars, disp, mpg)
    ##                plotfit(mtcars, disp, mpg, vlines=c(100, 400))
    ##                plotfit(mtcars, disp, mpg, vlines=c(100, 400), interval='none')

    ## extract the names of xx and yy for plot axis labels if not defined
    if (missing(xlabel)) xlabel <- deparse(substitute(xx))
    if (missing(ylabel)) ylabel <- deparse(substitute(yy))

    ## extract xx and yy from df
    xx <- eval(substitute(xx),df)   # need to recognize name passed into function as xx
    yy <- eval(substitute(yy),df)   # need to recognize name passed into function as yy
    
    ## if xx and yy were passed in without quotes, xx and yy will be vectors to be plotted
    ## if xx and yy were passed in with quotes, xx and yy will be name of vector to be plotted
    if (typeof(xx) == 'character') {
        xxcol <- which(grepl(xx, names(df)))  
        xx    <- df[, xxcol]
    }
    if (typeof(yy) == 'character')  {
        yycol <- which(grepl(yy, names(df)))  
        yy    <- df[, yycol]
    }
    
    ## put xx and yy into dataframe
    newdf <- data.frame(xx,yy)

    ## set min and max for plot if not specified
    xmin <- min(newdf$xx, xlimspec, na.rm=TRUE)
    xmax <- max(newdf$xx, xlimspec, na.rm=TRUE)
    ## xmax_plot <- xmax + (xmax-xmin)*1.2  # did this to make room for a legend
    ymin <- min(newdf$yy,na.rm=TRUE)
    ymax <- max(newdf$yy,na.rm=TRUE)
    if (missing(xlimspec)) { xlimspec <- c(xmin,xmax) }
    if (missing(ylimspec)) { ylimspec <- c(ymin,ymax) }

    ## change background of plot to specified color
    par(bg=bg)  

    ## perform regrssion
    new.xx <- seq(min(newdf$xx, vlines, na.rm=TRUE),max(newdf$xx, vlines, na.rm=TRUE), len=100)
    fit <- lm(yy~xx,data=newdf)
    ## estbound(fit)
    pred   <- predict(fit, new=data.frame(xx=new.xx), interval=interval, level=1-0.05/sided)
    intercept <- fit$coefficients[[1]]
    slope <- fit$coefficients[[2]]
    ## calculate rise of fit over range bounds if supplied, otherwise calculate rise over range of data
    if ( !is.null(vlines) ) {
        rise <- (vlines[2] - vlines[1]) * slope
    } else {
        rise <- (xmax - xmin) * slope
    }        

    ## plot points and fit
    if (main == 'equation') main = paste0("y = ", signif(slope,4), "* x + ", signif(intercept,4),
                                          ", ", expression(Delta), " = ", signif(rise,4))
    plot(newdf$xx,newdf$yy,xlim=xlimspec,ylim=ylimspec,
         xlab=xlabel,ylab=ylabel,
         main=main)
    grid(col='gray70')
    lines(new.xx,pred[,"fit"],lwd=1, col='red')

    ## if request confidence bounds, add to plot
    if (interval != 'none') {
        lines(new.xx,pred[,"lwr"],lty=3, col='red')
        lines(new.xx,pred[,"upr"],lty=3, col='red')
    }

    ## if request upper/lower range bounds, add to plot
    if ( !is.null(vlines) ) {
        abline(v=vlines[1], col='black', lty='dashed')
        abline(v=vlines[2], col='black', lty='dashed')
    }        

    if (suppress != 'yes')
        return( list(intercept = intercept, slope = slope, rise = rise, pred = as_tibble(pred)) )
    
}
testplots <- function() {
    source('/home/dlhjel/GitHub_repos/R-setup/setup.r')
    plotfit(mtcars, cyl, 'mpg')
    plotfit(mtcars, cyl, 'mpg', vlines=c(2,9), xlimspec=c(0, 10))
    plotfitcol(mtcars, cyl, 'mpg', byvar='cyl', ncol=3)
    plotfitcold(mtcars, cyl, 'mpg', byvar=cyl)
}
