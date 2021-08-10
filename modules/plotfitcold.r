plotfitcold <- function(df,xx,yy,byvar,xlimspec=NULL,ylimspec=NULL, vlines=NULL,
                        main = 'equation', xlabel = NULL, ylabel = NULL,
                        interval='conf', alpha=0.05, sided=2, bg="grey90",
                        suppress='no', outputfile=NULL) {
    ## create scatter plot with linear regression line and 95% confidence lines
    ## xx and yy are variables within dataframe df
    ## color datapoints by DISCRETE variable byvar  <-- This is the distinction from plotfitcol
    ## required inputs: df, xx, yy, byvar
    ## example usage: plotfitcold(df,hsub_in,power_ratio_fb,TS_ID)                        # to write to screen
    ##                plotfitcold(df,hsub_in,power_ratio_fb,TS_ID,outputfile="trend.jpg") # to write to jpg file
    ##                plotfitcold(df,hsub_in,power_ratio_fb,TS_ID,ncol=4)                 # 4 colors for byvar

    ## uncomment following breakpoit to enter debugger when source in rstudio
    ## browser()
    ## source main program then use F10 to step through after break point
    
    xlabel  <- deparse(substitute(xx))
    ylabel  <- deparse(substitute(yy))
    bylabel <- deparse(substitute(byvar))
    xx    <- eval(substitute(xx),df)    # need to recognize name passed into function as xx
    yy    <- eval(substitute(yy),df)    
    byvar <- eval(substitute(byvar),df) 

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
    if (typeof(byvar) == 'character')  {
        bycol <- which(grepl(byvar, names(df)))  
        byvar <- df[, bycol]
    }
    
    ## put xx, yy, and byvar into dataframe
    newdf <- data.frame(xx,yy,byvar)

    ## determine x and y range for plot
    xmin <- min(newdf$xx, xlimspec, na.rm=TRUE)
    xmax <- max(newdf$xx, xlimspec, na.rm=TRUE)
    xmax_plot <- xmax + (xmax-xmin)*0.2  # did this to make room for a legend
    ymin <- min(newdf$yy,na.rm=TRUE)
    ymax <- max(newdf$yy,na.rm=TRUE)
    if (missing(xlimspec)) { xlimspec <- c(xmin,xmax_plot) }
    if (missing(ylimspec)) { ylimspec <- c(ymin,ymax) }

    ## perform regression
    fit       <- lm(yy~xx,data=newdf)
    intercept <- fit$coefficients[[1]]
    slope     <- fit$coefficients[[2]]
    ## calculate rise of fit over range bounds if supplied, otherwise calculate rise over range of data
    if ( !is.null(vlines) ) {
        rise <- (vlines[2] - vlines[1]) * slope
    } else {
        rise <- (xmax - xmin) * slope
    }        

    ## plot points and fit
    if (main == 'equation') main = paste0("y = ", signif(slope,4), "* x + ", signif(intercept,4),
                                          ", ", expression(Delta), " = ", signif(rise,4))

    ## setup for output jpeg file
    if (!missing(outputfile)) jpeg(filename=outputfile)

    ## add color definition for each point using byvar
    ##cols <- c("blue","red")                         # define colors
    ##cols <- brewer.pal(11,"Spectral")             # Spectral color pallette has ncol=11 colors
    ##cols <- c("black","blue","green","grey","orange","red")
    ##cols <- c("black","blue","darkorchid2","darkturquoise","darkgreen","green","deeppink","magenta","red")
    ##cols <- c("black","darkturquoise","pink","red")
    ##pal  <- colorRampPalette(cols)                  # define color pallette
    cols <- data.frame(byvar=unique(c(as.character(newdf$byvar))))
    ncol <- nrow(cols)
    pal <- rainbow(ncol)
    ## library(scales) # needed for show_col
    ## show_col(pal)
    cols$color <- pal
    newdf <- merge(newdf, cols, by="byvar")

    ## create plot of points
    par(bg=bg)  # changes background of plot to specified color
    plot(newdf$xx, newdf$yy, xlim=xlimspec, ylim=ylimspec,
         xlab=xlabel, ylab=ylabel, main=main,
         col=newdf$color)
    grid(col='gray70')

    ## add fit to plot
    new.xx <- seq(min(newdf$xx, vlines, na.rm=TRUE),max(newdf$xx, vlines, na.rm=TRUE), len=100)
    fit    <- lm(yy~xx,data=newdf)
    pred   <- predict(fit, new=data.frame(xx=new.xx), interval=interval, level=1-0.05/sided)
    lines(new.xx,pred[,"fit"],lwd=2)

    ## if request confidence bounds, add to plot
    if (interval != 'none') {
        lines(new.xx,pred[,"lwr"],lty=3)
        lines(new.xx,pred[,"upr"],lty=3)
    }

    ## if request upper/lower range bounds, add to plot
    if ( !is.null(vlines) ) {
        abline(v=vlines[1], col='black', lty='dashed')
        abline(v=vlines[2], col='black', lty='dashed')
    }        
    
    ## add legend
    legendcolor <- unique(newdf$color)
    legendlabel <- unique(newdf$byvar)
    ##legend("right",title=bylabel,col=pal(ncol), pch=19,legend=round(legendlabel))
    ##legend("bottomright",title=bylabel,col=newdf$byvar, pch=19,legend=legendlabel)
    legend("bottomright", title=bylabel, col=legendcolor, pch=19, legend=legendlabel)

    ## for output file
    if (!missing(outputfile)) dev.off()

    ## print results of fit
    estbound(fit)

    if (suppress == 'no')
        return( list(intercept = intercept, slope = slope, rise = rise, pred = as_tibble(pred)) )
}
## example: plotfitcold(df,gnom,average_g_ts,TS_ID,bg="white")
##          plotfitcold(mtcars, disp, mpg, byvar=cyl)

