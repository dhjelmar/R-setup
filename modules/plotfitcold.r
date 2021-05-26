plotfitcold <- function(df,xx,yy,byvar,xlimspec=NULL,ylimspec=NULL,bg="grey90",outputfile=NULL) {
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

    ## debug
    ##xlabel <- deparse("gnom")
    ##ylabel <- deparse("average_g_ts")
    ##bylabel <- deparse("TS_ID")
    ##xx    <- eval(substitute(gnom),df)    ## need to recognize name passed into function as xx
    ##yy    <- eval(substitute(average_g_ts),df)    
    ##byvar <- eval(substitute(TS_ID),df)

    ## create dataframe with parameters to be plotted
    newdf <- data.frame(xx,yy,byvar)

    ## perform regression
    fit       <- lm(yy~xx,data=newdf)
    intercept <- fit$coefficients[[1]]
    slope     <- fit$coefficients[[2]]
    eq = paste0("y = ", signif(intercept,4), "* x + ", signif(slope,4))

    ## determine x and y range for plot
    xmin <- min(newdf$xx,na.rm=TRUE)
    xmax <- max(newdf$xx,na.rm=TRUE)
    xmax_plot <- xmax + (xmax-xmin)*0.2  # did this to make room for a legend
    ymin <- min(newdf$yy,na.rm=TRUE)
    ymax <- max(newdf$yy,na.rm=TRUE)
    if (missing(xlimspec)) { xlimspec <- c(xmin,xmax_plot) }
    if (missing(ylimspec)) { ylimspec <- c(ymin,ymax) }
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
    plot(newdf$xx,newdf$yy,xlim=xlimspec,ylim=ylimspec,
         xlab=xlabel,ylab=ylabel,main=eq,
         col=newdf$color)

    ## add fit to plot
    new.xx <- seq(min(newdf$xx,na.rm=TRUE),max(newdf$xx,na.rm=TRUE),len=100)
    fit    <- lm(yy~xx,data=newdf)
    pred   <- predict(fit, new=data.frame(xx=new.xx), interval="conf", level=1-0.05)
    lines(new.xx,pred[,"fit"],lwd=2)
    lines(new.xx,pred[,"lwr"],lty=3)
    lines(new.xx,pred[,"upr"],lty=3)

    ## add legend
    legendcolor <- unique(newdf$color)
    legendlabel <- unique(newdf$byvar)
    ##legend("right",title=bylabel,col=pal(ncol), pch=19,legend=round(legendlabel))
    ##legend("bottomright",title=bylabel,col=newdf$byvar, pch=19,legend=legendlabel)
    legend("bottomright", title=bylabel, col=legendcolor, pch=19, legend=legendlabel)

    ## four output file
    if (!missing(outputfile)) dev.off()

    ## print results of fit
    estbound(fit)

    ## return variables for debugging
    return(list( newdf=newdf, cols=cols, legendcolor=legendcolor, legendlabel=legendlabel, fit=fit))
}
## example: plotfitcold(df,gnom,average_g_ts,TS_ID,bg="white")
##          plotfitcold(mtcars, disp, mpg, byvar=cyl)

