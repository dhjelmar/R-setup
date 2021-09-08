plotfitcol <- function(df,xx,yy,byvar,ncol=2,xlimspec=NULL,ylimspec=NULL,bg="grey90",outputfile=NULL) {
    ## create scatter plot with linear regression line and 95% confidence lines
    ## xx and yy are variables within dataframe df
    ## color datapoints by continuous variable byvar
    ## required inputs: xx, yy, df, byvar
    ## example usage: plotfitcol(df,hsub_in,power_ratio_fb,kin_fwd)                        # to write to screen
    ##                plotfitcol(df,hsub_in,power_ratio_fb,kin_fwd,outputfile="trend.jpg") # to write to jpg file
    ##                plotfitcol(df,hsub_in,power_ratio_fb,kin_fwd,ncol=4)                 # 4 colors for byvar

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
    if (typeof(byvar) == 'character' & length(byvar)==1)  {
        bycol <- which(grepl(byvar, names(df)))  
        byvar <- df[, bycol]
    }
    
    ## put xx, yy, and byvar into dataframe
    newdf <- data.frame(xx,yy,byvar)
    names(newdf) <- c('xx', 'yy', 'byvar')

    ## debug
    ##ncol <- 2
    ##xlabel <- deparse("psyst")
    ##ylabel <- deparse("power_ratio")
    ##bylabel <- deparse("gin")
    ##xx    <- eval(substitute(psyst),df)    ## need to recognize name passed into function as xx
    ##yy    <- eval(substitute(power_ratio),df)    
    ##byvar <- eval(substitute(gin),df)
    ##newdf <- data.frame(xx,yy,byvar)

    ## perform regression
    fit       <- lm(yy~xx,data=newdf)
    intercept <- fit$coefficients[[1]]
    slope     <- fit$coefficients[[2]]
    eq = paste0("y = ", signif(slope,4), " * x + ", signif(intercept,4))

    ## set min and max limits for plot ignoring NaN
    xmin <- min(newdf$xx,na.rm=TRUE)
    xmax <- max(newdf$xx,na.rm=TRUE)
    xmax_plot <- xmax + (xmax-xmin)*0.2  # did this to make room for a legend
    ymin <- min(newdf$yy,na.rm=TRUE)
    ymax <- max(newdf$yy,na.rm=TRUE)
    if (missing(xlimspec)) { xlimspec <- c(xmin,xmax_plot) }
    if (missing(ylimspec)) { ylimspec <- c(ymin,ymax) }
    if (!missing(outputfile)) jpeg(filename=outputfile)

    ## set color pallet for plot byvar
    ##cols <- c("blue","red")                         ## define colors
    ##cols <- brewer.pal(ncol,"Spectral")             # Spectral color pallette has ncol=11 colors
    ##cols <- c("black","blue","green","grey","orange","red")
    ##cols <- c("black","blue","darkorchid2","darkturquoise","darkgreen","green","deeppink","magenta","red")
    cols <- c("black","darkturquoise","pink","red")
    pal  <- colorRampPalette(cols)                  # define color pallette

    ## changes background of plot to specified color
    par(bg=bg)

    ## plot data points
    plot(newdf$xx,newdf$yy,xlim=xlimspec,ylim=ylimspec,
         xlab=xlabel,ylab=ylabel, main=eq,
         col=pal(ncol)[as.numeric(cut(newdf$byvar,breaks=ncol))])
    grid(col='gray70')

    ## add fit lines
    new.xx <- seq(min(newdf$xx,na.rm=TRUE),max(newdf$xx,na.rm=TRUE),len=100)
    pred   <- predict(fit, new=data.frame(xx=new.xx), interval="conf", level=1-0.05)
    lines(new.xx,pred[,"fit"],lwd=2)
    lines(new.xx,pred[,"lwr"],lty=3)
    lines(new.xx,pred[,"upr"],lty=3)

    ## add legend
    low  <- min(newdf$byvar,na.rm=TRUE)
    high <- max(newdf$byvar,na.rm=TRUE)
    legendlabel <- seq( low, high, (high-low)/(ncol-1) )
    ##legend("right",title=bylabel,col=pal(ncol), pch=19,legend=round(legendlabel))
    legend("right",title=bylabel,col=pal(ncol), pch=19,legend=signif(legendlabel,digits=4))

    if (!missing(outputfile)) dev.off()
    estbound(fit)
}
## examples
##plotfitcol(df,kin_fwd,power_ratio_fb,byvar=kin_fwd,ncol=21)
##plotfitcol(df,gin,    power_ratio_fb,byvar=kin_fwd,ncol=4,outputfile="junk2.jpg")
##plotfitcol(df,gin,    power_ratio_fb,      kin_fwd)
##plotfitcol(df,psyst, power_ratio,    gin,  ylimspec=c(0.7,1.1))
##plotfitcol(df,psyst, power_ratio , gin, ncol=10, ylim=c(0.7,1.05))

## plotfitcol(mtcars, disp, mpg, byvar=cyl)
