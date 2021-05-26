## This is a work in progress

plotxy <- function(df,xx,yy,model=1,xlimspec=NULL,ylimspec=NULL,bg="grey90") {
    ## create scatter plot with linear regression line and 95% confidence lines
    ## xx and yy are variables within dataframe df
    ## example usage: plotfit(df,hsub_in,power_ratio,ylimspec=c(0,1.1))
    xlabel <- deparse(substitute(xx))
    ylabel <- deparse(substitute(yy))
    xx <- eval(substitute(xx),df)   # need to recognize name passed into function as xx
    yy <- eval(substitute(yy),df)   # need to recognize name passed into function as yy
    newdf <- data.frame(xx,yy)
    xmin <- min(newdf$xx,na.rm=TRUE)
    xmax <- max(newdf$xx,na.rm=TRUE)
    ## xmax_plot <- xmax + (xmax-xmin)*1.2  # did this to make room for a legend
    ymin <- min(newdf$yy,na.rm=TRUE)
    ymax <- max(newdf$yy,na.rm=TRUE)
    if (missing(xlimspec)) { xlimspec <- c(xmin,xmax) }
    if (missing(ylimspec)) { ylimspec <- c(ymin,ymax) }
    par(bg=bg)  # changes background of plot to specified color
    plot(newdf$xx,newdf$yy,xlim=xlimspec,ylim=ylimspec,xlab=xlabel,ylab=ylabel)
    ##plot(newdf$power_ratio~newdf$g_ave,type="n")
    ##wx    <- par("usr")[1:2]
    ##new.x <- seq(wx[1],wx[2],len=length(newdf$power_ratio))
    ##new.x <- seq(min(log10(newdf$g_ave)),max(log10(newdf$g_ave)),len=100)
    new.xx <- seq(min(newdf$xx,na.rm=TRUE),max(newdf$xx,na.rm=TRUE),len=100)
    if        (model==0) {
        fit <- lm(yy~1,                      data=newdf)
    } else if (model == 1) {
        fit <- lm(yy~xx,                     data=newdf)
    } else if (model == 2) {
        fit <- lm(yy~xx + I(xx^2),           data=newdf)
    } else {
        fit <- lm(yy~xx + I(xx^2) + I(xx^3), data=newdf)
    }
    pred   <- predict(fit, new=data.frame(xx=new.xx), interval="conf", level=1-0.05)
    lines(new.xx,pred[,"fit"],lwd=2)
    if (model==1) {
        lines(new.xx,pred[,"lwr"],lty=3)
        lines(new.xx,pred[,"upr"],lty=3)
    }
    estbound(fit)
}

