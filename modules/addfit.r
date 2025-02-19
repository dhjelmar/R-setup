addfit <- function(xx, yy, col='black', interval='conf', alpha=0.05, sided=2, lty=1, pch=1, vlines=NULL, addpoints=TRUE) {
    ## lty only used for interval = 'mean' and 'line'
    
    if (addpoints) {
        ## plot points
        points(xx, yy, pch=pch, col=col)
    }
    
    ## put points in dataframe and remove any pairs that do not have values
    df <- na.omit( data.frame(xx,yy) )

    ## write points back out to xx and yy
    xx <- df$xx
    yy <- df$yy

    ## perform fit
    fit <- lm(yy~xx)
    intercept <- fit$coefficients[[1]]
    slope     <- fit$coefficients[[2]]

    ## add vertical lines to plot
    if ( !is.nothing(vlines) ) abline(v=vlines, lty=3, col=col)

    ## add prediction to plot
    if ( !is.nothing(vlines) ) {
        ## restrict fit to within vlines
        new.xx <- seq(min(vlines, na.rm=TRUE), max(vlines, na.rm=TRUE), len=100)
    } else {
        ## extend data to range of data
        new.xx <- seq(min(xx, na.rm=TRUE), max(xx, na.rm=TRUE), len=100)
    }
    if (interval == 'mean') {
        lines(new.xx, as.numeric(slope) * as.numeric(new.xx) + as.numeric(intercept), lty=lty, lwd=2, col=col)
        pred <- NA
    } else if (interval == 'line') {
        newdf <- data.frame(xx,yy)
        newdf <- newdf[order(newdf$xx),]
        xx <- newdf$xx
        yy <- newdf$yy
        lines(xx, yy, lty=lty, lwd=2, col=col)
        pred <- NA
    } else if (interval == 'conf' | interval == 'pred') {
        pred <- predict(fit, new=data.frame(xx=new.xx), interval=interval, level=1-alpha/sided)
        pred <- cbind(new.xx, pred)
        lines(new.xx, pred[,"fit"], lwd=2, col=col)
        lines(new.xx, pred[,"lwr"], lty=3, col=col)
        lines(new.xx, pred[,"upr"], lty=3, col=col)
    } else {
        ## add no line
        pred <- NA
    }
        
    ## print results of fit
    ## estbound(fit)

    ## calculate rise of fit over range bounds if supplied,
    ## otherwise calculate rise over range of data
    if ( !is.nothing(vlines) ) {
        rise <- (vlines[2] - vlines[1]) * slope
    } else {
        rise <- (max(xx) - min(xx)) * slope
    }        
    equation = paste0("y = ", signif(intercept,4), " + ", signif(slope,4), " * x; ",
                      expression(Delta), " = ", signif(rise,4))
    return(list(equation=equation, slope=slope, intercept=intercept, rise=rise, pred=pred))

}
