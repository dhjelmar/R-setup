addfit <- function(xx, yy, col='black', interval='conf', alpha=0.05, sided=2, pch=1, vlines=NULL) {
    ## plot points
    points(xx, yy, pch=pch, col=col)

    ## perform fit
    fit <- lm(yy~xx)
    intercept <- fit$coefficients[[1]]
    slope     <- fit$coefficients[[2]]

    ## add vertical lines to plot
    if ( !is.null(vlines[1]) & !is.na(vlines[1]) ) abline(v=vlines, lty=3, col=col)

    ## add prediction to plot
    if ( !is.null(vlines[1]) & !is.na(vlines[1]) ) {
        ## restrict fit to within vlines
        new.xx <- seq(min(vlines, na.rm=TRUE), max(vlines, na.rm=TRUE), len=100)
    } else {
        ## extend data to range of data
        new.xx <- seq(min(xx, na.rm=TRUE), max(xx, na.rm=TRUE), len=100)
    }
    pred   <- predict(fit, new=data.frame(xx=new.xx), interval=interval, level=1-alpha/sided)
    lines(new.xx, pred[,"fit"], lwd=2, col=col)
    lines(new.xx, pred[,"lwr"], lty=3, col=col)
    lines(new.xx, pred[,"upr"], lty=3, col=col)

    ## print results of fit
    ## estbound(fit)

    ## calculate rise of fit over range bounds if supplied,
    ## otherwise calculate rise over range of data
    if ( !is.null(vlines[1]) & !is.na(vlines[1]) ) {
        rise <- (vlines[2] - vlines[1]) * slope
    } else {
        rise <- (max(xx) - min(xx)) * slope
    }        
    equation = paste0("y = ", signif(intercept,4), " + ", signif(slope,4), " * x; ",
                      expression(Delta), " = ", signif(rise,4))
    return(list(equation=equation, slope=slope, intercept=intercept, rise=rise))

}
