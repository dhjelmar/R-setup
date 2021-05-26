## This is a work in progress

plotres <- function(fit, df, x) {
    ## function to plot residuals
    fitsum <- summary(fit)
    residuals <- fitsum$residuals
    dffit <- cbind(df, dffit, residuals)
    #plotfit(dffit, logy, fit)
    #abline(a=0,b=1,col="red")   # y = a + b x
#
#    plotfit(dffit, fit, residuals)
#
#    plot(dffit$fit, dffit$residuals)
#    x <- dffit$fit
#    y <- dffit$residuals
#    smooth <- predict( lm(y ~ x + I(x^2) + I(x^3)) )
#    smooth <- data.frame(x,y,smooth)
#    smooth <- smooth[order(smooth$x),]
#    lines(smooth$x,smooth$smooth)
}


