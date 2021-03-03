qqplotqual <- function(x) {
    ## creates side by side, normal and Weibull qq plots
    par(mfrow=c(1,2))
    qqPlot(x,"normal", col='black')
    qqPlot(x,"Weibull", col='black')
}
## qqplotqual(mtcars$mpg)
