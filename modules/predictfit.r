predictfit <- function(df, fit) {
    dffit <- data.frame(predict(fit, df, interval="conf", level=1-.05))
    names(dffit) <- c("fit","lwr","upr")
    dffit <- cbind(df, dffit)
}

#-----------------------------

