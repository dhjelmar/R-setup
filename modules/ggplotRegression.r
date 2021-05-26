ggplotRegression <- function (fit) {
  require(ggplot2)
  titlesigfig <- 5
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, titlesigfig),
                       "Intercept =",signif(fit$coef[[1]],titlesigfig),
                       " Slope =",signif(fit$coef[[2]], titlesigfig),
                       " P =",signif(summary(fit)$coef[2,4], titlesigfig)))
}
#fit <- lm(hwy ~ displ, data=mpg)
#ggplotRegression(fit)

