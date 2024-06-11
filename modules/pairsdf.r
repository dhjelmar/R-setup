pairsdf <- function(df, var=names(df), upr='panel.line', diag='panel.hist', lwr='panel.cor',
                    col.points='black', col.line='red', main=NULL) {
    ## usage: pairsdf(mtcars, c('mpg', 'cyl', 'disp'))
    ##        pairsdf(iris, col.points=iris$Species, col.line='blue')
    ## options: upr/lwr = for upper (or lower) panels:
    ##                      - 'panel.line' (default) plots linear fit against data
    ##                      - 'panel.smooth' plots standard R 'panel.smooth' function against data
    ##                      - 'points' just plots the data points (this is the standard pairs() default)
    ##          diag = 'panel.hist' (default) plots histogram on each diagonal panel
    ##               = NULL only plots the variable name on each diagonal (this is the standard pairs() default)
    dfcompare <- subset(df,select=var)
    pairs(dfcompare,
          upper.panel=upr, 
          diag.panel=diag,
          lower.panel=lwr,
          col.points=col.points,
          col.line=col.line)
}

##----------------------------------------------------------------------------

panel.line <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.points = 'black', col.line = 'red', span = 2/3, iter = 3, ...) {
    points(x, y, pch = pch, col=col.points, col.line = 'red', bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    ## the following commented out is what panel.smooth does
    ##if (any(ok))
    ##    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
    ##        col = col.smooth, ...)
    if (any(ok)) {
        ## perform fit
        xx <- x[ok]
        yy <- y[ok]
        fit <- lm(yy~xx)
        intercept <- fit$coefficients[[1]]
        slope     <- fit$coefficients[[2]]
        lines(xx, as.numeric(slope) * as.numeric(xx) + as.numeric(intercept), lty=1, lwd=2, col=col.line)
    }
}

##----------------------------------------------------------------------------

panel.hist <- function(x, col='black', ...) {
    ## not sure how to get histogram for each color variable
    ## note can use python seaborn for this with R reticulate library
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, ...)
}

##----------------------------------------------------------------------------

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairsdf_test <- function() {
  pairs(iris,
        lower.panel=panel.cor,
        upper.panel=panel.line, 
        diag.panel=panel.hist)
  pairs(iris,
        lower.panel=panel.cor,
        upper.panel=panel.line, 
        diag.panel=panel.hist,
        col=iris$Species)
  pairsdf(iris)
  pairsdf(iris)
  pairsdf(iris, col.points=iris$Species, col.line='blue')
}
