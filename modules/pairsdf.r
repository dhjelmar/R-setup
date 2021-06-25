pairsdf <- function(df, var, hist='yes', color='black') {
    ## usage: pairsdf(mtcars, c('mpg', 'cyl', 'disp'))
    ##        pairsdf(iris, hist='yes', col=iris$Species)
    ##        pairsdf(iris, c('Sepal.Length', 'Petal.Width', 'Petal.Length'),
    ##                hist='yes', col=iris$Species)
    dfcompare <- subset(df,select=var)
    if (hist == 'no') {
        pairs(dfcompare, col=color)
    } else {
        ## add histograms, fits, and correlation coefficients
        if (color == 'black') {
            ## default is black
            ## drop the color so histogram is not solid black
            pairs(dfcompare,
                  lower.panel=panel.cor,
                  upper.panel=panel.smooth, 
                  diag.panel=panel.hist)
        } else {
            pairs(dfcompare,
                  lower.panel=panel.cor,
                  upper.panel=panel.smooth, 
                  diag.panel=panel.hist,
                  col=color)
        }
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

## pairs(iris,
##       lower.panel=panel.cor,
##       upper.panel=panel.smooth, 
##       diag.panel=panel.hist)
## pairs(iris,
##       lower.panel=panel.cor,
##       upper.panel=panel.smooth, 
##       diag.panel=panel.hist,
##       col=iris$Species)
