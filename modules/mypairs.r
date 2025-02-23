mypairs <- function(data) {
    
    # creates pairs plot with:
    #     - correlation coefficients in upper half 
    #     - liner fits in lower half  
  
    # drop non-numeric columns
    nums <- sapply(data, is.numeric)
    data <- data[,nums]
    
    pairs(data, 
          upper.panel = panel.cor,
          lower.panel = panel.fit,
          horOdd=TRUE,
          gap=0)
}

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.fit <- function(x,y){
  points(x,y)
  abline(lm(y~x), col='red')
}

# mypairs(mtcars[,c(1:6)])

