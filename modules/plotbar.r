plotbar <- function(seed, xnum, xcennum) {
  
  ## generate data
  set.seed(seed)
  x <- rnorm(xnum)
  x <- sort(x)
  xid <- round(runif(xcennum, 1, xnum))
  df1 <- data.frame(x.low=x, x.high=x, type='known')
  x.low  <- x[xid] * 0.50
  x.high <- x[xid] * 1.50
  df2 <- data.frame(x.low=x.low, x.high=x.high, type='censored')
  df <- rbind(df1, df2)
  df$x.med <- (df$x.low + df$x.high)/2
  
  ## create plot
  df <- df[order(df$x.med),]
  plot(df$x.med, ylim=range(df$x.low, df$x.high))
  lines(df$x.low)
  lines(df$x.high)
  return(df)
}

  plotspace(3,1)
  out <- plotbar(1, 10, 3)
  out <- plotbar(1, 100, 3)
  out <- plotbar(1, 1000, 3)


## Example of how to plot a dataset with error bars
plotbar_run <- function() {
  out[out$type=="censored",]
}
