plotbar <- function(y, ylow, yhigh, ylab='y', sorty=TRUE) {
  ## Given values y with low and high values representing uncertainty
  ## Plot the y values as open circles and add lines connecting low and high values

  ## put data into dataframe and sort if needed
  df <- data.frame(y=y, ylow=ylow, yhigh=yhigh)
  if (isTRUE(sorty)) df <- df[order(df$y),]
  df$index <- c(1:nrow(df))

  ## plot all points
  plot(df$index, df$y, ylim=range(ylow, yhigh), xlab='Index', ylab=ylab)
  
  ## plot known points
  known <- which(df$ylow==df$yhigh)
  points(df$index[known], df$y[known], pch=16)
  
  ## add "error bar" lines
  lines(df$index, df$ylow)
  lines(df$index, df$yhigh)
  
  ## add legend
  legend('topleft', 
         legend=c('Censored (lines for range)', 'Known'), 
         pch=c(1,16))
  
  return(df)
}

plotbar_run <- function() {
  plotspace(2,2)
  ylow  <- c(1,6,2,2.5,3.8,5)
  y     <- c(1,6,2,4,5,5)
  yhigh <- c(1,6,2,5,6,5)
  plotbar(y=y, ylow=ylow, yhigh=yhigh, sorty=NA)
  plotbar(y=y, ylow=ylow, yhigh=yhigh)
  plotbar(mtcars$mpg, mtcars$mpg*0.8, mtcars$mpg*1.2, 'mpg from mtcars database', sorty=NA)
  plotbar(mtcars$mpg, mtcars$mpg*0.8, mtcars$mpg*1.2, 'mpg from mtcars database')
}
