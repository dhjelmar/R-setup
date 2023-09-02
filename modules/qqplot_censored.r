qqplot_censored <- function(x, xcen, type='j', plotit='q') {
  ## x    = 'Known'
  ## xcen = dataframe of x.low, x.high pairs
  ## type = 'n' to fit normal distribution
  ##        'w' to fit Weibull distribution
  ##        'j' to fit Johnson SU distribution
  ## plotit = 'x' plots x vs. index
  ##          'h' plots histogram of known x values and pdf
  ##          'e' plots ECDF
  ##          'q' (default) plots only the qq plot
  ##          can contain any combination of q, x, h and e (e.g., plotit='xq')
  
  df1 <- data.frame(x.low =      x    , x.high =      x     , type='Known')
  df2 <- data.frame(x.low = xcen$x.low, x.high = xcen$x.high, type='Censored')
  df <- rbind(df1,df2)
  df$x <- df$x.high
  
  if (plotit == 'xheq') plotspace(2,2)
  
  if (grepl('x', plotit)) {
    ############################################################
    # PLOT THE DATA TO SEE CENSORED AND KNOWN
    ############################################################  plotspace(3,1)
    ## sort from lowest to highest with censored values before known
    df <- df[order(df$type, decreasing = FALSE),]
    df <- df[order(df$x, decreasing = FALSE),]
    df$index <- c(1:nrow(df))
    plotfit(df$index, df$x, interval='noline')
    addfit(df[df$type == 'Known',]$index, df[df$type == 'Known',]$x, pch=16, interval='noline')
    legend('topleft', legend=c('Censored', 'Known'), pch=c(1,16))
  }
  
  
  ############################################################
  # FIT THE DATA TO DETERMINE THEORETICAL PDF
  ############################################################
  if (grepl('h', plotit)) {
    out <- hist_nwj(x, xcen, type=type, tolerance=FALSE)
  } else {
    out <- hist_nwj(x, xcen, type=type, tolerance=FALSE, plot=FALSE)
  }
  nparms <- out$nparms
  wparms <- out$wparms
  jparms <- out$jparms
  
  
  ############################################################
  # ESTMATE CDF (ECDF)
  ############################################################
  ## sort from highest to lowest with censored values before known
  df <- df[order(df$type, decreasing = TRUE),]
  df <- df[order(df$x, decreasing = TRUE),]
  row.names(df) <- c(1:nrow(df))
  df$i <- c(nrow(df):1)
  
  ## known data
  df$known <- ifelse(df$type == 'Known', df$x, NA)
  ## df$lag <- lag(df$x)
  df$lag <- c(NA, df[1:(nrow(df)-1),'x'])
  df$lag[1] <- 0
  
  ## index of known and unique
  ku <- which(df$type == 'Known' & df$x != df$lag)
  df$lag <- NULL
  
  ## zj = Values of known and unique
  df[ku,'zj'] <- df[ku,'x']
  
  ## rj = known at risk
  df[ku,'rj'] <- df[ku,'i']
  
  ## nj = number of known values at zj
  df['count'] <- with(df,ave(known, known, FUN=function(x) length(x)))
  df[ku, 'nj'] <- df[ku, 'count']
  df['count'] <- NULL
  df
  
  ## dataframe (dfk) of unique, known data (Zj)
  dfk <- df[!is.na(df$zj),]
  dfk$j <- c(1:nrow(dfk))
  dfk <- dfk[,c('j', 'zj', 'rj', 'nj')]
  
  ## estimate the CDF (Fhat or ECDF)
  dfk$fraction <- (dfk$rj - dfk$nj)/dfk$rj
  dfk$Fhat <- 1
  for (j in c(2:nrow(dfk))) {
    dfk[j, 'Fhat'] <- dfk[j-1, 'Fhat'] * dfk[j-1, 'fraction']
  }
  if (grepl('e', plotit)) {
    plotfit(dfk$zj, dfk$Fhat, main='Empirical CDF', interval='line')
  }
  
  
  if (grepl('q', plotit)) {
    ############################################################
    # THEORETICAL CDF
    ############################################################
    
    ## sort data from lowest to highest
    ## put left censored data before known of same value
    ## where left censored have x.low = 0
    df <- df[order(df$type, decreasing = FALSE),]
    df <- df[order(df$x, decreasing = FALSE),]
    df$i <- c(1:nrow(df))
    
    ## theoretical CDF
    if (type == 'n') {
        dfk$cdf <- pnorm(dfk$zj, mean = nparms$mean, sd = nparms$sd)
    } else if (type == 'w') {
        dfk$cdf <- stats::pweibull(dfk$zj, shape=wparms$shape, scale=wparms$scale)
    } else if (type == 'j') {
        dfk$cdf <- SuppDists::pJohnson(dfk$z, jparms)
    }
    
    ## QQ Plot comparing ECDF to CDF
    with(dfk, plot(Fhat, cdf, xlab='Estimated CDF from Observations', ylab='Theoretical CDF from Fit', main='QQ Plot'))
    abline(0,1, col='red')
  }
  
  if (grepl('d', plotit)) {
    ## dlh view of how this would be done
    ## create QQ Plot
    df$xtheoretical <- with(df, stats::qweibull(ppoints(length(i)), shape=wparms$shape, scale=wparms$scale))
    with(df, plot(zj, xtheoretical, xlab='Observed value, x', ylab='Expected Value', main='DLH QQ Plot'))
    abline(0,1, col='red')
    
  }
}




qqplot_censored_test <- function() {

  ## setup
  os <- .Platform$OS.type
  if (os == 'windows') {
    ## load generic modules
    source("F:\\Documents\\01_Dave\\Programs\\GitHub_home\\R-setup\\setup.r")
    ## identify working folder
    path <- c("f:/Documents/01_Dave/Programs/GitHub_home/Finance/")
  } else {
    ## os == unix
    source('~/GitHub_repos/R-setup/setup.r')
    path <- c('~/GitHub_repos/Finance/')
  }

  
  ## method used in qqplot_censored() is based on 
  ## (1) Estimation of CDF for left censored data as described on p. 5 of:
  ##     https://pdixon.stat.iastate.edu/stat505/Chapter%2011.pdf
  ## (2) Plotting that against the theoretical CDF
  dfin <- '
x 	type
12	Known
6	Known
10	Censored
6	Known
5	Censored
5	Known
5	Censored
3	Known
'
  df <- readall(dfin)
  x      <- df[df$type == 'Known', 'x']
  x.low  <- 0
  x.high <- df[df$type == 'Censored', 'x']
  xcen <- data.frame(x.low=x.low, x.high=x.high)
  qqplot_censored(x,xcen,plotit='xheq')
  
  
  ## install.packages('EnvStats')
  df <- EnvStats::Olympic.NH4.df
  df$x <- log(df$NH4.mg.per.L)
  df$type <- ifelse(df$Censored == 'TRUE', 'Censored', 'Known')
  x      <- df[df$type == 'Known', 'x']
  x.low  <- log(1E-6)
  x.high <- df[df$type == 'Censored', 'x']
  xcen <- data.frame(x.low=x.low, x.high=x.high)
  qqplot_censored(x, xcen, type='n', plotit='xheq')    # good match to ECDF in reference

  ## unfortunately, above qq plot does not look like the one in the text
  ## but maybe that is OK because the above is CDF vs. ECDF while text has quantiles
  plotspace(1,2)
  qqplot_censored(x, xcen, type='n')
  with(df, EnvStats::qqPlotCensored(NH4.mg.per.L, Censored,
                                    distribution = "lnorm", add.line = TRUE, main = ""))
  
  
  
  ## install.packages('EnvStats')
  df <- EnvStats::Olympic.NH4.df
  df$x <- df$NH4.mg.per.L
  df$type <- ifelse(df$Censored == 'TRUE', 'Censored', 'Known')
  x      <- df[df$type == 'Known', 'x']
  x.low  <- 0
  x.high <- df[df$type == 'Censored', 'x']
  xcen <- data.frame(x.low=x.low, x.high=x.high)
  plotspace(2,2)
  qqplot_censored(x, xcen, type='n', plotit='q')  
  qqplot_censored(x, xcen, type='w', plotit='q')   
  qqplot_censored(x, xcen, type='j', plotit='q') 
  qqplot_censored(x, xcen, type='j', plotit='r') 
  
  
  ## R Package containing over 120 data sets used in the text Statistical Methods 
  ## for Reliability Data by Meeker and Escobar (1998) - Auburngrads/SMRD.data
  ## install.packages("remotes")
  ## remotes::install_github("Auburngrads/SMRD.data")
  
  ## shock absorber data w/ right censored failures
  df <- SMRD.data::shockabsorber   # p.630 of Meeker
  df$x <- df$miles
  df$type <- ifelse(df$event == 'Failure', 'Known', 'Censored')
  x      <- df[df$type == 'Known', 'x']
  x.low  <- 0
  x.high <- df[df$type == 'Censored', 'x']
  xcen <- data.frame(x.low=x.low, x.high=x.high)
  qqplot_censored(x, xcen, type='j', plotit='xheq')
}



