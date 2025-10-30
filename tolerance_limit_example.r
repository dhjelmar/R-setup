source('~/Programs/GitHub_home/R-setup/setup.r')
source('~/Documents/GitHub/R-setup/setup.r')

# create dataframe of weibull data
num   <- 1000
shape <- 2
scale <- 1
df <- rweibull(num, shape, scale)
df

# create 2x2 space for 4 plots
plotspace(2,2)

# create histogram with normal, Weibull, and Johnson SU fits and tolerance limits
out_hist <- hist_nwj(df, type='nwj', side='upper', sided=1, P=0.99, conf=0.99)

# create qq plots
out_n <- qqplot_nwj(df, type='n')
out_w <- qqplot_nwj(df, type='w')
out_j <- qqplot_nwj(df, type='j')

# distribution free statistics
out_non <- nonparametric.tol(num, conf=0.99, P=0.99)
out_non$message
out_non$df

highpoints <- 2
out_non <- nonparametric.tol(num, conf=0.99, tol.index=num-highpoints)
out_non$message
out_non$df

################################################################################
# 1-sided vs. 2-sided tolerance bounds

set.seed(100)
x <- rnorm(1000,mean=0,sd=1)

###########################

example_2sided <- function(seed=1, n=1000) {
    'requires n=1000 for tolerance::normtol.int() to match k-factor table example'
    set.seed(seed)
    x <- rnorm(n, mean=0, sd=1)
    
    # using tolerance factors (k)
    # tolerance factor table
    # https://ntrs.nasa.gov/api/citations/19670023646/downloads/19670023646.pdf
    # Notes: Uses sample mean and standard deviation rather than MLE fit to those parameters
    #        Considered less accurate then MLE especially for small data sets because more reliant on normality.
    k <- 4.383    # for n=25, conf=0.95, and P=0.99
    k <- 2.719    # n=1000, conf=0.99, P (reliabiilty) = 0.99
    k <- 2.036    # n=1000, conf=0.95, P (reliabiilty) = 0.95
    lower_tol_limit = mean(x) - k * sd(x) 
    upper_tol_limit = mean(x) + k * sd(x)
    lower_tol_limit                                              
    upper_tol_limit                                     
    df <- data.frame(sided=2, alpha=NA, conf=0.95, P=0.95, 
                     tol.lower=lower_tol_limit, tol.upper=upper_tol_limit, 
                     mean=mean(x), k=k, sd=sd(x),
                     description='***** 2-sided tolerance factor table   *****')

    # using tolerance factors (k) with MLE mean and sd
    # tolerance factor table
    # https://ntrs.nasa.gov/api/citations/19670023646/downloads/19670023646.pdf
    # Notes: Uses sample mean and standard deviation rather than MLE fit to those parameters
    #        Considered less accurate then MLE especially for small data sets because more reliant on normality.
    k <- 4.383    # for n=25, conf=0.95, and P=0.99
    k <- 2.719    # n=1000, conf=0.99, P (reliabiilty) = 0.99
    k <- 2.036    # n=1000, conf=0.95, P (reliabiilty) = 0.95
    out <- mle.normal(x)
    lower_tol_limit = out$parms$xbar - k * out$parms$sdev
    upper_tol_limit = out$parms$xbar + k * out$parms$sdev
    lower_tol_limit                                              
    upper_tol_limit                                     
    dfout <- data.frame(sided=2, alpha=NA, conf=0.95, P=0.95, 
                        tol.lower=lower_tol_limit, tol.upper=upper_tol_limit, 
                        mean=out$parms$xbar, k=k, sd=out$parms$sdev,
                        description='***** 2-sided tolerance factor table w/ MLE mean + sdev  *****')
    df <- rbind(df, dfout)
    
    # 2-sided using tolerance package
    P2 <- 0.95
    out <- tolerance::normtol.int(x, side=2, alpha=  1 - 0.95        , P=P2)
    out_normtol_2high <- out$`2-sided.upper`
    dfout <- data.frame(sided=2, alpha=out$alpha, conf=NA, P=out$P, 
                        tol.lower=out$`2-sided.lower`, tol.upper=out$`2-sided.upper`, 
                        mean=out$x.bar, k=NA, sd=NA,
                        description='***** 2-sided tolerance::normtol.int() *****')
    df <- rbind(df, dfout)
    
    # 2-sided using MLE
    sided <- 2
    conf_orig  <- 0.95
    conf <- conf_orig
    alpha <- (1-conf)/(2/sided) 
    out <- mle.normal.tol(x, sided=2, conf=conf, P=P2, plots=FALSE)    
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out_mle_2high <- out$tolerance$tol.upper
    out$tolerance$description <- '***** 2-sided MLE using conf           *****'
    df <- rbind(df, out$tolerance)
    #
    out <- mle.normal.tol(x, sided=2, alpha=alpha, P=P2)  
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- '      2-sided MLE using alpha               '
    df <- rbind(df, out$tolerance)
    #
    out <- mle.normal.tol(x, sided=2, alpha=alpha, P=P2, old=FALSE)  
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- '      2-sided MLE using alpha NEW              '
    df <- rbind(df, out$tolerance)
    
    # df <- df.init(rows=0, columns=c('mle.tol','P','normal.tol','diff','diff_percent'))
    # for (alpha in seq(0.05,0.2,0.01)) {
    # #alpha <- 0.05
    # #for (P2 in seq(0.94,0.96,0.001)) {
    #     tol <- mle.normal.tol(x, sided=2, alpha=alpha, P=P2)$tolerance$tol.upper
    #     df <- rbind(df, data.frame(alpha = alpha,
    #                                P = P2,
    #                                mel.tol = tol,
    #                                norm.tol = out_normtol_2high, 
    #                                diff = tol - out_normtol_2high, 
    #                                diff_percent = (tol / out_normtol_2high -1) * 100))
    # }
    # df
    
    # 1-sided using tolerance package
    out <- tolerance::normtol.int(x, side=1, alpha=  1 - 0.95        , P=0.95) 
    out_normtol_1high <- out$`1-sided.upper`
    dfout <- data.frame(sided=1, alpha=out$alpha, conf=NA, P=out$P, 
                        tol.lower=out$`1-sided.lower`, tol.upper=out$`1-sided.upper`, 
                        mean=out$x.bar, k=NA, sd=NA,
                        description='***** 1-sided tolerance::normtol.int()               *****')
    df <- rbind(df, dfout)
    
    # 1-sided at same, confidence and coverage
    conf  <- conf_orig
    out <- mle.normal.tol(x, side='lower',  sided=1, conf=conf, P=1-P2) 
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- '***** 1-sided MLE using conf (not as low as 2-sided) *****'
    df <- rbind(df, out$tolerance)
    out <- mle.normal.tol(x, side='upper',  sided=1, conf=conf, P=P2)
    out_mle_1high <- out$tolerance$tol.upper
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- '***** 1-sided MLE using conf (not as low as 2-sided) *****'
    df <- rbind(df, out$tolerance)
    
    
    # 1-sided at same, confidence and coverage
    conf  <- conf_orig
    out <- mle.normal.tol(x, side='lower',  sided=1, conf=conf, P=1-P2, old=FALSE) 
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- '***** 1-sided MLE_NEW using conf (not as low as 2-sided) *****'
    df <- rbind(df, out$tolerance)
    out <- mle.normal.tol(x, side='upper',  sided=1, conf=conf, P=P2, old=FALSE)
    out_mle_1high <- out$tolerance$tol.upper
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- '***** 1-sided MLE_NEW using conf (not as low as 2-sided) *****'
    df <- rbind(df, out$tolerance)
    
    # 1-sided at same, confidence and wider coverage since 100% on one side 
    conf  <- conf_orig
    P1L <-     (1 - P2)/2     # 0.025
    P1U <- 1 - (1 - P2)/2     # 0.975
    out <- mle.normal.tol(x, side='lower',  sided=1, conf=conf, P=P1L) 
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- 'MLE: 1-sided still not as low for 95% lower / 100% upper (P=0.025)'
    df <- rbind(df, out$tolerance)
    out <- mle.normal.tol(x, side='upper',  sided=1, conf=conf, P=P1U)
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- 'MLE: 1-sided still not as high for 100% lower / 95% upper (P=0.975'
    df <- rbind(df, out$tolerance)
    
    # to get identical results with 1-sided, use same alpha (which leads to different confidence for 1-sided)
    # also need higher coverage since 100% is covered on one side
    sided <- 1
    P1L <-     (1 - P2)/2 
    P1U <- 1 - (1 - P2)/2 
    out <- mle.normal.tol(x, side='lower',  sided=1, alpha=alpha, P=P1L) 
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- 'MLE: match 2-sided with 1-sided using alpha'
    df <- rbind(df, out$tolerance)
    out <- mle.normal.tol(x, side='upper',  sided=1, alpha=alpha, P=P1U)
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- 'MLE: match 2-sided with 1-sided using alpha'
    df <- rbind(df, out$tolerance)

    # if do with confidence, then:
    conf <- 1 - alpha*(sided/2)
    out <- mle.normal.tol(x, side='lower',  sided=1, conf=conf, P=P1L)
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- 'MLE: match 2-sided with 1-sided using conf'
    df <- rbind(df, out$tolerance)
    out <- mle.normal.tol(x, side='upper',  sided=1, conf=conf, P=P1U) 
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- 'MLE: match 2-sided with 1-sided using conf'
    df <- rbind(df, out$tolerance)

    return(list(df=df, mle_vs_normtol_2=out_mle_2high/out_normtol_2high-1, mle_vs_normtol_1=out_mle_1high/out_normtol_1high-1))
}
example_2sided(seed=1,n=1000)
