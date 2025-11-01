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
    df <- data.frame(sided=2, alpha.eff=NA, conf=0.95, P=0.95, 
                     tol.lower=lower_tol_limit, tol.upper=upper_tol_limit, 
                     mean=mean(x), k=k, sd=sd(x),
                     description='***** 2-sided tolerance factor table   *****')
    
    # 2-sided using tolerance package
    P     <- 0.95
    conf  <- 0.95
    alpha <- 1 - conf
    out <- tolerance::normtol.int(x, side=2, alpha=  alpha       , P=P)
    out_normtol_2high <- out$`2-sided.upper`
    dfout <- data.frame(sided=2, alpha.eff=out$alpha, conf=NA, P=out$P, 
                        tol.lower=out$`2-sided.lower`, tol.upper=out$`2-sided.upper`, 
                        mean=out$x.bar, k=NA, sd=NA,
                        description='***** 2-sided tolerance::normtol.int() *****')
    df <- rbind(df, dfout)
    
    # 2-sided using tolerance package but estimated with 1-sided
    out <- tolerance::normtol.int(x, side=1, alpha=  alpha/2        , P=P+(1-P)/2)
    out_normtol_2high <- out$`1-sided.upper`
    dfout <- data.frame(sided=1, alpha.eff=out$alpha, conf=NA, P=out$P, 
                        tol.lower=out$`1-sided.lower`, tol.upper=out$`1-sided.upper`, 
                        mean=out$x.bar, k=NA, sd=NA,
                        description='***** 2-sided conservative estimate with 1-sided tolerance::normtol.int() *****')
    df <- rbind(df, dfout)
    
    # 2-sided MLE
    out <- mle.normal.tol(x, sided=2, conf=conf, P=P, plots=FALSE)    
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out_mle_2high <- out$tolerance$tol.upper
    out$tolerance$description <- '***** 2-sided conservative estimate with MLE using conf           *****'
    df <- rbind(df, out$tolerance)
    
    # 1-sided using tolerance package
    out <- tolerance::normtol.int(x, side=1, alpha=  alpha        , P=P)
    dfout <- data.frame(sided=1, alpha.eff=out$alpha, conf=NA, P=out$P, 
                        tol.lower=out$`1-sided.lower`, tol.upper=out$`1-sided.upper`, 
                        mean=out$x.bar, k=NA, sd=NA,
                        description='***** 1-sided tolerance::normtol.int() *****')
    df <- rbind(df, dfout)
    
    # 1-sided MLE
    # upper
    out <- mle.normal.tol(x, sided=1, conf=conf, P=P, plots=FALSE)    
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- '***** 1-sided MLE using conf           *****'
    df <- rbind(df, out$tolerance)
    # lower
    out <- mle.normal.tol(x, side='lower', sided=1, conf=conf, P=1-P, plots=FALSE)    
    out$tolerance$mean=out$params$xbar
    out$tolerance$k=NA
    out$tolerance$sd=out$params$sdev
    out$tolerance$description <- '***** 1-sided MLE using conf           *****'
    df <- rbind(df, out$tolerance)
    
    return(df)
}
example_2sided(seed=1,n=1000)
