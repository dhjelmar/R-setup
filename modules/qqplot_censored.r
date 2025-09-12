# something else to try: EnvStats::qqPlotCensored()
# https://www.rdocumentation.org/packages/EnvStats/versions/2.3.1/topics/qqPlotCensored


qqplot_censored <- function(x=NA, xcen=NA, xcen0=NA, type='n', nfit='mle', wfit='mle', jfit='mle', 
                            algorithm='NLOPT_LN_COBYLA', onebyone=FALSE,
                            xtol_rel =1.0e-4, stopval=0.005, maxeval=500000,
                            mainadder=NULL, plotit='xqh', plotcen=FALSE) {
    ## x    = 'Known'
    ## xcen = dataframe of x.low, x.high pairs
    ## xcen0 = (xcen$x.low + xcen$x.high)/2 (default)
    ##       = vector of initial guess for plotting locations for censored points
    ## type = 'n' to fit normal distribution
    ##        'w' to fit Weibull distribution
    ##        'j' to fit Johnson SU distribution
    ##        'nwj' to fit all 3
    ## nfit = 'mle', 'standard', or vector
    ## wfit = 'mle', 'standard', or vector
    ## jfit = 'mle', 'standard', or vector
    ## mainadder = NULL (default)
    ##           = text to add to header
    ## plotit = 'x' plots x vs. index
    ##          'h' plots histogram of known x values and pdf
    ##          'e' plots ECDF
    ##          'c' plots only the ECDF vs CDF QQ plot
    ##          'q' plots only the standard qq plot
    ##          can contain any combination of x, h, e, c, and q (e.g., plotit='cq')
    ## plotcen = FALSE (default) does not include censored points on the qq plot
    
    
    ## ##########################################################
    ## FIT THE DATA TO DETERMINE THEORETICAL PDF PARAMETERS
    ## AND REMOVE KNOWN DATA IF ANY FROM XCEN
    ## ##########################################################
    out <- hist_nwj(x, xcen, type=type, nfit=nfit, wfit=wfit, jfit=jfit, tolerance=FALSE, plot=FALSE)
    if (type == 'n') {
        parms <- out$nparms
    } else if (type == 'w') {
        parms <- out$wparms
    } else {
        parms <- out$jparms
    }
    x <- out$known
    xcen <- out$censored
    
    
    ## Order xcen from smallest to largest interval with -inf and +inf
    ## set to user specified values (default min and max data values).
    ## sort from lowest to highest with censored values before known
    ## put x and xcen into single dataframe
    df1 <- data.frame(x.low =      x    , x.high =      x     , type='Known')
    df2 <- data.frame(x.low = xcen$x.low, x.high = xcen$x.high, type='Censored')
    df <- rbind(df1,df2)
    df <- df[order(df$type  , decreasing = FALSE),]
    df <- df[order(df$x.high, decreasing = FALSE),]
    df$index <- c(1:nrow(df))
    
    ## For each censored point in order from smallest to largest interval, 
    ## find x' in interval for which a new fit makes no change to fit parameters
    ## (optimization routine for 1 objective function with 4 constraints so just new - old for parameter = 0;
    ## minimize sum( (pnew-pold)^2 ) - error)
    
    ## create variables for range of censored values and xprime
    xcen$range <- xcen$x.high - xcen$x.low
    df$range <- df$x.high - df$x.low
    df$xprime <- (df$x.high + df$x.low)/2            # this is just a starting position

    ## sort by range (smallest to largest) and identify which are censored
    ## sorting xcen by smallest to largest range only important if setting values 1 at a time
    xcen <- xcen[order(xcen$range, decreasing = FALSE),]
    df <- df[order(df$range, decreasing = FALSE),]
    cenum <- which(df$type == 'Censored')
    

    
    ## #################################################
    ## #################################################
    
    ## find x' in interval for which a new fit makes no change to fit parameters
    ## install.packages('nloptr')

    ## options
    opts <- list('algorithm'   = algorithm,
                 'xtol_rel'    = xtol_rel,
                 'maxeval'     = maxeval,
                 'print_level' = 1 ,
                 'stopval'     = stopval)    # stop when objective function hits this

    
    if (isTRUE(onebyone)) {
        ## #################################################
        ## attempt to set xprime one xcen row at a time
        ## #################################################
        
        ## define function to minimize
        objfunc <- function(xtarget, x, xcen, type, parmsin) {
            x <- c(x, xtarget)
            if (type == 'n') {
                parmsin <- unlist(parmsin)
                out <- mle.normal(x, xcen)
                parms <- unlist(out$parms)
                ss <- sum((parms - parmsin)^2)
                debug <- TRUE
                if (isTRUE(debug)) cat('xtarget=', signif(xtarget,11),
                                       'mean=', signif(parms[['xbar']],11),
                                       'sd=', signif(parms[['sdev']],11),
                                       'ss=', signif(ss,11), "\n")
                
            } else if (type == 'j') {
                parmsin <- unlist(parmsin[1:4])
                out <- mle.johnsonsu(x, xcen)
                parms <- unlist(out$parms[1:4])
                ss <- sum((parms - parmsin)^2)
                debug <- TRUE
                if (isTRUE(debug)) cat('xtarget=', signif(xtarget,11),
                                       'gamma=', signif(parms[['gamma']],11),
                                       'delta=', signif(parms[['delta']],11),
                                       'xi   =', signif(parms[['xi']], 11),
                                       'lambda=', signif(parms[['lambda']],11),
                                       'ss=', signif(ss,11), "\n")
            }
        }
        
        x_plus <- df[-cenum]
        xcen_minus <- df[cnum]
        for (i in cennum) {
            ## move 1st row of dataframe to x.low, x.high, and xprime
            x.low  <- xcen_minus[1,'x.low']
            x.high <- xcen_minus[1,'x.high']
            xprime <- xcen_minus[1,'x.high']
            xcen_minus <- xcen_minus[-1,]
            ## objective function
            cat('x.low=', signif(x.low,11),
                'x.high=', signif(x.high,11),
                'parms=', signif(parms,11), "\n")
            out <- nloptr::nloptr(x0     = xprime,   # initial guess (mist be 1st argument here and of function)
                                  x      = x,       # data for function
                                  xcen   = xcen,    # data for function
                                  type   = type,     # either n, w, or j for normal, Weibull, Johnson SU
                                  parmsin = parms,
                                  eval_f = objfunc,     # objective function
                                  lb     = x.low,   # lower bound for x0
                                  ub     = x.high,      # upper bounds
                                  ## eval_grad_f =  # gradient function
                                  ## eval_g_ineq # inequality constraints
                                  ## eval_g_eq   # equality constraints
                                  opts   = opts)      # optimization options
        }

    } else {    
        
        
        ## #################################################
        ## attempt to set xprime for all xcen rows at once
        ## #################################################
        
        
        ## https://cran.r-project.org/web/packages/nloptr/nloptr.pdf
        ## 
        
        objfunc <- function(xtarget, x, type, parmsin) {
            if (type == 'n') {
                parmsin <- unlist(parmsin)
                out <- mle.normal(c(x, xtarget))
                parms <- unlist(out$parms)
                ss <- sum((parms - parmsin)^2)
            } else if (type == 'w') {
                parmsin <- unlist(parmsin)
                out <- mle.weibull(c(x, xtarget))
                parms <- unlist(out$parms)
                ss <- sum((parms - parmsin)^2)
                ## if (isTRUE(debug)) cat('xtarget=', signif(xtarget,11),
                ##                        #'shape=', signif(parms[['shape']],11),
                ##                        #'scale=', signif(parms[['scale']],11),
                ##                        'parms=', signif(parms,11),
                ##                        'ss=', signif(ss,11), "\n")

            } else if (type == 'j') {
                parmsin <- unlist(parmsin[1:4])
                out <- mle.johnsonsu(c(x, xtarget))
                parms <- unlist(out$parms[1:4])
                ss <- sum((parms - parmsin)^2)
                ## if (isTRUE(debug)) cat('xtarget=', signif(xtarget,11),
                ##                        #'gamma=', signif(parms[['gamma']],11),
                ##                        #'delta=', signif(parms[['delta']],11),
                ##                        #'xi   =', signif(parms[['xi']], 11),
                ##                        #'lambda=', signif(parms[['lambda']],11),
                ##                        'parms=', signif(parms[1:4],11),
                ##                        'ss=', signif(ss,11), "\n")
            }
            debug <- TRUE
            if (isTRUE(debug)) cat('parmsin=', signif(parmsin,11),
                                   'parms=', signif(parms,11),
                                   'ss=', signif(ss,11), "\n")
            return(ss)
        }

        cat('x.low=', signif(xcen$x.low,11),
            'x.high=', signif(xcen$x.high,11),
            'parms=', signif(unlist(parms[1:4]),11), "\n")
        if (is.na(xcen0)) {
            ## no initial guess for xcen locations
            xcen0 <- (xcen$x.low + xcen$x.high)/2
        }
        out <- nloptr::nloptr(x0     = xcen0,    # initial guess
                              x      = x,        # data for function
                              type   = type,     # either n, w, or j for normal, Weibull, Johnson SU
                              parmsin = parms,
                              eval_f = objfunc,     # objective function
                              lb     = xcen$x.low,   # lower bound for x0
                              ub     = xcen$x.high,      # upper bounds
                              ## eval_grad_f =  # gradient function
                              ## eval_g_ineq # inequality constraints
                              ## eval_g_eq   # equality constraints
                              opts   = opts)      # optimization options
        xcen$xprime <- out$solution
        xprime <- c(x, out$solution)
        df1 <- data.frame(x.low = x, x.high=x, type='Known', range=0, xprime=x)
        df2 <- xcen
        df2$type <- 'Censored'
        df2$range <- df2$x.high - df2$x.low
        df <- fastmerge(df1, df2)
        df$range <- NULL

        ## sort points in increasing order of xprime
        df <- df[order(df$xprime, decreasing = FALSE),]
        row.names(df) <- 1:nrow(df)
        df$index <- 1:nrow(df)

        if (grepl('x', plotit)) {
            ## PLOT THE DATA TO SEE CENSORED AND KNOWN
            ## plotfit(df$index, df$xplot, interval='noline')
            ## addfit(df[df$type == 'Known',]$index, df[df$type == 'Known',]$xplot, pch=16, interval='noline')
            ## legend('topleft', legend=c('Censored', 'Known'), pch=c(1,16))

            ## plot(df$xprime, ylim=range(df$x.low, df$x.high))
            ## points(df[df$type == 'Known',]$index, df[df$type == 'Known',]$xprime, pch=16)
            ## lines(df$x.low)
            ## lines(df$x.high)
            ## legend('topleft', legend=c('Censored (lines for range)', 'Known'), pch=c(1,16))
            
            plotbar(y=df$xprime, ylow=df$x.low, yhigh=df$x.high, ylab='x')
            
        }

        ## ##########################################################
        ## create qqplot
        ## option to plot censored points as open circles
        ## ##########################################################

        if (type == 'n') {
            ## make normal QQ plot
            main <- paste('Normal QQ Plot', mainadder, sep=" ")
            ##df$xtheoretical <- qnorm(ppoints(nrow(df)), mean = mean(x), sd = sd(x))
            df$xtheoretical <- qnorm(ppoints(nrow(df)), mean = parms$mean, sd = parms$sd)

        } else if (type == 'w') {
            ## Weibull
            main <- paste('Weibull QQ Plot', mainadder, sep=" ")
            df$xtheoretical <- stats::qweibull(ppoints(nrow(df)), shape=parms$shape, scale=parms$scale)
            
        } else {
            ## Johnson SU
            main <- paste('Johnson QQ Plot', mainadder, sep=" ")
            df$xtheoretical <- SuppDists::qJohnson(ppoints(nrow(df)), parms)
            
        }

        ## make the qqplot
        cenum <- which(df$type == 'Censored')
        plot(df[-cenum, 'xprime'], df[-cenum, 'xtheoretical'],
             xlim = range(df$xprime), ylim=range(df$xtheoretical),
             pch=16, xlab='Observed value, x', ylab='Expected Value', main=main)
        if (isTRUE(plotcen)) points(df[cenum, 'xprime'], df[cenum, 'xtheoretical'], pch=0)
        abline(0,1, col='red')
        
    }
    

    ## create hist for x and x' combined but plot both
    if (grepl('h', plotit)) {
        
        ## histogram
        col1 <- rgb(0,0,1,1/4)
        col2 <- rgb(1,0,0,1/4)
        hist(df$xprime           , col=col1, freq=FALSE)
        hist(df[-cenum, 'xprime'], col=col2, freq=FALSE, add = TRUE)
        legend('topright',
               legend=c('all data', 'known only', 'censored fit', 'positioned fit'),
               col=c(col1,col2, 'black', 'red'),
               lty=c( 1,  1, 1, 2),
               lwd=c(10, 10, 1, 1))
        
        ## add fit for known + censored
        ## add distributions
        xmax <- max(df$x.high)
        xmin <- min(df$x.low)
        chuncks <- (xmax-xmin)/1000
        xrange  <- seq(xmin,xmax,by=chuncks)
        if (grepl('j', type)) {
            xdensity <- SuppDists::dJohnson(xrange, parms)
        } else if (grepl('w', type)) {
            xdensity <- stats::dweibull(xrange, shape=parms$shape, scale=parms$scale)
        } else {
            xdensity <- stats::dnorm(xrange, parms$mean, parms$sd)
        }
        lines(x=xrange, y=xdensity, col='black', lty=1)
        
        
        ## add fit for xprime
        ## add distributions
        outprime <- hist_nwj(df$xprime, type=type, nfit=nfit, wfit=wfit, jfit=jfit, tolerance=FALSE, plot=FALSE)
        if (type == 'n') {
            parmsprime <- outprime$nparms
        } else if (type == 'w') {
            parmsprime <- outprime$wparms
        } else {
            parmsprime <- outprime$jparms
        }
        if (grepl('j', type)) {
            xdensity <- SuppDists::dJohnson(xrange, parmsprime)
        } else if (grepl('w', type)) {
            xdensity <- stats::dweibull(xrange, shape=parmsprime$shape, scale=parmsprime$scale)
        } else {
            xdensity <- stats::dnorm(xrange, parmsprime$mean, parmsprime$sd)
        }
        lines(x=xrange, y=xdensity, col='red', lty=2)
        
        
    }
    
    return(list(data = df, x0=out$solution, cenum=cenum, parms=parms, parmsprime=parmsprime))
}  

qqplot_censored_test <- function() {
    ## create Johnson SU dataset
    set.seed(1)
    jparms <- list(gamma=-3.3, delta=5.3, xi=1.8, lambda=1.9, type='SU')
    x <- ExtDist::rJohnsonSU(100, param=jparms)
    ## add very high known at x = 8
    x <- c(x,6)
    type <- 'n'
    ## add very high known at x = 8 and a censored point with a wide range
    ##xcen=data.frame(x.low=0, x.high=c(3.5, 3.5, 4, 4, min(x)))
    ## xcen=data.frame(x.low=0, x.high=rep(4,10))
    ## xcen=data.frame(x.low=0, x.high=rep(4,10))
    xcen=data.frame(x.low=c(rep(0, 10), 3.8), x.high=c(rep(4,10),8))   # works for type='j' but is slow
    xcen
    ## make plots
    plotspace(2,2)
    out <- qqplot_censored(x, xcen, xcen0=NA, type=type, plotit='xqh', plotcen=TRUE, 
                           algorithm='NLOPT_LN_COBYLA',  onebyone=FALSE,
                           xtol_rel=1.0e-4, stopval=1E-4, maxeval=1000)
    
    ## algorithm method choices: (G=Global; L=Local; D=Derivative Required)
    ## -------------------------
    ## NLOPT_GN_DIRECT
    ## NLOPT_GN_DIRECT_L
    ## NLOPT_GN_DIRECT_L_RAND
    ## NLOPT_GN_DIRECT_NOSCAL
    ## NLOPT_GN_DIRECT_L_NOSCAL
    ## NLOPT_GN_DIRECT_L_RAND_NOSCAL
    ## NLOPT_GN_ORIG_DIRECT
    ## NLOPT_GN_ORIG_DIRECT_L
    ## NLOPT_GN_CRS2_LM
    ## NLOPT_GN_MLSL
    ## NLOPT_GN_MLSL_LDS
    ## NLOPT_GN_ISRES
    ## NLOPT_GN_ESCH

    ## NLOPT_LN_PRAXIS
    ## NLOPT_LN_COBYLA
    ## NLOPT_LN_NEWUOA
    ## NLOPT_LN_NEWUOA_BOUND
    ## NLOPT_LN_NELDERMEAD
    ## NLOPT_LN_SBPLX
    ## NLOPT_LN_AUGLAG
    ## NLOPT_LN_AUGLAG_EQ
    ## NLOPT_LN_BOBYQA

    ## NLOPT_GD_STOGO
    ## NLOPT_GD_STOGO_RAND
    ## NLOPT_GD_MLSL_LDS
    ## NLOPT_GD_MLSL
    
    ## NLOPT_LD_MMA
    ## NLOPT_LD_CCSAQ
    ## NLOPT_LD_SLSQP
    ## NLOPT_LD_LBFGS_NOCEDAL
    ## NLOPT_LD_LBFGS
    ## NLOPT_LD_VAR1
    ## NLOPT_LD_VAR2
    ## NLOPT_LD_TNEWTON
    ## NLOPT_LD_TNEWTON_RESTART
    ## NLOPT_LD_TNEWTON_PRECOND
    ## NLOPT_LD_TNEWTON_PRECOND_RESTART
    ## NLOPT_LD_AUGLAG
    ## NLOPT_LD_AUGLAG_EQ
    
    
}













#####################################################################
#####################################################################
#####################################################################
#####################################################################


qqplot_censored_old <- function(x=NA, xcen=NA, type='nwj', nfit='mle', wfit='mle', jfit='mle', mainadder=NULL, 
                                plotit='xhecq', plotcen=FALSE) {
    
    
    ## put x and xcen into single dataframe
    df1 <- data.frame(x.low =      x    , x.high =      x     , type='Known')
    df2 <- data.frame(x.low = xcen$x.low, x.high = xcen$x.high, type='Censored')
    df <- rbind(df1,df2)
    df$x <- df$x.high
    
    
    if (plotit == 'xhecq') plotspace(2,3)
    
    if (grepl('x', plotit)) {
        ## ##########################################################
        ## PLOT THE DATA TO SEE CENSORED AND KNOWN
        ## ##########################################################  plotspace(3,1)
        ## sort from lowest to highest with censored values before known
        df <- df[order(df$type, decreasing = FALSE),]
        df <- df[order(df$x, decreasing = FALSE),]
        df$index <- c(1:nrow(df))
        plotfit(df$index, df$x, interval='noline')
        addfit(df[df$type == 'Known',]$index, df[df$type == 'Known',]$x, pch=16, interval='noline')
        legend('topleft', legend=c('Censored', 'Known'), pch=c(1,16))
    }
    
    
    ## ##########################################################
    ## FIT THE DATA TO DETERMINE THEORETICAL PDF PARAMETERS
    ## ##########################################################
    if (grepl('h', plotit)) {
        out <- hist_nwj(x, xcen, type=type, nfit=nfit, wfit=wfit, jfit=jfit, tolerance=FALSE)
    } else {
        out <- hist_nwj(x, xcen, type=type, nfit=nfit, wfit=wfit, jfit=jfit, tolerance=FALSE, plot=FALSE)
    }
    nparms <- out$nparms
    wparms <- out$wparms
    jparms <- out$jparms
    
    
    ## ##########################################################
    ## ESTMATE CDF (ECDF)
    ## ##########################################################
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
    if (grepl('c', plotit)) {
        plotfit(dfk$zj, dfk$Fhat, main='Empirical CDF', pch=16, interval='line')
    }
    
    
    if (grepl('e', plotit)) {
        ## ##########################################################
        ## THEORETICAL CDF
        ## ##########################################################
        
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
        main <- paste('CDF vs ECDF', mainadder, sep=" ")
        with(dfk, plot(Fhat, cdf, xlab='Estimated CDF from Observations', ylab='Theoretical CDF from Fit', pch=16, main=main))
        abline(0,1, col='red')
    }
    
    if (grepl('q', plotit)) {
        ## dlh view of how this would be done
        
        ## sort data from lowest to highest
        ## put left censored data before known of same value
        ## where left censored have x.low = 0
        df <- df[order(df$type, decreasing = FALSE),]
        df <- df[order(df$x, decreasing = FALSE),]
        df$i <- c(1:nrow(df))
        
        ## theoretical CDF
        if (type == 'n') {
            df$expected <- qnorm(ppoints(length(df$i)), mean = nparms$mean, sd = nparms$sd)
        } else if (type == 'w') {
            df$expected <- stats::qweibull(ppoints(length(df$i)), shape=wparms$shape, scale=wparms$scale)
        } else if (type == 'j') {
            df$expected <- SuppDists::qJohnson(ppoints(length(df$i)), jparms)
        }
        
        ## create QQ Plot but only plot the known values
        ## 'known' is NA for censored so that works easily
        main <- paste('QQ Plot', mainadder, sep=" ")

        if (isFALSE(plotcen)) {
            with(df, plot(known, expected, xlab='Observed value, x', ylab='Expected Value', pch=16, main=main))
            ## mtext('Censored data used to determine plotting positions but not plotted.', side=3, line=0.5, cex=.75, col='black')
            mtext('Censored data used to determine', side=3, line=0.75, cex=.75, col='black')
            mtext('plotting positions but not plotted.', side=3, line=0, cex=.75, col='black')
            abline(0,1, col='red')
            
        } else {
            ## plot all known and unknown points as open circles
            with(df, plot(x.high, expected, xlab='Observed value, x', ylab='Expected Value', pch=1, main=main))
            ## replot known points as closed circles
            with(df[df$type == 'Known',], addfit(x.high, expected, pch=16, interval='noline'))
            ## add legend
            legend('topleft', legend=c('Censored', 'Known'), pch=c(1,16))
            abline(0,1, col='red')
        }
    }
    
    return(list(df=df, dfk=dfk))
}







qqplot_censored_old_test <- function() {

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


    ## ####################################################################
    ## TEST CASE 1
    ## ####################################################################
    
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
    qqplot_censored(x,xcen,type='n',plotit='xhecq')
    qqplot_censored(x,xcen,type='n',plotit='q',plotcen=TRUE)
    
    
    
    ## ####################################################################
    ## TEST CASE 2 (from p. 11 of same reference chapter as test case 1)
    ## ####################################################################
    
    ## install.packages('EnvStats')
    df <- EnvStats::Olympic.NH4.df
    df$x <- log(df$NH4.mg.per.L)
    df$type <- ifelse(df$Censored == 'TRUE', 'Censored', 'Known')
    x      <- df[df$type == 'Known', 'x']
    x.low  <- log(1E-11)
    x.high <- df[df$type == 'Censored', 'x']
    xcen <- data.frame(x.low=x.low, x.high=x.high)
    qqplot_nwj(x, xcen, type='n', plotit='xhecq')    # good match to ECDF in reference
    qqplot_nwj(x, xcen, type='n', plotit='q', plotcen=TRUE)
    
    ## unfortunately, above qq plot does not look like the one in the text
    ## but maybe that is OK because the above is CDF vs. ECDF while text has quantiles
    with(df, EnvStats::qqPlotCensored(NH4.mg.per.L, Censored,
                                      distribution = "lnorm", add.line = TRUE, 
                                      main = "from text"))
    
    
    ## ####################################################################
    ## TEST CASE 3 - compare normal, weibull, and johnson su fits
    ## ####################################################################
    
    
    ## install.packages('EnvStats')
    df <- EnvStats::Olympic.NH4.df
    df$x <- df$NH4.mg.per.L
    df$type <- ifelse(df$Censored == 'TRUE', 'Censored', 'Known')
    x      <- df[df$type == 'Known', 'x']
    x.low  <- 0
    x.high <- df[df$type == 'Censored', 'x']
    xcen <- data.frame(x.low=x.low, x.high=x.high)
    plotspace(1,2)
    out <- qqplot_censored(x, xcen, type='nwj', plotit='xh')  
    plotspace(3,2)
    out <- qqplot_censored(x, xcen, type='n', plotit='eq', mainadder='normal')  
    out <- qqplot_censored(x, xcen, type='w', plotit='eq', mainadder='Weibull')  
    out <- qqplot_censored(x, xcen, type='j', plotit='eq', mainadder='Johnson SU')  
    

    ## what if I add high outliers as a known points?
    xplus <- c(x,.3,.5)
    out <- qqplot_censored(xplus, xcen, type='j', plotit='xheq') 
    
    
    ## ####################################################################
    ## TEST CASE 4 - Meeker example
    ## ####################################################################
    
    
    ## R Package containing over 120 data sets used in the text 
    ## "Statistical Methods for Reliability Data" by Meeker and Escobar (1998) - Auburngrads/SMRD.data
    
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
    plotspace(2,2)
    qqplot_censored(x, xcen, type='j', plotit='xheq')
    
    
    ## ####################################################################
    ## TEST CASE 5 - Johnson SU dataset
    ## #################################################################### 
    
    ## create Johnson SU dataset
    set.seed(1)
    jparms <- list(gamma=-3.3, delta=5.3, xi=1.8, lambda=1.9, type='SU')
    x <- ExtDist::rJohnsonSU(1000, param=jparms)
    
    
    ## ##############################
    ## TEST 5A - EFFECT OF ONLY 2 CONSISTENT CENSORED POINTS
    ## first plot histogram and qq plot for known data
    plotspace(3,2)
    hist_nwj(x, type='j', main='all known data')
    qqplot_nwj(x, type='j', mainadder = 'all known data')
    ## added two censored points
    qqplot_censored(x, xcen=data.frame(x.low=0, x.high=3.5), type='j', plotit='xeq')
    qqplot_censored(x, xcen=data.frame(x.low=0, x.high=3.5), type='j', plotit='q', plotcen=TRUE)
    ## as expected, no difference between known and censored qq plots
    
    
    ## ##############################
    ## TEST 5B - EFFECT OF MANY HIGH CENSORED POINTS
    plotspace(3,2)
    hist_nwj(x, type='j', main='all known data')
    outk <- qqplot_nwj(x, type='j', mainadder = 'all known data')
    ## added many censored points
    outc <- qqplot_censored(x, xcen=data.frame(x.low=0, x.high=seq(4.5,5,.01)), type='j', plotit='xheq')
    outc <- qqplot_censored(x, xcen=data.frame(x.low=0, x.high=seq(4.5,5,.01)), type='j', plotit='q', plotcen=FALSE)
    outc <- qqplot_censored(x, xcen=data.frame(x.low=0, x.high=seq(4.5,5,.01)), type='j', plotit='q', plotcen=TRUE)
    ## the CDF vs ECDF plot shows that the fit is still good but that may be because it is insensitive to very high observations
    ## as expected, qq plot shows that the fit is no good because there are too many high censored points
    
    
    ## ##############################
    ## TEST 5C - EFFECT OF ONE VERY HIGH CENSORED POINT
    plotspace(3,2)
    hist_nwj(x, type='j', main='all known data')
    outk <- qqplot_nwj(x, type='j', mainadder = 'all known data')
    ## added one high censored point
    outc <- qqplot_censored(x, xcen=data.frame(x.low=7.99999, x.high=8), type='j', plotit='xheq')
    

    ## ##############################
    ## TEST 5D - EFFECT OF ONE VERY HIGH KNOWN POINT
    plotspace(3,2)
    hist_nwj(x, type='j', main='all known data')
    outk <- qqplot_nwj(x, type='j', mainadder = 'all known data')
    ## added many censored points
    outc <- qqplot_censored(c(x,8), xcen=data.frame(x.low=0, x.high=3.5), type='j', plotit='xheq')

    
    ## #############################
    ## Comparison of very high censored to very high known
    plotspace(3,2)
    ## no high
    out <- qqplot_censored(x, xcen=data.frame(x.low=0, x.high=3.5), type='j', plotit='eq', plotcen=TRUE)
    ## very high censored <- no effect because it is treated as 0 to 8 even though it is really 8
    out <- qqplot_censored(x, xcen=data.frame(x.low=7.99999, x.high=8), type='j', plotit='eq', plotcen=TRUE)
    ## very high known
    out <- qqplot_censored(c(x,8), xcen=data.frame(x.low=0, x.high=3.5), type='j', plotit='eq', plotcen=TRUE)
    ## The CDF vs ECDF plot shows that even the fit with the known high point is good.
    ## That may be because it is insensitive to very high observations (i.e., max CDF = 1 regardless of how high the value is).
    ## As expected, 2nd and 3rd qq plot shows that the fit is no good (at least if censored are shown) 
    ## because the high point is too high censored points.
    
    
    ## ############################
    ## Another data set
    ## ############################
    df <- readall('F:\\Documents\\01_Dave\\Engineering\\STATISTICS\\data_censored.xlsx')
    x      <- as.numeric(df[df$Type == 'known'   ,]$Known)
    x.low  <- as.numeric(df[df$Type == 'censored',]$Low)
    x.high <- as.numeric(df[df$Type == 'censored',]$High)
    xcen <- data.frame(x.low=x.low, x.high=x.high)
    tol <- mle.johnsonsu.tol(x)
    tol <- mle.johnsonsu.tol(x,xcen)
    histout <- hist_nwj(x, xcen, type='j')
    plotspace(1,2)
    out <- qqplot_censored(x, xcen, type='j', plotit='eq', plotcen=FALSE)
    
    
    ## https://stackoverflow.com/questions/41968606/left-censoring-for-survival-data-in-r
    ## install.packages('survival')
    AF_at_baseline<-c(1,0,1,0,0) #where 1 indicates left censoring
    Followup_time <- c(NA, 3, NA, 15, 7)
    Followup_time2 <- c(11, NA, 8 ,15, 7)
    xcen <- data.frame(x.low = Followup_time, x.high = Followup_time2)
    Surv.Obj <- survival::Surv(xcen$x.low, xcen$x.high, type = 'interval2')
    Surv.Obj
    ## Then you can call survfit and plot the Kaplan-Meier curve:
    km <- survival::survfit(Surv.Obj ~ 1, conf.type = "none")
    km
    summary(km)
    plot(km, conf.int = FALSE, mark.time = TRUE)

    xcen$x <- apply(xcen, 1, max, na.rm=TRUE)
    xcen$type <- !mapply(identical, xcen$x.low, xcen$x.high)
    out <- NADA2::cfit(xcen$x, xcen$type)
    
    
    xcen1 <- data.frame(x.low = x, x.high = x)
    xcen <- rbind(xcen1, xcen)
    Surv.Obj <- survival::Surv(xcen$x.low, xcen$x.high, type = 'interval2')
    Surv.Obj
    NADA2::cenQQ(Surv.obj)
    
    
    xcen$x <- apply(xcen, 1, max, na.rm=TRUE)
    xcen$type <- ifelse(xcen$x.low == xcen$x.high, FALSE, TRUE)
    ## install.packages("EnvStats")
    library('EnvStats')
    NADA2::cenQQ(xcen$x, xcen$type, dist = "lnorm", Yname = 'plot name')
    

    ## #######################################  
    library('icenReg')

    set.seed(1)
    
    sim_data <- simIC_weib(n = 100, inspections = 5, inspectLength = 1)
    ph_fit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, 
                    data = sim_data)	
    ## Default fits a Cox-PH model
    
    summary(ph_fit)		
    ## Regression estimates close to true 0.5 and -0.5 values
    
    
    new_data <- data.frame(x1 = c(0,1), x2 = c(1, 1) )
    rownames(new_data) <- c('group 1', 'group 2')
    plot(ph_fit, new_data)
    ## plotting the estimated survival curves
    
    po_fit <- ic_sp(Surv(l, u, type = 'interval2') ~ x1 + x2, 
                    data = sim_data, model = 'po')
    ## fits a proportional odds model
    
    summary(po_fit)
    

    
    
    ## ####################################################################
    ## TEST CASE 6 - Johnson SU dataset
    ## #################################################################### 

    set.seed(1)
    jparms <- list(gamma=-3.3, delta=5.3, xi=1.8, lambda=1.9, type='SU')
    x <- ExtDist::rJohnsonSU(1000, param=jparms)
    
    ## add a few high points
    x <- c(x, 6, 7, 8)
    
    plotspace(2,2)
    hist_nwj(x, type='j', main='all known data')
    qqplot_nwj(x, type='j', mainadder='all known data')
    
    ## add censored points, redo fit, and remake plot
    xcen <- data.frame(x.low=0, x.high=c(3,4,5,6,7))
    out <- qqplot_censored(x, xcen, type='j', plotit='xq', plotcen=TRUE)

}



