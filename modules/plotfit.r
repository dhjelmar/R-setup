plotfit <- function(xx,
                    yy,
                    byvar      = NA,
                    xlabel     = NULL,
                    ylabel     = NULL,
                    bylabel    = NULL,
                    bynom      = NULL,
                    multifit   = FALSE,
                    grid       = TRUE,
                    color      = c('black', 'red', 'green3', 'cyan', 'blue',
                                   'darkorchid1', 'violet'),
                    pch        = 1,       # identify symbol style (1=open circle, 16=closed circle)
                    xlimspec   = NULL,
                    ylimspec   = NULL,
                    vlines     = NA,
                    main       = NULL,
                    equation   = TRUE,
                    interval   = 'conf',   # 'pred', 'mean', 'line', 'noline'
                    alpha      = 0.05,
                    sided      = 2,
                    bg         = "grey90",
                    legendloc  = NA,
                    legend.cex = 0.75,
                    outputfile = NULL,
                    suppress   = 'no') {

    ## usage: plotfit3(mtcars, 'mpg', 'disp', 'cyl')
    
    ## vlines   = vector where: 1st two entries correspond to 1st byvar
    ##                        : 2nd two entries correspond to 2nd byvar, etc.
    ## equation = TRUE writes equation as subtext under main title
    
    ## interval = 'conf' (default) plots mean fit and confidence limits
    ##          = 'pred' plots mean fit and prediction limits
    ##          = 'mean' plots mean fit only (and points)
    ##          = 'line' plots connect the dot lines
    ##          = 'noline' does not plot mean fit or limits (points only)
    
    ##-----------------------------------------------------------------------------
    ## pull label from the name of the vector if not specified otherwise
    if (is.null(xlabel)) xlabel   <- deparse(substitute(xx))
    if (is.null(ylabel)) ylabel   <- deparse(substitute(yy))
    if (is.null(bylabel)) bylabel <- deparse(substitute(byvar))

    ##-----------------------------------------------------------------------------
    # change xx, yy, and byvar to vector if not entered that way (more robust)
    xx <- unlist(xx)
    yy <- unlist(yy)
    byvar <- unlist(byvar)
    
    ##-----------------------------------------------------------------------------
    ## if bynom vector is specified, map byvar to closest bynom value for color and legend
    if (!is.null(bynom[1])) {

        if (length(bynom) == 1) {
            ## number of colors specified rather than a vector
            ## split range evenly into this number
            end1 <- min(byvar, na.rm=TRUE)
            end2 <- max(byvar, na.rm=TRUE)
            bynom <- seq(end1, end2, length.out = bynom)
        }

        
        ## find closest bynom value for each byvar value
        ## [1] used in function to select the 1st (lower) closest value if more than 1 found
        byvarnom <- unlist( purrr::pmap(list(x = byvar),
                                        function(x) DescTools::Closest(bynom, x)[1]) )
        if (length(byvarnom) > length(byvar)) {
            ## should no longer reach this error with addition of [1] in above function
            cat('\n##############################################\n')
            cat(  'FATAL ERROR: bynom results in one or more     \n')
            cat(  '             byvar values with multiple       \n')
            cat(  '             closesest points. Try again with \n')
            cat(  '             different bynom # or vector.)    \n')
            cat('\n##############################################\n\n')
            return()
        }

        ## replace byvar with byvarnom for plotting purposes
        byvar <- byvarnom

    }
    


    ##-----------------------------------------------------------------------------
    ## put data into dataframe and remove any pairs that do not have values
    xxcol  <- 1
    yycol  <- 2
    bycol  <- 3
    colcol <- 4
    pchcol <- 5
    if (is.na(byvar[1])) {
        df       <- na.omit( data.frame(xx, yy) )
        ## expand byvar=NA for every point
        df$byvar <- NA
    } else {
        df    <- na.omit( data.frame(xx, yy, byvar) )
        byvar <- df$byvar
    }
    ## write points back out
    xx <- df$xx
    yy <- df$yy
    
    ## sort the dataframe by byvar
    df <- df[gtools::mixedorder(df$byvar),]

    ## create color vector to identify color for every point
    cols  <- data.frame(byvar=unique(c(as.character(df$byvar))))
    if ((length(color) == 1) | (length(unique(byvar)) == 1)) {
        ## use same color for every point
        cols$color <- color[[1]]
    } else if (nrow(cols) <= length(color)) {
        ## user supplied color list has sufficient unique colors for plot
        cols$color <- color[1:nrow(cols)]
    } else {
        ## more colors needed than in user input so use color pallet
        pal   <- colorRampPalette(color)                  
        ncol  <- length(unique(byvar))
        color <- pal(ncol)
        ## to see colors
        ## scales::show_col(use)
        ## add column for color to dataframe
        cols$color <- color
    }
    ## merge with the original dataframe to add column for color
    df <- merge(df, cols, by="byvar")
    ## following was added to drop unwanted levels
    ## without it, the legend was wrong because it used the levels
    ## which had the colors in the wrong order
    df$color <- as.character(df$color) 
    color <- df$color
    
    ## create pch vector to identify symbol for every point
    symbols  <- data.frame(byvar=unique(c(as.character(df$byvar))))
    if ((length(pch) == 1) | (length(unique(byvar)) == 1)) {
      ## use same pch for every point
      symbols$pch <- pch[[1]]
    } else if (nrow(symbols) <= length(pch)) {
      ## input pch list has sufficient symbols for plot
      symbols$pch <- pch[1:nrow(symbols)]
    } else {
      ## more pch needed than defined
      npch  <- length(unique(byvar))
      if (npch < 26) {
          cat('WARNING: R only has 25 symbols but plot requires more than this.')
      }
      pch   <- seq(1, npch, 1)
      ## to see colors
      ## scales::show_col(use)
      ## add column for pch to dataframe
      symbols$pch    <- pch
    }
    ## merge with the original dataframe to add column for pch
    df <- merge(df, symbols, by="byvar")
    df$pch <- as.numeric(df$pch) # drops unwanted levels
    pch <- df$pch
    
    ## put order of columns back
    df <- data.frame(xx=df$xx, yy=df$yy, byvar=df$byvar, color=df$color, pch=df$pch)
    
    ## sort the dataframe by byvar
    df    <- df[gtools::mixedorder(df$byvar),]
    ## replace inputs with sorted order inputs
    xx    <- df$xx
    yy    <- df$yy
    byvar <- df$byvar

    ##-----------------------------------------------------------------------------
    ## handle byvar
    if (is.na(byvar[1])) {
        ## no byvar specified so single colored points
        colorpoints <- FALSE
        legendnames <- NA
        nfit <- 1
    } else if (isFALSE(multifit)) {
        ## byvar was specified so want colored points
        colorpoints <- TRUE
        ## but do not want multifit
        legendnames <- unique(byvar)
        nfit <- 1
    } else {
        ## byvar was specified so want colored points
        colorpoints <- TRUE
        ## want 1 fit for each unique byvar
        legendnames <- unique(byvar)
        nfit  <- length(legendnames)
    }

    
    ##-----------------------------------------------------------------------------
    ## define plot parameters
    if (is.null(xlimspec)) {
        xmax   <- max(xx, vlines, na.rm=TRUE)
        xmin   <- min(xx, vlines, na.rm=TRUE)
    } else {
        xmax   <- max(xlimspec)
        xmin   <- min(xlimspec)
    }
    
    if (isTRUE(colorpoints) & is.na(legendloc)) {
        ## make extra room for legend
        xmax_plot <- xmax + (xmax-xmin)*0.2
    } else {
        xmax_plot <- xmax
    }

    if (is.null(xlimspec)) xlimspec <- range(xmin, xmax_plot)
    if (is.null(ylimspec)) ylimspec <- range(ylimspec, yy, na.rm=TRUE)

    
    ##-----------------------------------------------------------------------------
    ## setup for output jpeg file
    if (!missing(outputfile)) jpeg(filename=outputfile)
    
    ## change background of plot to specified color
    par(bg=bg)

    ## create empty plot
    plot(xx, yy, type='n',
         xlim=xlimspec,  ylim=ylimspec,
         xlab=xlabel,  ylab=ylabel,
         main=main)

    if (isTRUE(grid)) grid(col='gray70')
    
    ##-----------------------------------------------------------------------------
    ## add points not used in fit
    ## not programmed

    ##-----------------------------------------------------------------------------
    ## add points, fits and vlines
    eq      <- NA
    pred    <- NA
    lineloc <- 1.5
    for (i in 1:nfit) {

        ## extract subset of data for ith byvar
        if (nfit == 1) {
            ## no byvar specified
            dfi <- df
        } else {
            ## byvar is specified
            dfi <- df[df[[bycol]] == legendnames[i],]
        }
        
        ## add points to plot
        points(df$xx, df$yy, col=df$color, pch=df$pch)
        #if (isTRUE(colorpoints) & nfit == 1) {
        #    ## color points even though there is only one fit
        #    points(xx, yy, col=color, pch=pch)
        #} else {
        #    ## add points to plot
        #    if (isTRUE(colorpoints)) {
        #        points(xx, yy, col=cols$color[i], pch=symbols$pch[i])
        #    } else {
        #        points(df$xx, df$yy, col=df$color, pch=df$pch)
        #    }
        #}

        ## add fit to plot
        ## fit is over range of data for current byvar
        ## lines extended to max(vlines, data)
        vpair <- c(vlines[i*2-1], vlines[i*2])
        if (interval != 'noline') {
            ## add fit to plot
            if (nfit ==1) {
                out   <- addfit(dfi[[xxcol]], dfi[[yycol]], col='black', vlines=vpair,
                                interval=interval, alpha=alpha, sided=sided, addpoints=FALSE)
            } else {
                out   <- addfit(dfi[[xxcol]], dfi[[yycol]], col=cols$color[i], vlines=vpair,
                                interval=interval, alpha=alpha, sided=sided, addpoints=FALSE)
            }
            eq[i] <- out$equation
            pred[i] <- list(out$pred)
            
            if (i == 1) {
                if (nfit == 1) {
                    fitname <- NA
                } else {
                    fitname <- legendnames[i]
                }
                fits <- data.frame(fit       = fitname,
                                   slope     = out$slope,
                                   intercept = out$intercept,
                                   rise      = out$rise)
            } else {
                fits  <- rbind(fits, data.frame(fit       = legendnames[i],
                                                slope     = out$slope,
                                                intercept = out$intercept,
                                                rise      = out$rise))
            }
        }
            
        ## add vertical lines if specified
        if ( !is.nothing(vlines) ) abline(v=vpair, lty=3, col=cols$color[i])
        
        
        if (isTRUE(equation) & interval != 'line') {
            ## add the fit equations under main title
            ## subtitle <- list(eq1, eq2)
            ## mtext(subtitle, side=3, line=c(0.75, 0), cex=.75, col=color)
            lineloc = lineloc - 0.75
            if (nfit ==1) {
                mtext(eq[i], side=3, line=lineloc, cex=0.75, col='black')
            } else {
                mtext(eq[i], side=3, line=lineloc, cex=0.75, col=cols$color[i])
            }
        }
    }

    ##-----------------------------------------------------------------------------
    ## add legend
    ## determine legend location by slope of all data
    if (isTRUE(colorpoints)) {
        if (is.na(legendloc)) {
            if (diff(range(xx)) > 0) {
                ## xx has more than 1 value so possible to define slope
                fitall       <- lm(yy~xx, data=df)
                fitall_slope <- fitall$coefficients[[2]]
                if (fitall_slope > 0) {
                    legendloc <- 'bottomright'
                } else {
                    legendloc <- 'topright'
                }
            } else {
                ## xx hsa only 1 value so not possible to define slope
                ## put legend at bottom right
                legendloc <- 'bottomright'
            }
        }
        legend(legendloc, title=bylabel, col=cols$color, legend=legendnames, 
               pch=symbols$pch, cex=legend.cex)
    }

    ## for output jpeg file
    if (!missing(outputfile)) dev.off()

    ## change background back to base R default
    par(bg="white")
    
    ##-----------------------------------------------------------------------------
    if (suppress == 'no') {
        ## return fitted equations    
        if (nfit > 1) {
            ## more than one equation so return as a list
            eq <- as.list(eq)
            names(eq) <- legendnames
        }
        legendprint <- data.frame(legendnames, color = cols$color, symbols = symbols$pch)
        if (interval == 'noline') {
            return(list(legend = legendprint))
        } else {
            return(list(fits = fits, legend = legendprint, pred=pred))
        }
    }

}

testplots <- function() {
    os <- .Platform$OS.type
    if (os == "windows") {
        #source('D:/Documents/01_Dave/Programs/GitHub_home/R-setup/setup.r')
        source('D:/Documents/01_Dave/Programs/GitHub_home/R-setup/modules/plotspace.r')
        source('D:/Documents/01_Dave/Programs/GitHub_home/R-setup/modules/addfit.r')
        source('D:/Documents/01_Dave/Programs/GitHub_home/R-setup/modules/is.nothing.r')
    } else {
        ## source('/home/dlhjel/GitHub_repos/R-setup/setup.r')
        source('~/Documents/GitHub/R-setup/modules/plotspace.r')
        source('~/Documents/GitHub/R-setup/modules/addfit.r')
        source('~/Documents/GitHub/R-setup/modules/is.nothing.r')
    }
    df <- mtcars
    plotspace(1,2)
    plotfit(df$hp, df$mpg, main='test 1: 1 fit; 1 color; 1 symbol')
    with(df, plotfit(hp, mpg, cyl, main='test 2: 1 fit with byvar for color'))
    
    plotspace(1,2)
    plotfit(df$hp, df$mpg, df$cyl, multifit=TRUE, main='test 4 multifit; no vlines')
    plotfit(df$hp,         df$mpg,         df$cyl, 
            'horsepower', 'miles per gal', 'cylinders',
            multifit=TRUE,
            vlines=c(50,120,  NA,NA,   90,410),
            xlimspec=c(0,500),
            main='test 5 multifit with vlines for 2 of 3 fits')
    
    df1 <- df
    df1$type <- 'type1'
    df2     <- df1
    df2$hp  <- df1$hp + 100
    df2$mpg <- df1$mpg * 1.1 + rnorm(1, 10, 1)
    df2$type <- 'type2'
    df <- rbind(df1, df2)
    
    plotspace(1,2)    
    plotfit(df$hp, df$mpg, df$type, main='test 6: 1 fit; color points with byvar')
    plotfit(df$hp, df$mpg, df$type, multifit=TRUE, vlines=c(30,350,  100,450), 
            main='test 7: multifit with vlines')
    
    df[58:64,]$type <- 'type is 3'
    plotfit(df$hp, df$mpg, df$type, main='test 8: 1 fit; color points with byvar')
    out <- plotfit(df$hp, df$mpg, df$type, multifit=TRUE, main='test 9: multifit')

    plotspace(1,3)
    plotfit(df$mpg, df$disp, df$drat, main='test 10: 1 fit; byvar')
    plotfit(df$mpg, df$disp, df$drat, bynom=5, main='test 11: auto 5 bynom')
    plotfit(df$mpg, df$disp, df$drat, bynom=c(2.76, 3.302, 3.845, 4.387, 4.93), main='test 12: specified bynom')

    plotspace(1,3)
    plotfit(mtcars$wt, mtcars$mpg, mtcars$qsec, main='test 13: 1 fit; byvar')
    with(mtcars, plotfit(wt, mpg, qsec, bynom=seq(14, 24, 2), main='test 14: sequence bynom'))
    with(mtcars, plotfit(wt, mpg, qsec, bynom=seq(14, 24, 2),
                         color=c('black', 'darkviolet', 'blue', 'green', 'red'),
                         main='test 15: spec color'))
    
    plotspace(1,3)
    plotfit(mtcars$wt, mtcars$mpg, mtcars$qsec, pch=c(16,11,4), main='test 16: 1 fit; auto color + pch')   # should overwrite pch with many different symbols
    with(mtcars, plotfit(wt, mpg, qsec, bynom=seq(14, 24, 2), pch=c(14,24,2), main='test 17: 1 fit; spec pch'))    # again, overwrite
    with(mtcars, plotfit(wt, mpg, qsec, bynom=seq(14, 24, 2), pch=c(0, 0, 16, 2, 3, 4, 5, 6),
                         color=c('black', 'darkviolet', 'blue', 'darkgreen', 'red'),
                         main='test 18: 1 fit; spec pch + color'))          # uses specificed pch because number >= unique(byvar)
}
#testplots()
#source('D:/Documents/01_Dave/Programs/GitHub_home/R-setup/modules/addfit.r')
#df <- mtcars
#plotspace(1,3)
#plotspace(1,2)

