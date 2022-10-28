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
        
    if (length(color) == length(xx)) {
        ## color provided for each datapoint
        ## create cols dataframe with byvar and color columns
        df$color <- color
        cols     <- select(df, byvar, color) 
    
    } else {
        ## color palette provided rather than defining one color for each point
        ## sort the dataframe by byvar
        df <- df[order(df$byvar),]
        ## expand color vector to identify color for every point
        cols  <- data.frame(byvar=unique(c(as.character(df$byvar))))
        if (nrow(cols) <= length(color)) {
            ## palette has sufficient unique colors for plot
            cols$color <- color[1:nrow(cols)]
        } else {
            ## too many colors for palette so ramp
            ## define color palette
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
        ## put order of columns back
        df <- data.frame(xx=df$xx, yy=df$yy, byvar=df$byvar, color=df$color)
    }
    ## sort the dataframe by byvar
    df    <- df[order(df$byvar),]
    ## replace inputs with sorted order inputs
    xx    <- df$xx
    yy    <- df$yy
    byvar <- df$byvar
    ## following was added to drop unwanted levels
    ## without it, the legend was wrong because it used the levels
    ## which had the colors in the wrong order
    df$color <- as.character(df$color) 
    color <- df$color
    
    
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
    xmax   <- max(xlimspec, xx, xx, vlines, na.rm=TRUE)
    xmin   <- min(xlimspec, xx, xx, vlines, na.rm=TRUE)
    ylimspec <- range(ylimspec, yy, yy, na.rm=TRUE)
    if (isTRUE(colorpoints) & is.na(legendloc)) {
        ## make extra room for legend
        xmax_plot <- xmax + (xmax-xmin)*0.2
    } else {
        xmax_plot <- xmax
    }
    xlimspec <- range(xmin, xmax_plot)
        
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
        
        ## fit is over range of data for current byvar
        ## lines extended to max(vlines, data)
        vpair <- c(vlines[i*2-1], vlines[i*2])
        if (interval != 'noline') {
            ## add points and fit to plot
            out   <- addfit(dfi[[xxcol]], dfi[[yycol]], col=cols$color[i], pch=pch, vlines=vpair,
                            interval=interval, alpha=alpha, sided=sided)
            eq[i] <- out$equation
            
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
        } else {
            ## only add points to plot
            points(xx, yy, col=color, pch=pch)
        }
            
        ## color points even though there is only one fit
        ## browser()
        if (isTRUE(colorpoints) & nfit == 1) {
            points(xx, yy, col=color, pch=pch)
        }

        if (isTRUE(equation) & interval != 'line') {
            ## add the fit equations under main title
            ## subtitle <- list(eq1, eq2)
            ## mtext(subtitle, side=3, line=c(0.75, 0), cex=.75, col=color)
            lineloc = lineloc - 0.75
            mtext(eq[i], side=3, line=lineloc, cex=0.75, col=cols$color[i])
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
               pch=pch, cex=legend.cex)
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
        legendprint <- data.frame(legendnames, color = cols$color)
        if (interval == 'nolne') {
            return(list(legend = legendprint))
        } else {
            return(list(fits = fits, legend = legendprint))
        }
    }

}

testplots <- function() {
    source('/home/dlhjel/GitHub_repos/R-setup/setup.r')
    df <- mtcars
    plotspace(1,3)
    plotfit(df$hp, df$mpg)
    plotfit(df$hp, df$mpg, df$cyl)
    with(df, plotfit(hp, mpg, cyl))
    
    plotspace(1,2)
    plotfit(df$hp, df$mpg, df$cyl, multifit=TRUE)
    plotfit(df$hp,         df$mpg,         df$cyl, 
            'horsepower', 'miles per gal', 'cylinders',
            multifit=TRUE,
            vlines=c(50,120,  NA,NA,   100,350),
            xlimspec=c(0,500),
            main='my title')
    
    df1 <- df
    df1$type <- 'type1'
    df2     <- df1
    df2$hp  <- df1$hp + 100
    df2$mpg <- df1$mpg * 1.1 + rnorm(1, 10, 1)
    df2$type <- 'type2'
    df <- rbind(df1, df2)
    
    plotspace(1,2)    
    plotfit(df$hp, df$mpg, df$type)
    plotfit(df$hp, df$mpg, df$type, multifit=TRUE, vlines=c(30,350,  100,450))
    
    df[58:64,]$type <- 'type is 3'
    plotfit(df$hp, df$mpg, df$type)
    out <- plotfit(df$hp, df$mpg, df$type, multifit=TRUE)

    plotspace(1,3)
    plotfit(df$mpg, df$disp, df$drat)
    plotfit(df$mpg, df$disp, df$drat, bynom=5)
    plotfit(df$mpg, df$disp, df$drat, bynom=c(2.76, 3.302, 3.845, 4.387, 4.93))

    plotspace(1,3)
    plotfit(mtcars$wt, mtcars$mpg, mtcars$qsec)
    with(mtcars, plotfit(wt, mpg, qsec, bynom=seq(14, 24, 2)))
    with(mtcars, plotfit(wt, mpg, qsec, bynom=seq(14, 24, 2),
                         color=c('black', 'darkviolet', 'blue', 'green', 'red')))
}
