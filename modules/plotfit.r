plotfit <- function(xx,
                    yy,
                    byvar      = NA,
                    xlabel     = NULL,
                    ylabel     = NULL,
                    bylabel    = NULL,
                    multifit   = FALSE,
                    color      = palette(),
                    xlimspec   = NULL,
                    ylimspec   = NULL,
                    vlines     = NA,
                    main       = NULL,
                    equation   = TRUE,
                    interval   = 'conf', alpha=0.05, sided=2,
                    bg         = "grey90",
                    outputfile = NULL,
                    suppress   = 'no') {

    ## usage: plotfit3(mtcars, 'mpg', 'disp', 'cyl')
    
    ## vlines   = vector where: 1st two entries correspond to 1st byvar
    ##                        : 2nd two entries correspond to 2nd byvar, etc.
    ## equation = TRUE writes equation as subtext under main title
    ## df       = if provided, then xx, yy, byvar, and color are in dataframe
    ##            and need to specify variables in quotes when call function

    
    ##-----------------------------------------------------------------------------
    ## pull label from the name of the vector if not specified otherwise
    if (is.null(xlabel)) xlabel   <- deparse(substitute(xx))
    if (is.null(ylabel)) ylabel   <- deparse(substitute(yy))
    if (is.null(bylabel)) bylabel <- deparse(substitute(byvar))

    
    ##-----------------------------------------------------------------------------
    ## put data into dataframe
    xxcol  <- 1
    yycol  <- 2
    bycol  <- 3
    colcol <- 4
    df <- data.frame(xx, yy, byvar)
    
    if (length(color) == length(xx)) {
        ## color provided for each datapoint
        df$color <- color
    
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
            ## define color pallette
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
        ## pull out new color vector
        color <- df$color
    }
    ## sort the dataframe by byvar
    df    <- df[order(df$byvar),]
    ## replace inputs with sorted order inputs
    xx    <- df$xx
    yy    <- df$yy
    byvar <- df$byvar
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
    if (isTRUE(colorpoints)) {
        ## make extra room for legend
        xmax_plot <- xmax + (xmax-xmin)*0.2
    } else {
        xmax_plot <- xmax
    }
    xlimspec <- range(xmin, xmax_plot)
        
    ##-----------------------------------------------------------------------------
    ## setup for output jpeg file
    if (!missing(outputfile)) jpeg(filename=outputfile)
    
    ## create empty plot
    plot(xx, yy, type='n',
         xlim=xlimspec,  ylim=ylimspec,
         xlab=xlabel,  ylab=ylabel,
         main=main)

    ## change background of plot to specified color
    par(bg=bg)  

    ##-----------------------------------------------------------------------------
    ## add nofit points
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
        out   <- addfit(dfi[[xxcol]], dfi[[yycol]], col=cols$color[i], vlines=vpair,
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
        
        ## color points even though there is only one fit
        if (isTRUE(colorpoints) & nfit == 1) {
            points(xx, yy, col=color)
        }
        
        ## add the fit equations under main title
        ## subtitle <- list(eq1, eq2)
        ## mtext(subtitle, side=3, line=c(0.75, 0), cex=.75, col=color)
        lineloc = lineloc - 0.75
        mtext(eq[i], side=3, line=lineloc, cex=0.75, col=cols$color[i])
    }

    ##-----------------------------------------------------------------------------
    ## add legend
    ## determine legend location by slope of all data
    if (isTRUE(colorpoints)) {
        fitall       <- lm(yy~xx, data=df)
        fitall_slope <- fitall$coefficients[[2]]
        if (fitall_slope > 0) {
            legendloc <- 'bottomright'
        } else {
            legendloc <- 'topright'
        }
        legend(legendloc, title=bylabel, col=cols$color, legend=legendnames, pch=1)
    }

    ## for output jpeg file
    if (!missing(outputfile)) dev.off()

    ##-----------------------------------------------------------------------------
    if (suppress == 'no') {
        ## return fitted equations    
        if (nfit > 1) {
            ## more than one equation so return as a list
            eq <- as.list(eq)
            names(eq) <- legendnames
        }
        return(fits)
    }

}

testplots <- function() {
    source('/home/dlhjel/GitHub_repos/R-setup/setup.r')
    df <- mtcars
    plotfit(df$hp, df$mpg)
    plotfit(df$hp, df$mpg, df$cyl)
    with(df, plotfit(hp, mpg, cyl))
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

    plotfitcol(df, mpg, disp, drat, ncol=5)
    plotfit(df$mpg, df$disp, df$drat)
}
