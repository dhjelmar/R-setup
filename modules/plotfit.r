plotfit <- function(df, xx, yy,
                     byvar      = NULL,
                     xlimspec   = NULL,
                     ylimspec   = NULL,
                     vlines     = NA,
                     main       = NULL,
                     equation   = TRUE,
                     xlabel     = NULL,
                     ylabel     = NULL,
                     interval   = 'conf', alpha=0.05, sided=2,
                     bg         = "grey90",
                     color      = palette(),
                     outputfile = NULL,
                     suppress   = 'no') {

    ## usage: plotfit3(mtcars, 'mpg', 'disp', 'cyl')
    
    ## vlines   = vector where: 1st two entries correspond to 1st byvar
    ##                        : 2nd two entries correspond to 2nd byvar, etc.
    ## equation = TRUE writes equation as subtext under main title
    
    ##-----------------------------------------------------------------------------
    ## extract column locations and values from dataframe
    if (is.null(xlabel)) xxlabel <- xx      ## name of xx variable
    xxcol <- which(grepl(xx, names(df)))    ## xx column
    xx    <- df[, xxcol]                    ## xx values
    if (is.null(ylabel)) yylabel <- yy
    yycol <- which(grepl(yy, names(df)))  
    yy    <- df[, yycol]
    if (is.null(byvar)) {
        ## no byvar specified
        legendnames <- NULL
        nfit <- 1
    } else {
        bycol <- which(grepl(byvar, names(df)))  
        by    <- df[, bycol]
        legendnames <- unique(df[[bycol]])
        nfit  <- length(legendnames)
    }

    ##-----------------------------------------------------------------------------
    ## define plot parameters
    xmax   <- max(xlimspec, xx, xx, vlines, na.rm=TRUE)
    xmin   <- min(xlimspec, xx, xx, vlines, na.rm=TRUE)
    ylimspec <- range(ylimspec, yy, yy, na.rm=TRUE)
    if (nfit > 1) {
        ## more than 1 fit so need a legend
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
    ## determine number of unique fits needed

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
        eq[i] <- addfit(dfi[[xxcol]], dfi[[yycol]], col=color[i], vlines=vpair,
                        interval=interval, alpha=alpha, sided=sided)

        ## add the fit equations under main title
        ## subtitle <- list(eq1, eq2)
        ## mtext(subtitle, side=3, line=c(0.75, 0), cex=.75, col=color)
        lineloc = lineloc - 0.75
        mtext(eq[i], side=3, line=lineloc, cex=0.75, col=color[i])
    }

    ##-----------------------------------------------------------------------------
    ## add legend
    ## determine legend location by slope of all data
    if (nfit > 1) {
        fitall       <- lm(yy~xx, data=df)
        fitall_slope <- fitall$coefficients[[2]]
        if (fitall_slope > 0) {
            legendloc <- 'bottomright'
        } else {
            legendloc <- 'topright'
        }
        legend(legendloc, title=byvar, col=color, legend=legendnames, pch=1)
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
        return(eq)
    }

}

testplots <- function() {
    source('/home/dlhjel/GitHub_repos/R-setup/setup.r')
    plotfit(mtcars, 'cyl', 'mpg')
    plotfit(mtcars, 'cyl', 'mpg', vlines=c(5,7))
    plotfit(mtcars, 'cyl', 'mpg', vlines=c(3.9,8.1))
    plotfit(mtcars, 'cyl', 'mpg', vlines=c(2,9), xlimspec=c(0, 10))
    plotfit(mtcars, 'cyl', 'mpg', outputfile='junk.jpg')
    plotfitcol(mtcars, 'cyl', 'mpg', byvar='cyl', ncol=3)
    plotfitcold(mtcars, 'cyl', 'mpg', byvar=cyl)
    plotfitcold(mtcars, 'cyl', 'mpg', byvar=cyl,                xlimspec=c(0,10))
    plotfitcold(mtcars, 'cyl', 'mpg', byvar=cyl, vlines=c(2,9), xlimspec=c(0,10))
    plotfitcold(mtcars, 'cyl', 'mpg', byvar=cyl, vlines=c(5,7))

    df1 <- mtcars
    df1$type <- 'type1'
    df2     <- df1
    df2$mpg <- df1$mpg * 1.1 + rnorm(1, 10, 1)
    df2$type <- 'type2'
    df <- rbind(df1, df2)
    df[58:64,]$type <- 'type is 3'

    plotfit(df, 'mpg', 'disp')

    plotfit(df, 'mpg', 'disp', 'type')

    out <- plotfit(df, 'mpg', 'disp', 'type',
                   xlimspec=c(0, 50), ylimspec=c(0, 500),
                   vlines=c(10,35,    NA,NA,   25,50),
                   color = c('black', 'green', 'red'),
                   main='specified black/green/red and no vlines for green fit')

    plotspace(1,2)    
    plotfitcol(df, mpg, disp, cyl, ncol=3)
    plotfit(df, 'mpg', 'disp', 'cyl')
}
