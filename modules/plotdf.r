plotdf <- function(df, xx, yy,
                   byvar      = NA,  
                   by.as.text = FALSE,
                   ncol       = NA,
                   nofit      = FALSE,
                   multifit   = FALSE,
                   vlines     = NA,
                   xlimspec   = NA,
                   ylimspec   = NA,
                   main       = NULL,
                   printeq    = 'yes',
                   xlabel     = NA,
                   ylabel     = NA,
                   interval   = 'conf',
                   alpha      = 0.05,
                   sided      = 2,
                   cols       = c('black', 
                                  '#9900CC', 'violet', 
                                  'blue', '#6699FF', '#33FFF6', 
                                  'green3', '#33FF66', 
                                  'orange', '#FF6600', 
                                  'red', 'firebrick'),
                   color      = NULL,
                   bg         = 'white',   # "grey90",  # grey90 resulted in hiding data
                   suppress   = 'yes',
                   outputfile = NA) {
    ## create scatter plot with linear regression line and 95% confidence lines
    ## xx and yy are variables within dataframe df

    ## required inputs: df = dataframe
    ##                  xx = x-axis parameter
    ##                  yy = y-axis parameter

    ## options (defaults shown in function definition above):
    ##    byvar    = variable to use to color code datapoints (NA is default)
    ##    by.as.text = FALSE to treat byvar values as continuous numerical values (default but requires byvar to be numeric)
    ##               = TRUE to treat byvar values as text (default if byvar is not numeric)
    ##                 (therefore, by.as.text is only needed if want to treat numeric byvar as a non-numeric label)
    ##    ncol     = NA sets ncol to the number of unique byvar values (default)
    ##             = # to specify number of colors (only used for numeric data)
    ##    multifit = FALSE specifies single fit to all data (default)
    ##             = TRUE  specifies separate fit for each byvar
    ##    vlines   = vector of 2 vertical lines to be added (e.g., c(0.5, 2.0))
    ##               (fits will extend to extent of vlines)
    ##    xlimspec = vector defining xlimits of plot (e.g., c(-2, 11.5))
    ##    ylimspec = vector defining ylimits of plot
    ##    main     = title for plot
    ##    printeq  = use the fit equation as subtext on the title (max of 2 equations)
    ##    xlabel   = NA = x-variable name used for x-axis
    ##             = user specified for label for x-axis
    ##    ylabel   = NA = y-variable name used for y-axis
    ##             = user specified for label for y-axis
    ##    inverval = 'conf' adds confidence limits to the plot
    ##    alpha    = 1 - confidence
    ##    sided    = 2 = 2-sided confidence limits
    ##             = 1 = 1-sided confidence limits (but both sides plotted)
    ##    color    = NULL uses cols to automatically set colors
    ##             = column of colors to use for data and, if multifit=TRUE, fits (possible future option)
    ##    cols     = vector of colors for use in automatically setting colors
    ##    bg       = background color for plot
    ##    suppress =
    ##    outputfile = name for PDF to be created
    ##               = NA does not create a PDF (only output is to screen)
    
    ## example usage: plotdf(mtcars, disp, mpg, byvar=cyl)

    ##-----------------------------------------------------------------------------
    ## PREPARE DATAFRAME AND EXTRACT LABELS FOR X, Y, AND BYVAR
    
    ## extract names of xx, yy, and byvar if not passed in with quotes
    xlabel  <- deparse(substitute(xx))
    ylabel  <- deparse(substitute(yy))
    bylabel <- deparse(substitute(byvar))

    ## attempt to replace xx, yy, and byvar with vector values
    xx    <- eval(substitute(xx),df)
    yy    <- eval(substitute(yy),df)    
    byvar <- eval(substitute(byvar),df) 

    ## if xx and yy were passed in without quotes, xx and yy will be vectors to be plotted
    ## if xx and yy were passed in with quotes, xx and yy will be name of vector to be plotted
    if (typeof(xx) == 'character') {
        xxcol <- which(grepl(xx, names(df)))  
        xx    <- df[, xxcol]
    }
    if (typeof(yy) == 'character')  {
        yycol <- which(grepl(yy, names(df)))  
        yy    <- df[, yycol]
    }
    if (typeof(byvar) == 'character' & length(byvar)==1)  {
        if (!is.na(byvar)) {
            ## byvar was specified
            byvar <- which(grepl(byvar, names(df)))  
            byvarv <- df[, byvar]
        } else {
            ## no byvar was specified
            ## create vector of values for byvar equal in length to xx and yy
            byvar <- rep(NA, length(xx))
        }
    }

    ## put xx, yy, and byvar into dataframe
    newdf <- data.frame(xx, yy, byvar=byvar)
    names(newdf) <- c('xx', 'yy', 'byvar')

    
    ##-----------------------------------------------------------------------------
    ## DETERMINE X AND Y RANGE FOR NEEDED FOR DATA AND REQUESTED LINES
    if (is.null(vlines[1])) {
        ## vlines not specified
        xmin <- min(newdf$xx, vlines, na.rm=TRUE)
        xmax <- max(newdf$xx, vlines, na.rm=TRUE)
    } else {
        ## vlines specified so use thoe for xmin and xmax
        xmin <- min(vlines, na.rm=TRUE)   # dlh may not be good ifused for fit and axes?
        xmax <- max(vlines, na.rm=TRUE)
    }
    ymin <- min(newdf$yy, na.rm=TRUE)
    ymax <- max(newdf$yy, na.rm=TRUE)


    ##-----------------------------------------------------------------------------
    ## DETERMINE NUMBER OF DATAPOINT COLORS AND FITS NEEDED
    if (is.na(byvar[1])) {
        ## no byvar specified so single color and single fit group
        ncol      <- 1
        if (isTRUE(nofit)) {
            ## no fit is to be put on the plot
            fitgroups <- 0
        } else {
            fitgroups <- 1
        }
        
    } else {
        ## byvar variable used to determine datapoint colors
        if (is.na(ncol)) {
            ## number of colors not specified
            ncol <- length(unique(newdf$byvar))
        }

        if (isTRUE(nofit)) {
            ## no fit is to be put on the plot
            fitgroups <- 0
        } else if (isFALSE(multifit)) {
            ## default is a single fit group
            fitgroups <- 1
        } else {
            ## if anything other than NA is specified, then separate fit per color
            fitgroups <- ncol
        }
    }    

    ##-----------------------------------------------------------------------------
    ## DEFINE COLOR PALLETTE FOR PLOT
    ## first determine whether by.as.text is applicable if selected
    if (isFALSE(by.as.text) & !is.numeric(newdf$byvar)) by.as.text <- TRUE
    ## define color pallette
    pal  <- colorRampPalette(cols)                  
    if (ncol == 1) {
        newdf$color <- pal(ncol)   # first color in cols
    } else {
        if (isFALSE(by.as.text)) {
            ## continuous variable for color
            newdf$color <- pal(ncol)[as.numeric(cut(newdf$byvar,breaks=ncol))]
            
            ## set colors and labels for legend
            legendcolor <- pal(ncol)
            low  <- min(newdf$byvar,na.rm=TRUE)
            high <- max(newdf$byvar,na.rm=TRUE)
            legendlabel <- seq( low, high, (high-low)/(ncol-1) )
            legendlabel <- signif(legendlabel,digits=4)
    
        } else {
            ## discrete variable (e.g., label) for color
            cols  <- data.frame(byvar=unique(c(as.character(newdf$byvar))))
            ncol  <- nrow(cols)
            color <- rainbow(ncol)
            color <- pal(ncol)
            ## to see colors
            ## scales::show_col(use)
            cols$color <- color
            newdf <- merge(newdf, cols, by="byvar")

            ## set colors and labels for legend
            legendcolor <- unique(newdf$color)
            legendlabel <- unique(newdf$byvar)            
        }
    }

    ##-----------------------------------------------------------------------------
    ## SPLIT DATAFRAME BY BYVAR IF MULTIPLE FITS ARE NEEDED
    if (fitgroups <= 1) {
        ## put newdf into list so can handle the same as in more complicated case for ncol > 1
        newdfl <- list(newdf)
    } else {
        ncol <- length(unique(byvar))
        newdfl    <- split(newdf, f = newdf$byvar)
        ## rename newdf list items to be sequential
        names(newdfl) <- seq(1:length(newdfl))
    }


    ##-----------------------------------------------------------------------------
    ## ADJUST X AND Y RANGE FOR PLOT IF NEEDED 
    if (ncol == 1) {
        ## do not need extra room for legend since no colors on plot needing legend
        xmax_plot <- xmax
    } else {
        ## make room for legend
        xmax_plot <- xmax + (xmax-xmin)*0.2  
    }
    if (missing(xlimspec)) { xlimspec <- c(xmin,xmax_plot) }
    if (missing(ylimspec)) { ylimspec <- c(ymin,ymax) }


    ##-----------------------------------------------------------------------------
    ## CREATE PLOT SPACE AND PLOT POINTS

    ## setup for output jpeg file
    if (!missing(outputfile)) jpeg(filename=outputfile)

    ## draw plot without points
    plot(newdf$xx,
         newdf$yy,
         col=newdf$col,
         ## type='n',
         xlim=xlimspec, ylim=ylimspec,
         xlab=xlabel, ylab=ylabel,
         main=main)
    
    ## changes background of plot to specified color
    par(bg=bg)

    ## add grid
    grid(col='gray70')

 ##    ## add points
##     for (i in 1:fitgroups) {
##         points(newdfl[[i]]$xx, newdfl[[i]]$yy, col=newdfl[[i]]$color)
##     }

    
    ## if request upper/lower range bounds, add to plot
    if ( !is.na(vlines[1]) ) {
        abline(v=vlines[1], col='black', lty='dashed')
        abline(v=vlines[2], col='black', lty='dashed')
    }        
    


    ##-----------------------------------------------------------------------------
    ## PERFORM REGRESSION AND PLOT fit FOR EACH FITGROUP

    if (fitgroups > 0) {

        fit <- NA
        fitcol <- NA
        intercept <- NA
        slope <- NA
        rise <- NA
        equation <- NA
        legendn <- NA
        legendc <- NA
        
        for (ifit in 1:fitgroups) {
            fit             <- lm(yy~xx, data=newdfl[[ifit]])
            ## print results of fit
            estbound(fit)
            intercept[ifit] <- fit$coefficients[[1]]
            slope[ifit]     <- fit$coefficients[[2]]
            ## calculate rise of fit over range bounds if supplied,
            ## otherwise calculate rise over range of data
            if ( !is.na(vlines[1]) ) {
                rise[ifit] <- (vlines[2] - vlines[1]) * slope[ifit]
            } else {
                rise[ifit] <- (xmax - xmin) * slope[ifit]
            }        
            equation[ifit] = paste0("y = ", signif(slope[ifit],4), " * x + ",
                                    signif(intercept[ifit],4),
                                    ", ", expression(Delta), " = ", signif(rise[ifit],4))

            ## legend info
            if (ncol == 1) {
                ## no legend
                legendn <- NULL
                legendc <- NULL
            } else if (fitgroups == 1 & ncol > 1) {
                ## all legend names and colors are defined in this newdfl[[1]] 
                legendn <- list( unique(newdfl[[ifit]]$byvar) )
                legendc <- list( unique(newdfl[[ifit]]$color) )
            } else if (fitgroups > 1 & ncol > 1) {
                ## byvar has a single color and subsequent ifit numbers will 
                legendn[ifit] <- newdfl[[ifit]]$byvar[1]
                legendc[ifit] <- newdfl[[ifit]]$color[1]
            }

            ## add fit to plot
            new.xx <- seq(min(newdfl[[ifit]]$xx, vlines, na.rm=TRUE),
                          max(newdfl[[ifit]]$xx, vlines, na.rm=TRUE), len=100)
            pred   <- predict(fit, new=data.frame(xx=new.xx), interval=interval, level=1-0.05/sided)
            if (fitgroups == 1) {
                fitcol[ifit] <- 'black'
            } else {
                ## colorcode lines
                fitcol[ifit] <- newdfl[[ifit]]$color[1]
            }
            lines(new.xx, pred[,"fit"], lwd=2, col=fitcol[ifit])

            ## if request confidence bounds, add to plot
            if (interval != 'none') {
                lines(new.xx, pred[,"lwr"], lty=3, col=fitcol[ifit])
                lines(new.xx, pred[,"upr"], lty=3, col=fitcol[ifit])
            }


        }
    }
    
    ##-------------------------------------------------------------------------------------     
    ## Add title and subtitle with equation if requested
    if (printeq == 'yes') {
        if (fitgroups == 1) {
            ## add fit equation under main title for the only fit
            mtext(equation[1], side=3, line=0, cex=1)
        } else if (fitgroups == 2) {
            ## add the 2 fit equations under main title
#            subtitle <- list('line color: red = normal, blue = Weibull, black = Johnson',
#                             'line type: solid = distribution or mean, dashed = 1-sided upper bound')
            subtitle <- list(equation[1], equation[2])
            ## following works but is a little clunky working with lists when do not need to
            ## textcol  <- c(newdfl[[1]]$color[1], newdfl[[2]]$color[1])
            textcol  <- unique(newdf$color)
            mtext(subtitle, side=3, line=c(0.75, 0), cex=.75, col=textcol)
        } else{
            ## too many equations put them all in the title so skip
        }
    }
    

    ##-------------------------------------------------------------------------------------     
    ## add legend
    if (ncol > 1) {
        ## no legend is needed if there is only one color on the plot
        
        ## determine legend location by slope of all data
        fitall       <- lm(yy~xx, data=newdf)
        fitall_slope <- fitall$coefficients[[2]]
        if (fitall_slope > 0) {
            legendloc <- 'bottomright'
        } else {
            legendloc <- 'topright'
        }
        legend(legendloc, title=bylabel, col=legendcolor, pch=19,
               legend=legendlabel)
    }

    
    ##-----------------------------------------------------------------------------
    ## WRAPUP AND RETURN

    ## create output file if requested
    if (!missing(outputfile)) dev.off()
    
    if (suppress == 'no')
        return( list(intercept = intercept, slope = slope, rise = rise, pred = as_tibble(pred)) )
}




plotdf_test <- function() {
    ## example: plotdf(df,gnom,average_g_ts,TS_ID,bg="white")
    df <- mtcars
    df$size[df$cyl <= 4] <- 'small'
    df$size[df$cyl >  4] <- 'large'
    df$test <- df$cyl
    df$test[1] <- 6.01
##     plotdf(df, disp, mpg,                               main='test 1: plotfit')
##     plotdf(df, disp, mpg, byvar=wt,    ncol=2,          main='test 2: plotfitcol, 2 colors')
##     plotdf(df, disp, mpg, byvar=wt,    ncol=3,          main='test 3: plotfitcol, 3 colors')
##     plotdf(df, disp, mpg, byvar=test,  by.as.text=TRUE, main='test 4a: plotfitcold; want to 4, 6, 6.01, and 8 as text')
##     plotdf(df, disp, mpg, byvar=wt,                     main='test 4: plotfitcold too many numeric')
##     plotdf(df, disp, mpg, byvar=wt,    ncol=8,          main='test 5: plotfitcol  reduce numeric')
##     plotdf(df, disp, mpg, byvar=size,                   main='test 6: plotfitcold few text')
##     plotdf(df, disp, mpg, byvar=size,                   main='test 7: plotfitcold multiple fits')
##     plotdf(df, disp, mpg, byvar=cyl,  by.as.text=FALSE, multifit=TRUE, main='test 8: plotfitcol  multiple fits')
##     plotdf(df, disp, mpg, byvar=cyl, by.as.text=FALSE, multifit=TRUE,
##            vlines=c(50, 500),
##            main='test 8a: plotfitcol  multiple fits; extend to vlines')
    df$type <- 'type1'
    df2     <- df
    df2$mpg <- df$mpg * 1.1 + rnorm(1, 10, 1)
    df2$type <- 'type2'
    df2 <- rbind(df, df2)
##     plotdf(df2, disp, mpg, byvar=type, multifit=TRUE, vlines=c(50, 500), main='test9: 2 fits by labels')


##     ## add another set of points and fit
##     plotdf(df2, disp, mpg, byvar=type, multifit=TRUE, vlines=c(50, 500), main='test10: add 3rd fit by hand')
##     xx <- seq(150, 400, 10)
##     rand <- rnorm(length(xx), 10, 1)
##     yy <- 5/100 * xx + rand
##     newcol <- 'darkviolet'
##     vlines <- c(120, 450)
## ##     addfit(xx, yy, col='darkviolet', vlines=c(120, 450))
##     addfit(xx, yy, col='darkviolet')

    
    ## separate fits for 2 variables with separate vlines and equations
    plotdf(df2, disp, mpg, byvar=type, nofit=TRUE, main='test10: add 3rd fit by hand')
    with(df2[df2$type == 'type1',], addfit(disp, mpg, vlines=c(80, 400), col='green'))
    with(df2[df2$type == 'type2',], addfit(disp, mpg, vlines=c(70, 550), col='red'))

    
    
}
