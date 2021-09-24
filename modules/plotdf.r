plotdf2 <- function() {
    fits <- lme4::lmList(mpg ~ disp | cyl, data=mtcars)    
    coefs <- data.frame(coef(fits)); names(coefs) = c("Int", "slopes")
    plot(mpg ~ disp, data = mtcars)
    lapply(fits, abline)
}


plotdf <- function(df, xx, yy,
                   bycol      = NA,
                   byfit      = NA,
                   type       = NA,
                   ncol       = NA,
                   vlines     = NA,
                   xlimspec   = NA, ylimspec = NA,
                   main       = NULL,
                   printeq    = 'yes',
                   xlabel     = NA,
                   ylabel     = NA,
                   interval   = 'conf', alpha = 0.05, sided = 2,
                   bg         = "grey90",
                   suppress   = 'yes',
                   outputfile = NA) {
    ## create scatter plot with linear regression line and 95% confidence lines
    ## xx and yy are variables within dataframe df

    ## required inputs: df = dataframe
    ##                  xx = x-axis parameter
    ##                  yy = y-axis parameter

    ## options (defaults shown in function definition above):
    ##    bycol    = variable to use to color code datapoints (NA is default)
    ##    type     = NA is the default and behaves as if bycol = NA
    ##             = 'continuous'     = data colored using bycol as continuous parameter
    ##                                  single fit for all data combined
    ##             = 'discrete combined' = data colored using bycol as discrete label
    ##                                  single fit for all data combined
    ##             = 'discrete separate' = data colored using bycol as discrete label
    ##                                  separate fits for each discrete bycol

    ##    bycol    = variable to use to color code datapoints (NA is default)
    ##    byfit    = 'combined' specifies single fit to all data (default)
    ##             = 'bycol' specifies separate fit for each bycol

    ##    vlines   = vector of vertical lines to be added (e.g., c(0.5, 2.0))
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
    ##    bg       = background color for plot
    ##    suppress =
    ##    outputfile = name for PDF to be created
    ##               = NA does not create a PDF (only output is to screen)
    
    ## example usage: plotdf(mtcars, disp, mpg, bycol=cyl)

    ##-----------------------------------------------------------------------------
    ## PREPARE DATAFRAME AND EXTRACT LABELS FOR X, Y, AND BYCOL
    
    ## extract names of xx, yy, and bycol if not passed in with quotes
    xlabel  <- deparse(substitute(xx))
    ylabel  <- deparse(substitute(yy))
    bylabel <- deparse(substitute(bycol))

    ## attempt to replace xx, yy, and bycol with vector values
    xx    <- eval(substitute(xx),df)
    yy    <- eval(substitute(yy),df)    
    bycol <- eval(substitute(bycol),df) 

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
    if (typeof(bycol) == 'character' & length(bycol)==1)  {
        if (!is.na(bycol)) {
            ## bycol was specified
            bycol <- which(grepl(bycol, names(df)))  
            bycolv <- df[, bycol]
        } else {
            ## no bycol was specified
            ## create vector of values for bycol equal in length to xx and yy
            bycol <- rep(NA, length(xx))
        }
    }

    ## put xx, yy, and bycol into dataframe
    newdf <- data.frame(xx, yy, bycol=bycol)
    names(newdf) <- c('xx', 'yy', 'bycol')

    ##-----------------------------------------------------------------------------
    ## DETERMINE X AND Y RANGE FOR NEEDED FOR DATA AND REQUESTED LINES
    xmin <- min(newdf$xx, vlines, na.rm=TRUE)
    xmax <- max(newdf$xx, vlines, na.rm=TRUE)
    ymin <- min(newdf$yy, na.rm=TRUE)
    ymax <- max(newdf$yy, na.rm=TRUE)


    ##-----------------------------------------------------------------------------
    ## DETERMINE NUMBER OF DATAPOINT COLORS AND FITS NEEDED
    if (is.na(bycol[1])) {
        ## no bycol specified so single color and single fit group
        ncol      <- 1
        fitgroups <- 1
    } else {
        ## variable specified for datapoint colors
        
        if (is.na(ncol)) {
            ## number of colors not specified
            ncol <- length(unique(newdf$bycol))
        }

        if (is.na(byfit)) {
            ## default is a single fit group
            fitgroups <- 1
        } else {
            ## if anything other than NA is specified, then separate fit per color
            fitgroups <- ncol
        }
    }    

    ##-----------------------------------------------------------------------------
    ## DEFINE COLOR PALLETTE FOR PLOT
    cols <- c("black","darkturquoise","pink","red")
    pal  <- colorRampPalette(cols)                  # define color pallette
    if (ncol == 1) {
        newdf$color <- 'black'
    } else {
        if (is.numeric(newdf$bycol)) {
            ## continuous variable for color
            newdf$color <- pal(ncol)[as.numeric(cut(newdf$bycol,breaks=ncol))]
        } else {
            ## discrete variable (e.g., label) for color
            cols  <- data.frame(bycol=unique(c(as.character(newdf$bycol))))
            ncol  <- nrow(cols)
            color <- rainbow(ncol)
            color <- pal(ncol)
            ## to see colors
            ## scales::show_col(use)
            cols$color <- color
            newdf <- merge(newdf, cols, by="bycol")
        }
    }

    ##-----------------------------------------------------------------------------
    ## SPLIT DATAFRAME BY BYCOL IF MULTIPLE FITS ARE NEEDED
    if (fitgroups == 1) {
        ncol <- 1
        ## put newdf into list so can handle the same as in more complicated case for ncol > 1
        newdfl <- list(newdf)
    } else {
        ncol <- length(unique(bycol))
        newdfl    <- split(newdf, f = newdf$bycol)
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
    plot(xx, yy, type='n',
         xlim=xlimspec, ylim=ylimspec,
         xlab=xlabel, ylab=ylabel,
         main=main)
    
    ## changes background of plot to specified color
    par(bg=bg)  # changes background of plot to specified color

    ## add grid
    grid(col='gray70')

    ## if request upper/lower range bounds, add to plot
    if ( !is.na(vlines) ) {
        abline(v=vlines[1], col='black', lty='dashed')
        abline(v=vlines[2], col='black', lty='dashed')
    }        
    


    ##-----------------------------------------------------------------------------
    ## PERFORM REGRESSION AND PLOT FOR EACH BYCOL

    fit <- NA
    intercept <- NA
    slope <- NA
    rise <- NA
    equation <- NA
    legendn <- NA
    legendc <- NA
    
    for (ifit in 1:fitgroups) {
        fit             <- lm(yy~xx, data=newdfl[[ifit]])
        intercept[ifit] <- fit$coefficients[[1]]
        slope[ifit]     <- fit$coefficients[[2]]
        ## calculate rise of fit over range bounds if supplied,
        ## otherwise calculate rise over range of data
        if ( !is.na(vlines) ) {
            rise[ifit] <- (vlines[2] - vlines[1]) * slope
        } else {
            rise[ifit] <- (xmax - xmin) * slope
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
            legendn <- list( unique(newdfl[[ifit]]$bycol) )
            legendc <- list( unique(newdfl[[ifit]]$color) )
        } else if (fitgroups > 1 & ncol > 1) {
            ## bycol has a single color and subsequent ifit numbers will 
            legendn[ifit] <- newdfl[[ifit]]$bycol[1]
            legendc[ifit] <- newdfl[[ifit]]$color[1]
        }

        ## add fit to plot
        new.xx <- seq(min(newdfl[[ifit]]$xx, vlines, na.rm=TRUE),
                      max(newdfl[[ifit]]$xx, vlines, na.rm=TRUE), len=100)
        pred   <- predict(fit, new=data.frame(xx=new.xx), interval=interval, level=1-0.05/sided)
        lines(new.xx, pred[,"fit"], lwd=2)

        ## if request confidence bounds, add to plot
        if (interval != 'none') {
            lines(new.xx, pred[,"lwr"], lty=3)
            lines(new.xx, pred[,"upr"], lty=3)
        }


    }
    
    ##-------------------------------------------------------------------------------------     
    ## plot points and fit
    
    
    ## determine title for plot
    if (printeq == 'yes') {
        if (fitgroups == 1) {
            ## add fit equation under main title for the only fit
            mtext(equation[1], side=3, line=0, cex=1)
        } else if (fitgroups == 2) {
            ## add the 2 fit equations under main title
            subtitle <- list('line color: red = normal, blue = Weibull, black = Johnson',
                             'line type: solid = distribution or mean, dashed = 1-sided upper bound')
            mtext(subtitle, side=3, line=c(0.75, 0), cex=.75, col='black')
        } else{
            ## too many equations put them all in the title so skip
        }
    }
    
    ## add legend
    if (ncol > 1) {

## dlh stopped here
        
        if (bycol_type == 'discrete') {
            ## define color palette for discrete bycol
            legendcolor <- unique(newdfl$color)
            legendlabel <- unique(newdfl$bycol)
            ##legend("right",title=bylabel,col=pal(ncol), pch=19,legend=round(legendlabel))
            ##legend("bottomright",title=bylabel,col=newdfl$bycol, pch=19,legend=legendlabel)
            legend("bottomright", title=bylabel, col=legendcolor, pch=19, legend=legendlabel)

        } else {
            low  <- min(newdfl$bycol,na.rm=TRUE)
            high <- max(newdfl$bycol,na.rm=TRUE)
            legendlabel <- seq( low, high, (high-low)/(ncol-1) )
            ##legend("right",title=bylabel,col=pal(ncol), pch=19,legend=round(legendlabel))
            legend("right",title=bylabel,col=pal(ncol), pch=19,legend=signif(legendlabel,digits=4))
        }
    }
    
    ## for output file
    if (!missing(outputfile)) dev.off()
    
    ## print results of fit
    estbound(fit)
    
    if (suppress == 'no')
        return( list(intercept = intercept, slope = slope, rise = rise, pred = as_tibble(pred)) )
}


plotdf_test <- function() {
    ## example: plotdf(df,gnom,average_g_ts,TS_ID,bg="white")
    df <- mtcars
    df$size[df$cyl <= 4] <- 'small'
    df$size[df$cyl >  4] <- 'large'
    plotdf(df, disp, mpg,                            main='test 1')
    plotdf(df, disp, mpg, bycol=cyl,                 main='test 2')
    plotdf(df, disp, mpg, bycol=cyl,                 main='test 3')
    plotdf(df, disp, mpg, bycol=cyl,  byfit='bycol', main='test 4')
    plotdf(df, disp, mpg, bycol=size, byfit='bycol', main='test 5')
}

