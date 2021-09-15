plotdf <- function(df, xx, yy,
                   byvar      = NA,
                   byvar_type = 'discrete',
                   fit_type   = 'combined',
                   ncol       = NULL,
                   vlines     = NULL,
                   xlimspec   = NULL, ylimspec = NULL,
                   main       = 'equation', xlabel = NULL, ylabel = NULL,
                   interval   = 'conf', alpha = 0.05, sided = 2,
                   bg         = "grey90",
                   suppress   = 'yes',
                   outputfile = NULL) {
    ## create scatter plot with linear regression line and 95% confidence lines
    ## xx and yy are variables within dataframe df

    ## required inputs: df = dataframe
    ##                  xx = x-axis parameter
    ##                  yy = y-axis parameter

    ## options (defaults shown in function definition above):
    ##    byvar = variable to use to color code datapoints
    ##    byvar_type = 'discrete' or 'continuous'; only used if byvar != NA
    ##    fit_type   = 'combined' = combines all data in single fit
    ##               = 'separate' = performs separate fit for each byvar (only for byvar_type = discrete)
    ##    vlines   = vector of vertical lines to be added (e.g., c(0.5, 2.0))
    ##    xlimspec = vector defining xlimits of plot (e.g., c(-2, 11.5))
    ##    ylimspec = vector defining ylimits of plot
    ##    main     = use the equation yy = m * xx + b for the title
    ##             = option for user specified title
    ##    xlabel   = NULL = x-variable name used for x-axis
    ##             = user specified for label for x-axis
    ##    ylabel   = NULL = y-variable name used for y-axis
    ##             = user specified for label for y-axis
    ##    inverval = 'conf' adds confidence limits to the plot
    ##    alpha    = 1 - confidence
    ##    sided    = 2 = 2-sided confidence limits
    ##             = 1 = 1-sided confidence limits (but both sides plotted)
    ##    bg       = background color for plot
    ##    suppress =
    ##    outputfile = name for PDF to be created
    ##               = NULL does not create a PDF (only output is to screen)
    
    ## example usage: plotdf(mtcars, disp, mpg, byvar=cyl)
    
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
            bycol <- which(grepl(byvar, names(df)))  
            byvar <- df[, bycol]
        } else {
            ## no byvar was specified
            ## create vector of values for byvar equal in length to xx and yy
            byvar <- rep(NA, length(xx))
        }
    }
    
    ## put xx, yy, and byvar into dataframe
    newdf <- data.frame(xx,yy,byvar)
    names(newdf) <- c('xx', 'yy', 'byvar')

    ## perform regression
    if (byvar[1] == 'NA' | fit_type == 'combined') {
        fits <- 1
    } else {
        fits <- length(unique(byvar))
    }
    
    intercept <- NA
    slope <- NA
    rise <- NA
    title <- NA
    for (ifit in 1:fits) {
        fit             <- lm(yy~xx, data=newdf)
        intercept[ifit] <- fit$coefficients[[1]]
        slope[ifit]     <- fit$coefficients[[2]]
        ## calculate rise of fit over range bounds if supplied, otherwise calculate rise over range of data
        if ( !is.null(vlines) ) {
            rise[ifit] <- (vlines[2] - vlines[1]) * slope
        } else {
            rise[ifit] <- (xmax - xmin) * slope
        }        
        title[ifit] = paste0("y = ", signif(slope[ifit],4), " * x + ", signif(intercept[ifit],4),
                             ", ", expression(Delta), " = ", signif(rise[ifit],4))
    }
    
    ##-------------------------------------------------------------------------------------     
    ## plot points and fit
    
    ## determine x and y range for plot
    xmin <- min(newdf$xx, xlimspec, vlines, na.rm=TRUE)
    xmax <- max(newdf$xx, xlimspec, vlines, na.rm=TRUE)
    if (is.na(byvar[1])) {
      ## do not need extra room for legend since no colors on plot needing legend
      xmax_plot <- xmax
    } else {
      ## make room for legend
      xmax_plot <- xmax + (xmax-xmin)*0.2  
    }
    ymin <- min(newdf$yy, na.rm=TRUE)
    ymax <- max(newdf$yy, na.rm=TRUE)
    if (missing(xlimspec)) { xlimspec <- c(xmin,xmax_plot) }
    if (missing(ylimspec)) { ylimspec <- c(ymin,ymax) }
    
    ## determine title for plot
    if (fits == 1 & main == 'equation') {
        main = title[1]
    } else if (fits > 1) {
        ## potentially too many equations put them all in the title
        main = NULL
    } else {
        ## keep user supplied title
        main = main
    }
    ## setup for output jpeg file
    if (!missing(outputfile)) jpeg(filename=outputfile)

    ## changes background of plot to specified color
    par(bg=bg)  # changes background of plot to specified color

    ## add color definition for each point using byvar
    ##cols <- c("blue","red")                         # define colors
    ##cols <- brewer.pal(11,"Spectral")             # Spectral color pallette has ncol=11 colors
    ##cols <- c("black","blue","green","grey","orange","red")
    ##cols <- c("black","blue","darkorchid2","darkturquoise","darkgreen","green","deeppink","magenta","red")
    ##cols <- c("black","darkturquoise","pink","red")
    ##pal  <- colorRampPalette(cols)                  # define color pallette
    cols <- data.frame(byvar=unique(c(as.character(newdf$byvar))))
    if (nrow(cols) > 1) {

        if (byvar_type == 'discrete') {
            ## define color palette for discrete byvar
            ncol <- nrow(cols)
            pal <- rainbow(ncol)
            ## library(scales) # needed for show_col
            ## show_col(pal)
            cols$color <- pal
            newdf <- merge(newdf, cols, by="byvar")
            ## plot data points
            plot(newdf$xx, newdf$yy, xlim=xlimspec, ylim=ylimspec,
                 xlab=xlabel, ylab=ylabel, main=main,
                 col=newdf$color)

        } else {
            ## define color palette for continuous byvar
            if (is.null(ncol)) ncol <- nrow(cols)
            cols <- c("black","darkturquoise","pink","red")
            pal  <- colorRampPalette(cols)                  # define color pallette
            newdf$color <- pal(ncol)[as.numeric(cut(newdf$byvar,breaks=ncol))]
            ## plot data points
            plot(newdf$xx,newdf$yy,xlim=xlimspec,ylim=ylimspec,
                 xlab=xlabel,ylab=ylabel, main=main,
                 col=pal(ncol)[as.numeric(cut(newdf$byvar,breaks=ncol))])
        }
        
    } else {
        ## use black if no byvar
        cols$color <- 'black'
        newdf <- merge(newdf, cols, by="byvar")
        ## plot data points
        plot(newdf$xx, newdf$yy, xlim=xlimspec, ylim=ylimspec,
             xlab=xlabel, ylab=ylabel, main=main)
}
    
    ## add grid
    grid(col='gray70')

    ## add fit to plot
    new.xx <- seq(min(newdf$xx, vlines, na.rm=TRUE),max(newdf$xx, vlines, na.rm=TRUE), len=100)
    fit    <- lm(yy~xx,data=newdf)
    pred   <- predict(fit, new=data.frame(xx=new.xx), interval=interval, level=1-0.05/sided)
    lines(new.xx,pred[,"fit"],lwd=2)

    ## if request confidence bounds, add to plot
    if (interval != 'none') {
        lines(new.xx,pred[,"lwr"],lty=3)
        lines(new.xx,pred[,"upr"],lty=3)
    }

    ## if request upper/lower range bounds, add to plot
    if ( !is.null(vlines) ) {
        abline(v=vlines[1], col='black', lty='dashed')
        abline(v=vlines[2], col='black', lty='dashed')
    }        
    
    ## add legend
    if (ncol > 1) {

        if (byvar_type == 'discrete') {
            ## define color palette for discrete byvar
            legendcolor <- unique(newdf$color)
            legendlabel <- unique(newdf$byvar)
            ##legend("right",title=bylabel,col=pal(ncol), pch=19,legend=round(legendlabel))
            ##legend("bottomright",title=bylabel,col=newdf$byvar, pch=19,legend=legendlabel)
            legend("bottomright", title=bylabel, col=legendcolor, pch=19, legend=legendlabel)

        } else {
            low  <- min(newdf$byvar,na.rm=TRUE)
            high <- max(newdf$byvar,na.rm=TRUE)
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
## example: plotfitcold(df,gnom,average_g_ts,TS_ID,bg="white")
##          plotfitcold(mtcars, disp, mpg, byvar=cyl)

