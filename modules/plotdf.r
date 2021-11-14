plotdf <- function(df, xx, yy, multifit=TRUE, interval='conf', ...) {

    ## extract desired parts of df
    ## not quite sure why all_of() is needed; it was suggested in a warning when run without it
    dfnew <- select(df, c(all_of(xx), all_of(yy)))
    
    
    ## melt creates dataframe with only 3 paramters: id.vars,         variable.name, and "value"
    ##                                              "duration_years", "series",      and "value"
    dflong <- reshape2::melt(dfnew, id.vars=xx, variable.name='series')
    ## with(dflong, plotfit(xx, value, series, multifit=multifit, interval=interval, ...))
    plotfit(dflong[[1]], dflong[[3]], dflong[[2]], 
            xlabel=xx, ylabel='value', bylabel='series',
            multifit=multifit, interval=interval, ...)
    
}
## plotspace(2,3)
## with(mtcars, plotfit(mpg, drat, ylimspec=c(0, 4)))
## with(mtcars, plotfit(mpg, cyl))
## plotdf(mtcars, 'mpg', c('cyl', 'drat'))
## plotdf(mtcars, 'mpg', c('cyl', 'drat'), multifit=FALSE)
## plotdf(mtcars, 'mpg', c('cyl', 'drat'), multifit=FALSE, interval='line')
## plotdf(mtcars, 'mpg', c('cyl', 'drat'), multifit=TRUE,  interval='line')

