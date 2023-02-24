plotspace <- function(rows, columns, nplots=NA) {
    ## input: rows    = number of plots per row
    ##        columns = number of plots per column
    ##        nplots  = NA (default) ignors this input
    ##                = number overwrites any specificed number of rows and columns
    ##                  if provided to instead best fit all plots

    if (!is.na(nplots)) {
        ## figure out the best arrangement for plots
        rows    <- ceiling(nplots^0.5)
        columns <- ceiling(nplots / rows)
    }

    ## set number of rows and columns
    par(mfrow=c(rows,columns))
}

## plotspace(2,2)
## plot(mtcars$cyl , mtcars$mpg)
## plot(mtcars$disp, mtcars$mpg)
## plot(mtcars$wt  , mtcars$mpg)
## 
## plotspace(nplots=11)
## plot(mtcars$cyl , mtcars$mpg)
## plot(mtcars$disp, mtcars$mpg)
## plot(mtcars$wt  , mtcars$mpg)
## plot(mtcars$gear, mtcars$mpg)
## plot(mtcars$carb, mtcars$mpg)
## plot(mtcars$hp  , mtcars$mpg)
## 
## plotspace(2,2,nplots=5)
## plot(mtcars$cyl , mtcars$mpg)
## plot(mtcars$disp, mtcars$mpg)
## plot(mtcars$wt  , mtcars$mpg)
## plot(mtcars$gear, mtcars$mpg)
## plot(mtcars$carb, mtcars$mpg)
## plot(mtcars$hp  , mtcars$mpg)

