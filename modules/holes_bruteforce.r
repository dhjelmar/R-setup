holes_bruteforce <- function(candidatepoints, database, select='all', max.to.min=TRUE) {
    ## calculate the minimum distance from the test point to the nearest point in the database.

    ## candidatepoints = dataframe of parameters to use in evaluation
    ## database        = dataframe of database parameters (can include extra parameters)
    ## select          = vector of parameters which can be a subset of what is in database
    ##                 = 'all' (default) uses all paramters in candidatepoints
    ## max.to.min      = TRUE (default) sorts results from max to min distance

    if (select[1] == 'all') {
        ## check to make sure database parameters are consistent with candidatepoints
        if (identical(candidatepoints, database)) {
            ## same names exactly in both dataframes so no adjustment is needed
        } else {
            ## extract and/or reorder parameters needed from database
            database <- subset(database, select=names(candidatepoints))
        }

    } else {
        ## select requesteed parameters from candidatepoints and database
        candidatepoints <- subset(candidatepoints, select=select)
        database        <- subset(database       , select=select)
    }
            
    ## calculate mininimum distance from each candidatepoints row to a database point
    ## in following apply statement, each row of candidatepoints is passed to the
    ## mindistance function as vector vec. Mindistance also requires a dataframe argument.
    ## candidatepoints$mindist <- apply(candidatepoints, 1, function(vec) mindistance(vec, database))
    ## potentially faster to execute as a matrix rather than a dataframe
    candidatepoints$mindist <- apply(as.matrix(candidatepoints), 1, function(vec) mindistance(vec, as.matrix(database)))

    if (isTRUE(max.to.min)) {
        ## sort results from max to min distance
        holes <- candidatepoints[order(candidatepoints$mindist, decreasing=TRUE),]
    } else {
        holes <- candidatepoints
    }
    
    return(holes)
}

## data <- mtcars
## testpoints <- data.frame(mpg=c(18, 20, 25),
##                          cyl=c(6,4,2),
##                          disp=c(200, 200, 200),
##                          hp=c(100, 120, 110), 
##                          qsec=c(15, 15, 15))
## holes_bruteforce(candidatepoints=testpoints, database=data, select='all', max.to.min=TRUE)
## holes_bruteforce(candidatepoints=testpoints, database=data, select=c('qsec', 'cyl'), max.to.min=TRUE)
