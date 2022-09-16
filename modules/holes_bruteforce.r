holes_bruteforce <- function(candidatepoints, database) {
    ## calculate the minimum distance from the test point to the nearest point in the database.

    ## candidatepoints = dataframe of parameters
    ## database        = dataframe of database parameters (i.e., datapointnames removed)
    
    ## calculate mininimum distance from each candidatepoints row to a database point
    ## in following apply statement, each row of candidatepoints is passed to the
    ## mindistance function as vector vec. Mindistance also requires a dataframe argument.
    candidatepoints$mindist <- apply(candidatepoints, 1, function(vec) mindistance(vec, database))

    ## sort results from max to min distance
    holes <- candidatepoints[order(candidatepoints$mindist, decreasing=TRUE),]

    return(holes)
}
