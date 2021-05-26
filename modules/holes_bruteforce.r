    holes_bruteforce <- function(candidatepoints, database) {
        ## calculate the minimum distance from the test point to the nearest point in the database.
        ## candidatepoints = dataframe of parameters
        ## database        = dataframe of database parameters (i.e., datapointnames removed)
        dfxs <- database   # dfxs is needed for mindistance function
        candidatepoints  <- candidatepoints
        nparam <- ncol(candidatepoints)
        ## calculate mininimum distance to new dfxs data set
        candidatepoints$mindist <- apply(candidatepoints, 1, mindistance)
        ## sort results from max to min distance
        holes <- candidatepoints[order(candidatepoints$mindist, decreasing=TRUE),]
        return(holes)
    }
