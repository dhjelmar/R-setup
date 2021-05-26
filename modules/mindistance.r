mindistance <- function(vec, df=dfxs) {
    ## Function to find distance from vector vec to nearest point in dataframe df
    # subtract vec from dataframe column
    dif   <- sweep(df, 2, vec, "-")
    # head(dif,2)
    # sumsquare each row of dataframe of diferences
    sq    <- function(x) { sum(x^2) }
    difsq <- apply(dif, 1, sq)
    # return the minimum distance
    mindist <- sqrt( min(difsq) )
    return(mindist)
}

negative_mindistance <- function(vec, df=dfxs){
    mindistance(vec,df) * -1
}

## Tests of function

## ## add center point to database for testing
## testdata <- dfxs
## vec <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
## testdata <- rbind(vec, testdata)
## head(testdata)
## 
## ## test to return 0.1
## vec1 <- vec + c(0, 0.1, 0, 0, 0, 0)
## ##vec1
## testanswer <- mindistance(vec1, testdata)
## cat("test to return 0.1 = ", testanswer,"\n")
## 
## ## test to return 0.14142
## vec1 <- vec + c(0, 0, 0.1, 0, 0.1, 0)
## ##vec1
## testanswer <- mindistance(vec1, testdata)
## cat("test to return 0.14142 = ", testanswer,"\n\n")
## 
## ## how about if supply vec from a row in a dataframe and look for min distance from any point in dfxs?
## cat("test to use function on vector from dfxs dataframe (same as above but no null row)\n")
## vecdf <- data.frame( t( vec ))
## cat("vecdf =\n")
## vecdf
## cat("convert to vector\n")
## vec <- as.numeric( vecdf[1,] )
## vec
## mindistance(vec, dfxs)
## 
## ## can also use apply to run mindistance function on each row of default dataframe, dfxs
## cat("\ntest to use function on all points in one dataframe to the dfxs dataframe\n")
## vecdf <- data.frame(rbind(vec, vec1))
## vecdf
## apply(vecdf, 1, mindistance)

