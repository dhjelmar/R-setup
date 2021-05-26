closestdata <- function(target, database, weight=NULL) {
    ## Function to sort database based on distance from target datapoint
    ## target   = dataframe containing only the parameters to be searched for
    ## database = dataframe containing database
    ##            can list more parameters than in target
    ## weight   = vector of weights for target parameters
    ##          = NULL results in value of 1 for all parameters
    ##          = 0 for a particular parameter, results in the value not mattering at all
    ##          = any large number for a particular parameter weight relative to scale of 
    ##            parameters can result in parameter having to be matched exactly

    ## create dataframe from database with only target parameters
    df <- subset(database, select=names(target))

    ## reorder subset dataframe to have the same order as target
    df <- df[names(target)]
    
    ## convert target to a vector and subtract from dataframe
    vec <- as.numeric(as.vector(target))
    dif <- sweep(df, 2, vec, "-")

    ## multiply each difference by the weights
    if (is.null(weight)) weight <- rep(1, ncol(target))
    difw  <- sweep(dif, 2, weight, "*")

    ## rename columns to indicate these are differneces
    colnames(difw) <- paste("dif", colnames(difw), sep = "_")
    
    ## sumsquare each row of dataframe of diferences and add to difference dataframe
    sq       <- function(x) { sum(x^2) }
    distance <- sqrt( apply(difw, 1, sq) )
    difw     <- cbind( distance, difw )
    
    ## combine the dataframe of weighted differences with the original dataframe then sort
    dfnew <- cbind(difw, database)
    dfnew <- dfnew[order(dfnew$distance),]
    
    ## return sorted dataframe
    return(dfnew)
}

## ## test
## df <- '
## name  param3 param1 param2
## point1     1      1      1
## point2     0      1      2
## point3    10     10     10
## point4     3      2      4
## '
## df <- readall(df)
## 
## closestdata(data.frame(t(c(param1=1.2, param2=2.1, param3=-0.1))), df)   # correctly returns points 2, 1, 4, 3
## 
## closestdata(data.frame(t(c(param1=2.5, param2=4.5, param3=10))), df)   # 4,3,1,2
## closestdata(data.frame(t(c(param1=2.5, param2=4.5, param3=10))), df,
##             weight=c(1, 1, 1000))   # weight correctly moves point 3 to the closest
## 
## closestdata(data.frame(t(c(param1=1, param2=2, param3=10))), df) # 4, 1, 2, 3
## closestdata(data.frame(t(c(param1=1, param2=2, param3=10))), df,
##             weight=c(1, 1, 0))   # weight correctly ignores param3 and correctly identifies 2 as closest
