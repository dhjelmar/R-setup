## # # Write insert, delete, and containment search functions

insert <- function(t, h) {
    ## inserts any MHR h into collection of MHRs t
    t <- rbind(t, h)
    return(t)
}
## ## test
## testh <- h*2
## test <- insert(t, testh)
## test <- insert(test, h)
## test

delete <- function(t, h, end) {
    ## deletes MHR h from collection of MHRs t
    matchingid <- which(apply(t[1:end], 1, function(x) all(x == h[1:end])))
    ## may not be present if previously deleted so check for matchingid > 0
    if (length(matchingid) > 0) t <- t[-c(matchingid),]
    return(t)
}
## delete(test, testh)

## return TRUE if x is a number and FALSE otherwise
is.number <- function(x) suppressWarnings(!is.na(as.numeric(as.character(x))))

## find parameter values for h based on location names
hlocation <- function(h, df, space) {
    nparam <- ncol(space) - 1
    end    <- nparam * 2
    hloc   <- h[1:end]
    for (i in 1:end) {
        boundname <- as.character(h[[i]])
        j <- i + 1
        if (i > nparam) j <- i - nparam + 1
        if (h[i] == 'SL' | h[i] == 'SU') {
            ## boundary is set to original space of interest
            ##print(i)
            ##print(as.character(h[[i]]))
            hloc[i] <- space[space$name == boundname, j]
        ## } else if (typeof(h[[i]]) == 'double') {
        } else if (is.number(h[[i]]) == 'TRUE') {
            ## implies boundary is set based on location and no boundary point yet
            ##print(i)
            hloc[i] <- h[i]
        } else {
            ## boundary is location of a datapoint on the boundary
            hloc[i] <- dfs[dfs$name == boundname, j]
        }
    }
    ## make sure value passed back is numeric
    hloc <- data.frame( t( sapply(hloc, function(x) as.numeric(as.character(x))) ) )
    return(hloc)
}
## hloc <- hlocation(h, dfs, space)

## # identifies where point x is relative to MHR h (i.e.,x is inside the bounds of the MHR)
ContainmentSearch <- function(h, pname, df, space) {
    ## # if pname is between every upper and lower bound for every parameter, then it is inside MHR
    ## # if pname is equal to upper or lower bound for any parameter, then it is on the boundary
    
    ## ## debug
    ## browser()
    
    ## what are the parameter values, x, for datapoint pname
    nparam <- ncol(space) - 1
    end    <- nparam*2
    datapoint <- df[df$name == pname,]
    x <- datapoint[2:(nparam+1)]
    if (nrow(x) == 0) return(list(location='ContainmentSearch: no matching point found',
                                  distance='0', mindist='0', maxdist='0'))
    if (nrow(x) > 1) cat('****** fatal error: duplciate datapoint names found ******\n')
    x <- as.numeric(x)
    
    ## what are the parameter values, hloc, for MHR h
    hloc <- hlocation(h[1:end], df, space)
    
    ## how far is pname from the low and high extremes of the MHR for each parameter
    ## low  <- hloc[1         :nparam]
    ## high <- hloc[(nparam+1):end   ]
    low  <- as.numeric( hloc[1         :nparam] )
    high <- as.numeric( hloc[(nparam+1):end   ] )
    distance <- (x - low)/(high - low)
    
    ## where is pname in relation to the MHR
    maxdist <- max(distance)
    mindist <- min(distance)
    if ( (mindist < 0) | (maxdist > 1) ) {
        location <- 'outside'
    } else if ( (mindist == 0) | (maxdist == 1) ) {
        location <- 'on boundary'
    } else {
        location <- 'inside'
    }
    
    return(list(location=location, distance=distance, mindist=mindist, maxdist=maxdist))
    ## return(location)
}
## h
## pname <- 'point1'
## contains <- ContainmentSearch(h, pname, dfs, spaces)
## contains

## ## following finds all matching points but that was not what was needed
## # ContainmentSearch <- function(t,x) which(apply(t, 1, function(xx) all(xx == x)))
## # ContainmentSearch(t,dfx)

size <- function(hloc) {
    ## use volume of hyper rectangle
    end    <- ncol(hloc)
    nparam <- end/2
    low  <- hloc[1         :nparam]
    high <- hloc[(nparam+1):end   ]
##    browser()
    volume <- apply(high - low, 1, prod)
    ## low  <- as.numeric( hloc[1         :nparam] )
    ## high <- as.numeric( hloc[(nparam+1):end   ] )
    ## volume <- prod(high-low)
    return(volume)
}


bigenough <- function(hloc, small=sizesmall) {
    ## return TRUE if MHR h is sufficiently large
    sizeh <- size(hloc)
    logical <- sizeh >= small
    return(logical)
}
## #sizesmall
## #size(h)
## #bigenough(h,(1e+6-1))
## #bigenough(h,1e+6)
## #bigenough(h,(1e+6+1))
## #bigenough(h,sizesmall)

bounded <- function(h, end) {
    ## return true if the MHR, h, is bounded on all sides by a data point or original MHR side
    ## irow is the row number for the current data point of interest in df
    unbounded <- sapply(h[1:end], is.numeric)
    if (ncol(h[unbounded]) == 0) {
        bounded <- TRUE
    } else {
        bounded <- FALSE
    }
    return(bounded)
}

pbound <- function(h, jdim, side, df, space) {
    ## for parameter jdim:
    ##     return the name of the boundary if bounded by SL, SU, or a datapoint
    ##     otherwise return the location of the boundary

    ## print(h)
    ## print(jdim)
    nparam <- ncol(space) - 1
    if (side == "L" ) {
        jside <- jdim
    } else {
        jside <- jdim + nparam
    }
    hchar <- as.character(h[[jside]])  ## pull h[jside] out of dataframe
    ## browser()

    ## check boundary points and replace if needed
    if (h[jside] == 'SL' | h[jside] == 'SU') {
        ## still bounded by original space
        boundedby <- hchar
        
    } else if (is.numeric(h[jside])) {
        ## parent was not bounded either just pass the same
        boundedby <- as.numeric(h[[jside]])
        
    } else {
        ## test that the parent point on the boundary is really on the new MHR boundary
        testloc <- ContainmentSearch(h, hchar, df, space)$location
        if (testloc == 'on boundary') {
            boundedby <- hchar
            
        } else {
            ## parent point not on boundary so reset boundary to location value
            hloc <- hlocation(h, df, space)
            boundedby <- hloc[jside]

            ## now look at additional boundary points
            ## and use if on the boundary of interest (jside)
            bounds  <- h$bounds[[1]]
            lbounds <- length(bounds)
            ## only look if bounds contains more than the 1 NA point
            if (lbounds > 1) {
                bounds <- bounds[2:lbounds]  # stripped out the NA
                for (boundpoint in bounds) {
                    cs <- ContainmentSearch(h, boundpoint, df, space)
                    if (cs$location == 'on boundary' &
                        df[df$name == boundpoint, 1+jdim] == hloc[[jside]] ) {
                        ## then point is on the same boundary
                        boundedby <- boundpoint
                        break
                    }
                }
            }
        }
    }
    return(boundedby)
}


inexisting <- function(t, h) {
    ## return no if h is not conatined within an existing MHR in t
    ## return yes if h is contained within an existing MHR in t
    ## t and h need to contain parameter values
    nparam <- ncol(h) / 2
    inexisting <- 'no'
    for (trow in 1:nrow(t)) {
        ## nwithin will be used to track the number of parameters for which h is within trow
        nwithin <- 0
        for (j in 1:nparam) {
            ## for parameter j, test to see if MHR h is above lower bound of trow
            ## and below upper bound of trow
            if (h[j] >= t[trow, j] & h[j+nparam] <= t[trow, j+nparam]) nwithin <- nwithin + 1
        }
        if (nwithin == nparam) {
            inexisting <- 'yes'
            break
        }
    }
    return(inexisting)
}
            


rglcube <- function(data, color='blue') {
    ## adds a rectangular prism to a scatter plot made using rgl::plot3d
    ## data = dateframe where 1st 6 columns are: xL, yL, zL, xU, yU, zU
    ##        additional columns, if present, are ignored
    ## polygon3d is part of rgl packages
    ## 
    liu <- data
    polygon3d(x=c(liu[[1]], liu[[4]], liu[[4]], liu[[1]]), # x-y plane at lower z
              y=c(liu[[2]], liu[[2]], liu[[5]], liu[[5]]),
              z=c(liu[[3]], liu[[3]], liu[[3]], liu[[3]]),
              coords=c(1,2), col=color, alpha=0.2)
    polygon3d(x=c(liu[[1]], liu[[4]], liu[[4]], liu[[1]]), # x-y plane at upper z
              y=c(liu[[2]], liu[[2]], liu[[5]], liu[[5]]),
              z=c(liu[[6]], liu[[6]], liu[[6]], liu[[6]]),
              coords=c(1,2), col=color, alpha=0.2)
    polygon3d(x=c(liu[[1]], liu[[1]], liu[[1]], liu[[1]]),  # y-z plane at lower x
              y=c(liu[[2]], liu[[2]], liu[[5]], liu[[5]]),
              z=c(liu[[6]], liu[[3]], liu[[3]], liu[[6]]),
              coords=c(2,3), col=color, alpha=0.2)
    polygon3d(x=c(liu[[4]], liu[[4]], liu[[4]], liu[[4]]),  # y-z plane at upper x
              y=c(liu[[2]], liu[[2]], liu[[5]], liu[[5]]),
              z=c(liu[[6]], liu[[3]], liu[[3]], liu[[6]]),
              coords=c(2,3), col=color, alpha=0.2)
    polygon3d(x=c(liu[[1]], liu[[4]], liu[[4]], liu[[1]]),  # z-x plane at lower y
              y=c(liu[[2]], liu[[2]], liu[[2]], liu[[2]]),
              z=c(liu[[3]], liu[[3]], liu[[6]], liu[[6]]),
              coords=c(3,1), col=color, alpha=0.2)
    polygon3d(x=c(liu[[1]], liu[[4]], liu[[4]], liu[[1]]),  # z-x plane at upper y
              y=c(liu[[5]], liu[[5]], liu[[5]], liu[[5]]),
              z=c(liu[[3]], liu[[3]], liu[[6]], liu[[6]]),
              coords=c(3,1), col=color, alpha=0.2)
}
