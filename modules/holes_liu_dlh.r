holes_liu_dlh <- function(space, data, sidesmall, reorder='yes') {
    ## Finds MHR (maximal empty hyper-rectangle) in a defined space
    ## Inputs should potentially be transposed and scaled if it is desired that:
    ## (1) an extent at a low value is equal to an extent at a high value
    ##     (e.g., pressure = 200, 400, 800, 1600 should be input as
    ##      log(pressure) if holes at 200 and 1600 should be considered
    ##      the same distance from 400 and 800, respectively)
    ## (2) an extent in 1 parameter is considered equal to the
    ##     same extent in another parameter (e.g., scale from 0 to 1)
    ##
    ## Input:  data      = dataframe containing datapoints
    ##                     format: 1st column = name of datapoint
    ##                             nth column = parameter values
    ##         space     = 2 row dataframe defining region of interest
    ##                     format: same as for data except name of "datapoint"
    ##                             in is SL in 1st row for the lower bound
    ##                             and   SU in 2nd row for the upper bound
    ##         sidesmall = MHR is too small to consider if
    ##                     volume < sidesmall^(# of parameters)
    ##         reorder   = reorder list of MHRs from largest to smallest
    ## Output: MHR         = sorted list of candidate MHRs from largest to smallest
    ##         MHRdropped  = sorted list of deleted because of being too small or unbounded
    ##         MHRreplaced = list of MHRs deleted because a datapoint exists inside it so not empty
    ##         volumesmall = sidesmall^(# of parameters)

    ## ## debug
    ## browser()

    ## ## debug
    ## space <- spaces
    ## data  <- dfs
    ## sidesmall <- 7
    ## reorder <- 'yes'

    ## rename data to df for use within function
    df <- data

    ## number of parameters
    nparam <- ncol(space) - 1
    end    <- nparam * 2
    
    ## Define first MHR h as the entire space and insert into the set of MHRs t
    ## MHR to be defined using lower and upper bounds
    ## hL <- data.frame(t(rep('SL', nparam)), stringsAsFactors=FALSE)
    ## names(hL) <- names(space[2:(nparam+1)])
    hL <- space[1, 2:(nparam+1)]
    names(hL) <- paste(names(hL), "_l", sep="")
    ## hU <- data.frame(t(rep('SU', nparam)), stringsAsFactors=FALSE)
    ## names(hU) <-  names(space[2:(nparam+1)])
    hU <- space[2, 2:(nparam+1)]
    names(hU) <- paste(names(hU), "_u", sep="")
    h <- cbind(hL,hU)
    ## add a list of additional bounding points on the MHR
    ## that may be needed to keep from inappropriately throwing out an MHR
    ## h$bounds <- list(bounds = as.character(h[1,]))
    ## h$bounds <- list(bounds = 'NA')
    t <- h

    ## define size for an MHR that is too small to consider
    volumespace <- size(hlocation(h, df, space))
    volumesmall <- sidesmall^nparam
    if (volumesmall > volumespace/4 ) {
        volumesmall <- volumespace/4
        cat('\n**** WARNING: sidesmall too big so reduced to 1/4 of space of interest ****\n')
    }
    
    ## initialize variables and dataframes
    hcount <- 1
    treplaced <- delete(t, t[1,], end)
    tdropped   <- treplaced

    ## For each point in database, search whether it is in t
    cat('number of datapoints', nrow(df),'\n')
    for (irow in 1:nrow(df)) {
        ##browser()
        ## if (irow >= 6) browser()
        pname <- as.character(df$name[irow])
        cat('starting datapoint', pname,'; number of MHRs =',nrow(t),'\n')
        ## print(pname)
        ## identify list of MHRs (each is an h) in the set of MHRs t that contain x
        ## for loop will modify t so need a temporary dataframe tprev with the previous values of t
        ## that the for loop can work through
        tprev <- t
        for ( hi in 1:nrow(tprev) ) {
            ## cat('irow=',irow,'; hi=',hi,'; jdim=',jdim,'\n')
            ## if (irow==12 & hi==130 & jdim==2) browser()
            h <- tprev[hi,]
            testloc <- ContainmentSearch(h, pname, df, space)$location
            if (testloc == 'on boundary') {
                ## insert x into set of bounding points in that surface
                ## t[hi,]$bounds[[1]] <- append(h$bounds[[1]], pname)
                ## ignore because point on boundary does not cause MHR to subdivide

            } else if (testloc == 'inside') {

                ## delete MHR h from collection of MHRs t
                t     <- delete(t, h, end)

                ## add MHR h to list of deleted MHRs for later tracing if needed
                treplaced <- insert(treplaced, h)

                ## break MHR h into 2 separate smaller MHRs along each parameter
                ## add the two new MHRs, ha and hb, to the set of MHRs t
                ## for each dimension
                for (jdim in 1:nparam) {
                    ## print(pname)
                    ## print(jdim)
                    ## start with the MHR of interest
                    ha <- h
                    hb <- h
                    ## print(h)

                    ## split h in each dimension, jdim
                    ## first assign location to the split location
                    datapoint <- df[df$name == pname,]
                    ha[jdim+nparam] <- datapoint[1 + jdim]
                    hb[jdim]        <- datapoint[1 + jdim]

                    ## identify name of the parent MHR
                    parent <- stri_split_fixed(str = row.names(h), pattern = "_", n = 2)[[1]][1]
                    
                    ## add new MHR ha if big enough and not completely contained in an existing MHR
                    haloc  <- hlocation(ha, df, space)
                    sizeha <- size(haloc)
                    ## if (pname == 'point3' & jdim == 2) browser()
                    ## name the MHR
                    hcount <- hcount+1
                    row.names(ha) <- paste(hcount, '_jdim', jdim, "L", "_parent_", parent, sep="")
                    ## add MHR to t if big enough and drop if too small
                    if ( sizeha >= volumesmall ) {
                        if ( nrow(t) == 0 ) {
                            ## t is empty so add ha
                            t <- rbind(t, ha)
                        } else if (  inexisting(t, haloc) == 'no' ) {
                            ## add ha because it is not fully contained inside another MHR
                            t <- rbind(t, ha)
                        }
                    } else {
                        ## drop the MHR
                        tdropped <- rbind(tdropped, ha)                   
                    }

                    ## add new MHR ha if big enough and not completely contained in an existing MHR
                    hbloc  <- hlocation(hb, df, space)
                    sizehb <- size(hbloc)
                    hcount <- hcount+1
                    row.names(hb) <- paste(hcount, '_jdim', jdim, "U", "_parent_", parent, sep="")
                    if ( sizehb >= volumesmall ) {
                        if ( nrow(t) == 0 ) {
                            t <- rbind(t, hb)
                        } else if (  inexisting(t, hbloc) == 'no' ) {
                            t <- rbind(t, hb)
                        }
                    } else {
                        tdropped <- rbind(tdropped, hb)                   
                    }
               
                }
            }        
            ## else testloc == 'outside' so look at the next MHR in t
        }
    }

    ## identify location of each side of MHRs
    ## strip off the boundary column for later lapply function
    tloc <- t[1:end]
    for (i in 1:nrow(tloc)) {
        ## tloc <- apply(t, 1, function(h) hlocation(h, df, space))
        h     <- tloc[i,]
        tloc[i,] <- as.numeric( hlocation(h, df, space) )
    }
    tloc
    ## tloc was still all character for some reason so need to convert to numeric
    tloc[] <- lapply(tloc, as.numeric)   

    ## add column of MHR sizes
    tloc$size <- size(tloc)
    t$size    <- tloc$size

    ## reorder by size
    if (reorder == 'yes') {
        tloc     <- tloc    [order(tloc$size,     decreasing=TRUE),]
        t        <- t       [order(t$size,        decreasing=TRUE),]
    }

    ## do the same for the dataframe of MHRs that were too small
    if (nrow(tdropped) > 0) {
        tlocdropped <- tdropped
        for (i in 1:nrow(tlocdropped)) {
            h     <- tlocdropped[i,]
            tlocdropped[i,] <- as.numeric( hlocation(h, df, space) )
        }
        tlocdropped[] <- lapply(tlocdropped, as.numeric)   # tlocdropped was all character
        tdropped$size <- size(tlocdropped)
        if (reorder == 'yes') tdropped <- tdropped[order(tdropped$size, decreasing=TRUE),]
    }

    return(list(MHRloc=tloc, MHR=t, MHRdropped=tdropped,
                MHRreplaced=treplaced, volumesmall=volumesmall))
            
}
