holes_liu <- function(space, data, sidesmall, requirebounded='yes', reorder='yes') {
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
    ##         requirebounded = yes removes unbounded MHRs during search as recommended in Liu paper
    ##                          this is the faster option but I am not convinced it will not drop a valid MHR
    ##                        = no removes unbounded MHRs after completion of search
    ## Output: MHR         = sorted list of candidate MHRs from largest to smallest
    ##         MHRdropped  = sorted list of deleted because of being too small or unbounded
    ##         MHRreplaced = list of MHRs deleted because a datapoint exists inside it so not empty
    ##         volumesmall = sidesmall^(# of parameters)

    ## debug
    ## browser()
    space <- spaces
    data  <- dfs
    sidesmall <- 0.2
    requirebounded <- 'yes'
    reorder <- 'yes'

    ## rename data to df for use within function
    df <- data

    ## number of parameters
    nparam <- ncol(space) - 1
    end    <- nparam * 2
    
    ## Define first MHR h as the entire space and insert into the set of MHRs t
    ## MHR to be defined using lower and upper bounds
    hL <- data.frame(t(rep('SL', nparam)), stringsAsFactors=FALSE)
    names(hL) <- names(space[2:(nparam+1)])
    names(hL) <- paste(names(hL), "_l", sep="")
    hU <- data.frame(t(rep('SU', nparam)), stringsAsFactors=FALSE)
    names(hU) <-  names(space[2:(nparam+1)])
    names(hU) <- paste(names(hU), "_u", sep="")
    h <- cbind(hL,hU)
    ## add a list of additional bounding points on the MHR
    ## that may be needed to keep from inappropriately throwing out an MHR
    h$bounds <- list(bounds = 'NA')
    t <- h

    ## define size for an MHR that is too small to consider
    ## sidesmall <- 0.1
    volumesmall <- sidesmall^nparam

    ## initialize variables and dataframes
    hcount <- 1
    treplaced <- delete(t, t[1,], end)
    tdropped   <- treplaced

    ## For each point in database, search whether it is in t
    for (irow in 1:nrow(df)) {
        ## browser()
        ## if (irow >= 6) browser()
        pname <- as.character(df$name[irow])
        ## print(pname)
        ## identify list of MHRs (each is an h) in the set of MHRs t that contain x
        ## for loop will modify t so need a temporary dataframe tprev with the previous values of t
        ## that the for loop can work through
        tprev <- t
        for ( hi in 1:nrow(tprev) ) {
            ## browser()


            
            if (irow >= 12 & hi > 50) browser()


            
            h <- tprev[hi,]
            testloc <- ContainmentSearch(h, pname, df, space)$location
            if (testloc == 'on boundary') {
                ## insert x into set of bounding points in that surface
                ##t <- delete(t, h, end)
                ##h$bounds[[1]] <- append(h$bounds[[1]], pname)
                ##t <- insert(t, h)
                t[hi,]$bounds[[1]] <- append(h$bounds[[1]], pname)
            } else if (testloc == 'inside') {
                ## line 9
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
                    ## if pname is on the boundary, swap location for pname
                    datapoint <- df[df$name == pname,]
                    ha[jdim+nparam] <- datapoint[1 + jdim]
                    hb[jdim]        <- datapoint[1 + jdim]
                    ## print(ha)
                    ## print(hb)
                    if (ContainmentSearch(ha, pname, df, space)$location == 'on boundary')
                        ha[jdim+nparam] <- pname
                    if (ContainmentSearch(hb, pname, df, space)$location == 'on boundary')
                        hb[jdim]        <- pname
                    ## print(ha)
                    ## print(hb)

                    ## delete additional boundary points and only add back if they are still on the new MHR
                    bounds    <- h$bounds[[1]]
                    numbounds <- length(bounds)
                    ha$bounds  <- list(bounds = 'NA')
                    hb$bounds  <- list(bounds = 'NA')
                    if (numbounds > 1) {
                        for (k in 2:numbounds) {
                            if (ContainmentSearch(ha, bounds[k], df, space)$location == 'on boundary')
                                ha$bounds[[1]] <- append(ha$bounds[[1]], bounds[k])
                            if (ContainmentSearch(hb, bounds[k], df, space)$location == 'on boundary')
                                hb$bounds[[1]] <- append(hb$bounds[[1]], bounds[k])
                        }
                    }

                    for (jjdim in 1:nparam) {
                        ## replace boundary points if subdivided MHR
                        ## no longer has boundary points
                        if (jjdim != jdim) {
                            ## already addressed jdim so can skip
                            ## if (irow==3 & hi==4 & jdim==1) browser()
                            ha[jjdim]        <- pbound(ha, jjdim, "L", df, space)
                            ha[jjdim+nparam] <- pbound(ha, jjdim, "U", df, space)
                            hb[jjdim]        <- pbound(hb, jjdim, "L", df, space)
                            hb[jjdim+nparam] <- pbound(hb, jjdim, "U", df, space)
                        }
                    }

                    ## identify name of the parent MHR
                    parent <- stri_split_fixed(str = row.names(h), pattern = "_", n = 2)[[1]][1]
                    
                    ## add new MHR ha if big enough and bounded
                    haloc  <- hlocation(ha, df, space)
                    sizeha <- size(haloc)
                    ## if (pname == 'point3' & jdim == 2) browser()
                    hcount <- hcount+1
                    row.names(ha) <- paste(hcount, '_jdim', jdim, "L", "_parent_", parent, sep="")
                    if ( sizeha >= volumesmall & (requirebounded == 'no' | bounded(ha,end) == TRUE) ) {
                        ## row.names(ha) <- paste(pname, '__parent_',row.names(h),
                        ##                        '__jdim_',jdim, '__L', sep="")
                        t <- rbind(t, ha)
                    } else {
                        ## too small or unbounded by datapoint
                        ## if the latter, the MHR should be picked up elsewhere
                        tdropped <- rbind(tdropped, ha)                   
                    }

                    ## add new MHR ha if big enough and bounded
                    hbloc  <- hlocation(hb, df, space)
                    sizehb <- size(haloc)
                    hcount <- hcount+1
                    row.names(hb) <- paste(hcount, '_jdim', jdim, "U", "_parent_", parent, sep="")
                    if ( sizehb >= volumesmall & (requirebounded == 'no' | bounded(hb,end) == TRUE) ) {
                        ## big enough and either do not care if bounded or it is bounded
                        t <- rbind(t, hb)
                    } else {
                        ## if (row.names(hb) == 'num_jdim1U_parent_171') browser()
                        ## if (jdim == 'L' & parent == '17') browser()
                        tdropped <- rbind(tdropped, hb)                                       
                    }

                }

                ## if (irow==6 & hi==8 & hcount>36) browser()

             
                ## delete MHR h from collection of MHRs t
                t     <- delete(t, h, end)
                ## add MHR h to list of deleted MHRs for later tracing if needed
                treplaced <- insert(treplaced, h)
            }        
            ## else testloc == 'outside' so look at the next MHR in t 
        }
    }

    if (requirebounded == 'no') {
        ## chose option to not remove unbounded MHRs above so doing it now
        for ( hi in 1:nrow(t) ) {
            ## browser()
            h <- t[hi,]
            if (bounded(h,end) == FALSE) t <- delete(t, h, end)
            browser()
            tdropped <- insert(tdropped, h)
        }
    }

    ## identify location of each side of MHRs
    tloc <- t
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

    ## do the same for the dataframe of too small or unbounded MHRs, if any
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
