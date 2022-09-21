df.duplicates <- function(df, tol=0.01, vector=NA, param=NA, tol.type='fraction') {
    ## input:   df       = numeric, 2D dataframe
    ##          tol      = tolerance
    ##          tol.type = 'fraction' indicates tolerance is the 
    ##                     specified fraction of value being compared 
    ##                     (i.e., tolerance = tol * max(vector + df[i,j], 1E-6)
    ##                   = 'absolute' indicates tolerance is the
    ##                      actual value of the tolerance on df and vec
    ##                      (i.e., tolerance = 2*tol)

    ## options: vector   = NA (default) looks for duplicates of all rows in df
    ##                   = vector looks for duplicates only of supplied vector
    ##          param    = vector of parameters to use from within df and vector

    ## renumber df
    rownames(df) <- 1:nrow(df)
    
    ##-----------------------------------------------------------------------------
    ## determine which columns of df and, if specified, vector to use in search
    if (!is.na(param[1])) {
        ## param is used to determine which columns of df and vector to use
        df.cols  <- names(df) %ino% param
        if (!is.na(vector[1])) {
            vec.cols <- names(vector) %ino% param
        } else {
            vec.cols <- df.cols
        }
            
    } else if (!is.na(vector[1])) {
        ## vector is used to determine which columns of df to use (all vector columns used)
        df.cols  <- names(df) %ino% names(vector)
        vec.cols <- 1:length(vector)

    } else {
        ## only df is specified, so use all parameters
        df.cols <- 1:ncol(df)
        vec.cols <- df.cols
    }
    nparam <- length(df.cols)
        
    ##-----------------------------------------------------------------------------
    ## search for near duplicates
    df.out <- df
    df$duplicate <- NA
    dups.to.remove.all <- NA
    for (row in 1:nrow(df)) {
        
        ## define vector, vec, to be compared to each row of df
        if (is.na(vector[1])) {
            vec <- df[row,]
        } else {
            vec <- vector
        }
        
        out <- data.frame(matrix(nrow=nrow(df.out), ncol=nparam))
        
        for (i in 1:nparam) {
            ## use lapply to look at every row in shrinking df.out for same value (within tol) in each column
            vec.col.i     <- vec[[vec.cols[i]]]
            df.out.col.i  <- df.out[,df.cols[i]]
            if (tol.type == 'fraction') {
                ## use fractional tolerance
                tol.col.i <- tol*vec.col.i
                out[,i] <- unlist(lapply(df.out.col.i, function(x) dplyr::near(x, vec.col.i,
                                                                               tol=tol*max(x+vec.col.i, 1E-6))))
            } else {
                ## use absolute tolerance
                tol.col.i <- tol
                out[,i] <- dplyr::near(df.out.col.i, vec.col.i, tol=2*tol.col.i)
            }
        }  
        
        ## duplcates are any rows in df that are the same in every column
        duplicates <- list(which(apply(out, 1, function(x) all(x, na.rm=TRUE) )))    # all tests entire vector for true
        dup.vec <- as.vector(duplicates[[1]])
        ## above is the row(s) in df.out; what row(s) is dup.vec in df?
        dup.vec.df <- as.numeric(rownames(df.out[dup.vec,]))

        if (!is.na(vector[1])) {
            ## only need this one time through
            df.out <- df.out[dup.vec,]
            if (!is.na(param[1])) {
                if (is.data.frame(vector)) {
                    ## if vector input was actually entered as a dataframe row
                    ## then can directly select paramters used in the searh using subset()
                    vector <- subset(vector, select=param)
                } else {
                    ## if vector was actually entered as a vactor
                    ## then need to transpose vector and turn it into a dataframe first
                    vector <- subset(as.data.frame(t(vector)), select=param)
                }
            }
            df.out <- fastmerge(df.out, vector)
            ## colnames(df.out) <- colnames(df)
            rownames(df.out)[nrow(df.out)] <- 'target'
            return(df.out)
        }
            
        ## add list of duplicates to dataframe column
        ## need to add to df, not df.out, because row corresponds to df
        ## purpose of df.out is only to facilitate search
        ## if (row == 3) browser()
        df$duplicate[row] <- list(dup.vec.df)
        
        ## reduce dataframe by removing duplicates as they are found
        ## if (row == 4) browser()
        if (length(dup.vec) > 1) {
            ## more than just the row being evaluated was ientified as a duplicate
            ## remove the row being evaluated from the list of rows to be removed
            ##    dups.to.remove <- dup.vec[-row]
            ## better yet, remove row being evaluated and earlier rows from list to be removed
            dups.to.remove <- dup.vec.df[dup.vec.df>row]
            ## find location in df.out
            loc <- which(rownames(df.out) %in% dups.to.remove)
            ## eliminate remainng rows
            if (length(loc) > 0) df.out <- df.out[-loc,]
            ## add dups.to.remove to growing list for later use with df
            dups.to.remove.all <- c(dups.to.remove.all, dups.to.remove)
        }
        ## cat('\n')
        ## cat('row ', row, 'has duplicates:', dups.to.remove, '\n')
        ## print(df.out)
  
    }
    
    ## remove initial NA from dups.to.remove.all
    dups.to.remove.all <- as.numeric(dups.to.remove.all[-1])
    ## remove duplicates from df
    df.kept    <- df[-dups.to.remove.all,]
    df.removed <- df[ dups.to.remove.all,]
    return(list(df.kept=df.kept, df.removed=df.removed))
    
}

df.duplicates_test <- function(tol.fraction=0.1, tol.absolute=3) {

    ##--------------------------------
    cat('\n\n##--------------------------------\n')
    cat('original dataframe\n')
    df <- data.frame(aa=c(    1,    1,    1,    2),
                     bb=c(  100,  101,   99,  110),
                     cc=c( 1000,  999, 1010, 1100),
                     dd=c(10000, 9999, 9999, 9999),
                     ee=c(  100,  102,  100,  100))
    print(df)
    
    
    ##--------------------------------
    cat('\n##--------------------------------\n')

    tol <- tol.fraction

    cat('\nidentify duplicates based on fractional tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol, tol.type='fraction')
    print(duplicates)

    cat('\nfind duplicates of provided vector based on fractional tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol, vector=c(aa=1,bb=100,cc=1008,dd=10000,ee=100),
                                tol.type='fraction')
    print(duplicates)


    ##--------------------------------
    cat('\n##--------------------------------\n')
    
    tol <- tol.absolute
    
    cat('\nidentify duplicates based on absolute tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol, tol.type='absolute')
    print(duplicates)

    cat('\nfind duplicates of provided vector based on absolute tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol, vector=c(aa=1,bb=100,cc=1008,dd=10000,ee=100),
                                tol.type='absolute')
    print(duplicates)

    cat('\nfor specified parameters, find duplicates of provided vector based on absolute tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol,
                                param=c('aa', 'dd'),
                                vector=c(dd=10000,aa=1,bb=100,ee=100),
                                tol.type='absolute')
    print(duplicates)

}
## df.duplicates_test(0.1, 3)
