df.duplicates <- function(df, tol=0.01, tol.type='fraction', remove=TRUE) {
    ## input: df       = numeric, 2D dataframe
    ##        tol      = tolerance
    ##        tol.type = 'fraction' indicates tolerance is the 
    ##                   specified fraction of value being compared 
    ##                   (i.e., tolerance = tol * df[i,j])
    ##                 = 'absolute' indicates tolerance is the
    ##                    actual value of the tolerance
    ##                    (i.e., tolerance = tol)
    ##        remove   = TRUE (default) removes duplicates from df
    ##                 = FALSE just identifies duplicates in a new "duplicate" parameter
    
    df.out <- df
    if (isFALSE(remove)) df.out$duplicate <- NA
    for (row in 1:nrow(df)) {
        ## define vector, vec, to be compared to each row of df
        vec <- df[row,]
        
        if (isFALSE(remove)) {
            out <- data.frame(matrix(nrow=ncol(df), ncol=nrow(df)))
        } else {
            out <- data.frame(matrix(nrow=ncol(df.out), ncol=nrow(df.out)))
        }
        
        for (col in 1:ncol(df)) {
            ## use lapply to look at every row for same value (within tol) in each column
            if (isFALSE(remove)) {
                if (tol.type == 'fraction') {
                    out[col,] <- unlist(lapply(df[,col], function(x) dplyr::near(x, vec[col], tol=tol*x)))
                } else {
                    out[col,] <- unlist(lapply(df[,col], function(x) dplyr::near(x, vec[col], tol=tol)))
                }
            } else {
                ## insetad of above, look in shrinking df.nodup as duplicates are removed
                if (tol.type == 'fraction') {
                    out[col,] <- unlist(lapply(df.out[,col], function(x) dplyr::near(x, vec[col], tol=tol*x)))
                } else {
                    out[col,] <- unlist(lapply(df.out[,col], function(x) dplyr::near(x, vec[col], tol=tol)))
                }
            }
        }  
        
        ## duplcates are any rows in df that are the same in every column
        duplicates <- list(which(apply(out, 2, all)))    # all tests entire vector for true
        dup.vec <- as.vector(duplicates[[1]])
        
        if (isFALSE(remove)) {
            ## add list of duplicates to dataframe column
            df.out$duplicate[row] <- list(dup.vec)
            
        } else {
            ## reduce dataframe by removing duplicates as they are found
            ## if (row == 4) browser()
            if (length(dup.vec) > 1) {
                ## more than just the row being evaluated was ientified as a duplicate
                ## remove the row being evaluated from the list of rows to be removed
                ##    dups.to.remove <- dup.vec[-row]
                ## better yet, remove row being evaluated and earlier rows from list to be removed
                dups.to.remove <- dup.vec[dup.vec>row]
                ## find location in df.out
                loc <- which(rownames(df.out) %in% dups.to.remove)
                ## eliminate remainng rows
                if (length(loc) > 0) df.out <- df.out[-loc,]
            }
            ## cat('\n')
            ## cat('row ', row, 'has duplicates:', dups.to.remove, '\n')
            ## print(df.out)
        }
        
    }
    
    return(df.out)
    
}

df.duplicates_test <- function() {
    cat('original dataframe\n')
    df <- data.frame(aa=c(    1,    1,    1,    2),
                     bb=c(  100,  101,   99,  110),
                     cc=c( 1000,  999, 1010, 1100),
                     dd=c(10000, 9999, 9999, 9999),
                     ee=c(  100,  102,  100,  100))
    print(df)
    
    tol=0.1
    tol=3
    
    cat('\nidentify duplicates based on fractional tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol, tol.type='fraction', remove=FALSE)
    print(duplicates)
    
    cat('\nremove duplicates based on fractional tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol, tol.type='fraction', remove=TRUE)
    print(duplicates)
    
    cat('\nidentify duplicates based on absolute tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol, tol.type='absolute', remove=FALSE)
    print(duplicates)
    
    cat('\nremove duplicates based on absolute tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol, tol.type='absolute', remove=TRUE)
    print(duplicates)
}

