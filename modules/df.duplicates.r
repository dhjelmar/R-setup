df.duplicates <- function(df, tol=0.01, tol.type='fraction', target=NA, param=NA) {
    ## input:   df       = numeric, 2D dataframe
    ##          tol      = tolerance
    ##          tol.type = 'fraction' indicates tolerance is the 
    ##                     specified fraction of value being compared 
    ##                     (i.e., tolerance = tol * max(target + df[i,j], 1E-6)
    ##                   = 'absolute' indicates tolerance is the
    ##                      actual value of the tolerance on df and target
    ##                      (i.e., tolerance = 2*tol)

    ## options: target   = NA (default) looks for duplicates of all rows in df
    ##                   = vector looks for duplicates only of supplied target
    ##          param    = vector of parameters to use from within df and target

    ## renumber df
    rownames(df) <- 1:nrow(df)

    
    ##-----------------------------------------------------------------------------
    ## if target is entered as vector, convert to dataframe
    if (is.na(target[1])) {
        ## target is not specified
        target.specified <- FALSE
    } else {
        ## target is specified
        target.specified <- TRUE
        if (!is.data.frame(target)) {
            ## convert to dataframe
            target <- as.data.frame(t(target))
        }
        ## convert any parameters of type "factor" to "character"
        target <- droplevels.all(target)
    }

    
    ##-----------------------------------------------------------------------------
    ## determine which columns of df and, if specified, target to use in search
    if (!is.na(param[1])) {
        ## param is used to determine which columns of df and target to use
        df.cols  <- names(df) %ino% param
        if (!is.na(target[1])) {
            target.cols <- names(target) %ino% param
        } else {
            target.cols <- df.cols
        }
            
    } else if (!is.na(target[1])) {
        ## target is used to determine which columns of df to use (all target columns used)
        df.cols  <- names(df) %ino% names(target)
        target.cols <- 1:length(target)

    } else {
        ## only df is specified, so use all parameters
        df.cols <- 1:ncol(df)
        target.cols <- df.cols
    }
    nparam <- length(df.cols)
        
    ##-----------------------------------------------------------------------------
    ## search for near duplicates
    df.out <- df
    df$duplicate <- NA
    dups.to.remove.all <- NA
    for (row in 1:nrow(df)) {
        
        ## define target to be compared to each row of df
        if (isFALSE(target.specified)) {
            target <- droplevels.all(df[row,])
        }
        
        out <- data.frame(matrix(nrow=nrow(df.out), ncol=nparam))
        
        for (i in 1:nparam) {
            ## use lapply to look at every row in shrinking df.out for same value (within tol) in each column
            target.col.i     <- target[[target.cols[i]]]
            df.out.col.i  <- df.out[,df.cols[i]]
            if (is.character(target.col.i)) {
                ## need identical match because value is a character
                out[,i] <- unlist(lapply(df.out.col.i, function(x) x == target.col.i))
            } else if (tol.type == 'fraction') {
                ## use fractional tolerance
                tol.col.i <- tol*target.col.i
                out[,i] <- unlist(lapply(df.out.col.i, function(x) dplyr::near(x, target.col.i,
                                                                               tol=tol*max(x+target.col.i, 1E-6))))
            } else {
                ## use absolute tolerance
                tol.col.i <- tol
                out[,i] <- dplyr::near(df.out.col.i, target.col.i, tol=2*tol.col.i)
            }
        }  
        
        ## duplcates are any rows in df that are the same in every column
        duplicates <- list(which(apply(out, 1, function(x) all(x, na.rm=TRUE) )))    # all tests entire target for true
        dup.target <- as.vector(duplicates[[1]])
        ## above is the row(s) in df.out; what row(s) is dup.target in df?
        dup.target.df <- as.numeric(rownames(df.out[dup.target,]))

        if (isTRUE(target.specified)) {
            ## only need this one time through
            df.out <- df.out[dup.target,]
            if (!is.na(param[1])) target <- subset(target, select=param)
            df.out <- droplevels.all(df.out)
            df.out <- fastmerge(df.out, target)
            ## colnames(df.out) <- colnames(df)
            rownames(df.out)[nrow(df.out)] <- 'target'
            return(df.out)
        }
            
        ## add list of duplicates to dataframe column
        ## need to add to df, not df.out, because row corresponds to df
        ## purpose of df.out is only to facilitate search
        ## if (row == 3) browser()
        df$duplicate[row] <- list(dup.target.df)
        
        ## reduce dataframe by removing duplicates as they are found
        ## if (row == 4) browser()
        if (length(dup.target) > 1) {
            ## more than just the row being evaluated was ientified as a duplicate
            ## remove the row being evaluated from the list of rows to be removed
            ##    dups.to.remove <- dup.target[-row]
            ## better yet, remove row being evaluated and earlier rows from list to be removed
            dups.to.remove <- dup.target.df[dup.target.df>row]
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
    
    if (length(dups.to.remove.all) == 1) {   #  == 1 because 1st field will be NA
        ## no duplicates were found
        df.kept    <- df
        df.removed <- NA
    } else {
        ## remove initial NA from dups.to.remove.all
        dups.to.remove.all <- as.numeric(dups.to.remove.all[-1])
        ## remove duplicates from df
        df.kept    <- df[-dups.to.remove.all,]
        df.removed <- df[ dups.to.remove.all,]
    }
    
    return(list(df.kept=df.kept, df.removed=df.removed, df.all=df.out))
    
}

df.duplicates_test <- function(tol.fraction=0.1, tol.absolute=3) {

    ## expandFunctions::reset.warnings()
  
    ##--------------------------------
    cat('\n\n##--------------------------------\n')
    cat('original dataframe\n')
    df <- data.frame(aa=c(    1,    1,    1,    2),
                     bb=c(  100,  101,   99,  110),
                     cc=c( 1000,  999, 1010, 1100),
                     dd=c(10000, 9999, 9999, 9999),
                     char=c('one', 'two', 'one', 'two'),
                     ee=c(  100,  102,  100,  100))
    print(df)
    
    
    ##--------------------------------
    cat('\n##--------------------------------\n')
    ## source("F:\\Documents\\01_Dave\\Programs\\GitHub_home\\R-setup\\setup.r")
    ## tol.fraction=0.1; tol.absolute=3
    tol <- tol.fraction

    cat('\nidentify duplicates based on fractional tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol, tol.type='fraction')
    print(duplicates)

    cat('\nfind duplicates of provided target based on fractional tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol, target=c(aa=1,bb=100,cc=1008,dd=10000,ee=100),
                                tol.type='fraction')
    print(duplicates)


    ##--------------------------------
    cat('\n##--------------------------------\n')
    
    tol <- tol.absolute
    
    cat('\nidentify duplicates based on absolute tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol, tol.type='absolute')
    print(duplicates)

    cat('\nfind duplicates of provided target based on absolute tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol, target=c(aa=1,bb=100,cc=1008,dd=10000,ee=100),
                                tol.type='absolute')
    print(duplicates)

    cat('\nfor specified parameters, find duplicates of provided target based on absolute tolerance =', tol, '\n')
    duplicates <- df.duplicates(df, tol=tol,
                                param=c('aa', 'dd'),
                                target=c(dd=10000,aa=1,bb=100,ee=100),
                                tol.type='absolute')
    print(duplicates)

    cat('\nsame as above but with dataframe for target\n')
    duplicates <- df.duplicates(df, tol=tol,
                                param=c('aa', 'dd'),
                                target=data.frame(dd=10000,aa=1,bb=100,ee=100),
                                tol.type='absolute')
    print(duplicates)
    
}
## df.duplicates_test(0.1, 3)
