droplevels.all <- function(df) {
    ## drop all levels from dataframe, df
    df.num <- df
    for (i in 1:ncol(df)) {
        df.num[,i] <- as.numeric(as.character(df.num[,i]))
    }
    ## df.num
    df.sum <- colSums(df.num)
    for (i in 1:ncol(df)) {
        if (isFALSE(is.na(df.sum[i]))) {
            ## column is all numeric
            df[,i] <- df.num[,i]
        } else {
            ## column is not all numeric
            df[,i] <- as.character(df[,i])
        }
    }
    return(df)
}
