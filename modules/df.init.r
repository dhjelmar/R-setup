df.init <- function(columns) {
    ## initialize empty dataframe with named columns
    ## columns = vector of column names
    df <- data.frame(matrix(nrow=0, ncol=length(columns)))
    colnames(df) <- columns
    return(df)
}
