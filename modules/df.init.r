df.init <- function(rows, columns, value=NA) {
    ## initialize empty dataframe with named columns
    ## rows = number of rows
    ## columns = vector of column names
    ##         = number of columns
    if (is.numeric(columns) & length(columns) == 1) {
        df <- data.frame(matrix(data=value, nrow=rows, ncol=columns))
    } else {
        df <- data.frame(matrix(data=value, nrow=rows, ncol=length(columns)))
        colnames(df) <- columns
    }
    return(df)
}
