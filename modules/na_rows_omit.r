na_rows_omit <- function(df) {
    ## removes any row in dataframe that is all NA
    df <- df[rowSums(is.na(df)) != ncol(df),]
 }
