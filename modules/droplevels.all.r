droplevels.all <- function(df) {
    ## drop all levels from dataframe, df
    fac_cols <- sapply(df, is.factor)                   # Identify all factor columns
    df[fac_cols] <- lapply(df[fac_cols], as.character)  # Convert all factors to characters
    return(df)
}


droplevels.all.test <- function() {
    df <- data.frame(x1 = letters[1:5],                           # Create example data
                     x2 = 1:5,
                     x3 = "XXX",
                     x4 = c('one', 'two', 'three', 'four', 'five'),
                     x5 = c(1.0, 2, 3, 4, 5.2))
    df$x4 <- as.character(df$x4)
    print(as_tibble(df))
    cat('\n')
    out <- droplevels.all(df)
    print(as_tibble(out))
    cat('\n')
    print(out[,4])
}
