'%ino%' <- function(x, table) {
    ## operator that preserves matching order unlike %in%
    ## credit to https://stackoverflow.com/questions/10586652/r-preserve-order-when-using-matching-operators-in
    ## can use similar to %in% but do not need to put it in which()
    xSeq <- seq(along = x)
    names(xSeq) <- x
    Out <- xSeq[as.character(table)]
    Out[!is.na(Out)]
}
## both of the following return the positions of the 2nd vector as a vector
##
## the %in% operator does not preserve the orer and returns a simple vector
## which(c('a', 'b', 'c', 'd', 'e') %in% c('d', 'b'))
## # returns: [1] 2 4
##
## the %ino% operator does preseerve the order and returns a column vector with row names d and b
## c('a', 'b', 'c', 'd', 'e') %ino% c('d', 'b')
## # returns: d b 
## #          4 2
