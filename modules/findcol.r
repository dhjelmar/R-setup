findcol <- function(mylist, target, index=FALSE) {
    '
    Identifies parameters or indices matching target substring in mylist.
    
    Input: mylist = list or can be a dataframe
                    (if a dataframe, will look in column names for substring)
           target = target substring
           index  = FALSE (default) returns values
                    TRUE returns location numbers of target
    Output: 
    '
    if (is.data.frame(mylist)) {
        mylist <- names(mylist)
    }
    location <- grep(target, mylist)
    if (isTRUE(index)) {
        return(location)
    } else {
        values <- mylist[location]
        return(values)
    }
}

#mylist <- c('asdfred', 2, 'jkl;', 'fred')
#findcol(mylist, 'jk')
#findcol(mylist, 2)
#findcol(mylist, 'l;', index=TRUE)
#findcol(mylist, 'fred')
#findcol(mylist, 'fred', index=TRUE)

