# data structures example

# vector parts are coerced to type that will work for all
v1 <- c(1.1, 'a', 3, 4)    # all character
typeof(v1[1])
v1
v2 <- c(1, 2.01, 3, 4)     # all double but shown as 1.00 2.01 3.00 4.00
v2
v3 <- c(1, 2, 3, 4)        # all double but shown as 1 2 3 4
v3
v4 <- 6:9
v4

# create dataframe from vectors where each vector is a column
df <- data.frame(v1, v2, v3)
df

# create dataframe from vectors where each vector is a row
df <- data.frame(rbind(v1, v2, v3))
df
typeof(df[1,4])     # character because, although column is all numbers, one of them came from v1
typeof(df[2,1])     # same here
names(df) <- c('a', 'b', 'c', 'd')  # rename columns
df
row.names(df) <- seq(1,nrow(df))    # replace row names with numbers
df
df[, c('a', 'c', 'd')] <- lapply(df[, c('a', 'c', 'd')], as.numeric)
typeof(df[1,4])     # double
typeof(df[2,1])     # double
typeof(df[2,2])     # character
df

# create dataframe from list of vectors
my_list <- list(1:4,
                c('koala', 'hedgehog', 'sloth', 'panda'), 
                c('Australia', 'Italy', 'Peru', 'China'),
                c(21, 18, 17, 10))
df <- data.frame(my_list)
names(df) <- c('a', 'b', 'c', 'd')
df

# create empty dataframe then popoulate row by row
df <- data.frame(rating=numeric(),
                 animal=character(),
                 country=character(),
                 avg_sleep_hours=numeric())
df[1,] <- c(1, 'a', 'b', 10)
df[2,] <- c(2, 'c', 'd', 20)
df
df[3,] <- c('a', 'c', 'd', 20) # specifying rating as numeric not useful
df                             # because it did not impact ability to add character

# dictionary
sounds <- c("cat"="meow", "dog"="woof", "horse"="neigh")
sounds["cat"]

sounds <- c(a="meow", b="woof", c="neigh")
sounds['a']


# dictionary with m ore structure
install.packages('Dict')
ages <- Dict::Dict$new(
    Charlie = 40L,
    Alice = 30L,
    Bob = 25L,
    .class = "integer",
    .overwrite = TRUE
)
ages["Bob"]


# create dataframe from list that was populated row by row
mylist <- list() #create an empty list
for (i in 1:5) {
    # create list of 5 vectors where each will eventually be a row in dataframe
    vec <- numeric(5) #preallocate a numeric vector
    for (j in 1:5) { #fill the vector
        vec[j] <- i^j 
    }
    mylist[[i]] <- vec #put all vectors in the list
}
mylist
mymatrix <- do.call("rbind",mylist) #combine all vectors into a matrix
mymatrix
df <- data.frame(mymatrix)          # assigns column names as X1, X2, X3...
df
names(df) <- c('a', 'b', 'c', 'd', 'e') # reassign column names
df

mylist_of_lists <- list() #create an empty list
mylist <- list()
mylist <- append(mylist, c('a', 'b', 1, 2))
mylist
mylist <- append(mylist, list('c', 'd', 3, 4))
mylist


for (i in 1:5) {
    # create list of 5 lists where each will eventually be a row in dataframe
    vec <- numeric(5) #preallocate a numeric vector
    for (j in 1:5) { #fill the vector
        vec[j] <- i^j 
    }
    mylist[[i]] <- vec #put all vectors in the list
}
mylist
mymatrix <- do.call("rbind",mylist) #combine all vectors into a matrix
mymatrix
df <- data.frame(mymatrix)          # assigns column names as X1, X2, X3...
df
names(df) <- c('a', 'b', 'c', 'd', 'e') # reassign column names
df

