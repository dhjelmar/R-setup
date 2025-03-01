# data structures example

################################################################################
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

################################################################################
# create dataframe from vectors where each vector is a column
# BEST way to create a dataframe
df <- data.frame(v1, v2, v3)
df

# create empty dataframe then populate row by row
# BEST way to not have to redo column types
source('~/Documents/GitHub/R-setup/modules/df.init.r')
df <- df.init(rows=0, columns=c('a', 'b', 'c', 'd')) # initialize data frame
df
df[1,] <- list(1, 'a', 'b', 10)
df[2,] <- list(2, 'c', 'd', 20)
df[3,] <- list('a', 'c', 'd', 20)
df
typeof(df[1,4])   # double because list preserved numbers as numbers

###########################
# OTHER WAYS TO CREATE DATAFRAMES

# create dataframe from vectors where each vector is a row
# reasonable way to create dataframe but may need to fix column types
df <- NULL
df <- data.frame(rbind(v1, v2, v3))
names(df) <- c('a', 'b', 'c', 'd')                                   # rename columns
row.names(df) <- seq(1,nrow(df))                                     # replace row names with numbers
df[, c('a', 'c', 'd')] <- lapply(df[, c('a', 'c', 'd')], as.numeric) # fix column types
df                  # looks great
typeof(df[2,1])     # double
typeof(df[2,2])     # character

# I thought the following might work and the df looks good but the types are odd
L1 <- list(1.1, 'a', 3, 4)    # all character
L2 <- list(1, 2.01, 3, 4)     # all double but shown as 1.00 2.01 3.00 4.00
L3 <- list(1, 2, 3, 4)        # all double but shown as 1 2 3 4
L4 <- list(6, 7, 8, 9)
df1 <- data.frame(rbind(L1, L2, L3))
df1                # looks great but...
typeof(df1[1,1])   # type 'list'??????
typeof(df1[1,4])   # type 'list'??????

# create dataframe from list of vectors for each column
my_list <- list(1:4,
                c('koala', 'hedgehog', 'sloth', 'panda'), 
                c('Australia', 'Italy', 'Peru', 'China'),
                c(21, 18, 17, 10))
df <- data.frame(my_list)
names(df) <- c('a', 'b', 'c', 'd')
df

################################################################################

# dictionary
sounds <- c(cat="meow", dog="woof", horse="neigh")
sounds['cat']

# dictionary with more structure
# install.packages('Dict')
ages <- Dict::Dict$new(
            Charlie = 40L,
            Alice = 30L,
            Bob = 25L,
            .class = "integer",
            .overwrite = TRUE
)
ages["Bob"]   # returns 25
ages$keys     # returns 'Charlie' 'Alice' 'Bob
ages$values
ages$items
# to add an item
ages['Fred'] <- 30L
ages$keys
ages['Fred']
