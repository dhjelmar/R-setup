plotspace <- function(rows,columns) par(mfrow=c(rows,columns))

#----------------------------------------------------------------------------

# printdf(df,nrows,var) where var is a vector of variable names.
# var is a vector list of all variables in dataframe, df, to be printed.
# Each variable needs to be in quotes, e.g.:
#     printdf(df,18,c('ID','Work Item Type','Effort'))
