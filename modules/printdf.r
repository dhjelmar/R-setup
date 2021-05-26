printdf <- function(df,nrows,var) {
  df <- df[,var]
  print(df,n=nrows)
}

## following is an alternate to the above
#printdf_merge <- function(df,nrows,var) {
#df <- df[,var]
##print(df,na.print=NULL,n=nrows)
#entries <- ncol(df) * nrows
#print(df,na.print=NULL,max=entries)
#}
##printdf(df,5,c('TS_ID','Spacing'))
