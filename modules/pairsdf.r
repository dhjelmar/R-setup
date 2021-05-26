pairsdf <- function(df,var) {
  # usage: pairsdf(df,c('power_ratio','gin_fb','hsub_in','beta'))
  dfcompare <- subset(df,select=var)
  pairs(dfcompare)
}

#----------------------------------------------------------------------------

