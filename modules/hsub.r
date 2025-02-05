hsub <- function(h, p) {
  # returns subcooling enthalpy given enthalpy and pressure
  #
  # given:  h    = enthalpy in Btu/lbm
  #         p    = pressure in psia
  #
  # return: hsub = subcooling in Btu/lbm
  
  Tsat_K <- unlist( purrr::pmap(list(p = p / 145.038),
                                function(p) IAPWS95::TSatp(p)))
  hsat   <- unlist( purrr::pmap(list(T = Tsat_K),
                                function(T) IAPWS95::hfT(T) * 0.429923))
  hsub <- hsat - h
  return(hsub)
}

hsub_test <- function() {
  df <- data.frame(h=c(170, 400, 600), p=c(1000, 2000, 2000))
  df$hsub <- hsub(df$h, df$p)
  return(df)
}
hsub_test()
