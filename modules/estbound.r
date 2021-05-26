estbound <- function(model) {
  # summary(model) prints estimates which can be accessed as model$coefficient[#]
  # confint(model) prints 95% upper / lower bounds on estimates
  # econfint(model) prints estimates and 95% upper / lower bounds on estimates
  bounds <- confint(model)
  estbound <- cbind( c(model$coefficient),
                     bounds[,1:ncol(bounds)],
                     summary(model)$coefficients[,"Pr(>|t|)"]
  )
  colnames(estbound)[1] <- "estimates"
  colnames(estbound)[4] <- "Pr(>|t|)"
  estbound
}


#-----------------------------

