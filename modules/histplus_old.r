histplus_old <- function(x, breaks=NULL, sided=2, tol_coverage=0.95, confidence=0.95, nsigma=2) {
    ## HISTOGRAM WITH MEAN AND LIMITS (confidence, #sigma, prediction, tolerance)
    ## plot histogram and normal distribution

    #### debug
    ##x <- mtcars$mpg
    ##sided  <- 2
    ##tol_coverage <- 0.95
    ##confidence <- 0.95
    ##nsigma <- 2
    ##breaks  <- NULL

    ## set paramters
    pvalue <- tol_coverage
    alpha  <- 1-confidence
    if (sided==1) {
        alpha <- alpha/2    # also need to halve pvalue?                       # dlh to check
    }
    xmean <- mean(x)
    xsd   <- sd(x)
    
    ## calculate limits
    student_t <- abs(qt(alpha*sided, length(x)-1))
    upper_conf_limit <- xmean + student_t*xsd*(  1/length(x))^0.5
    upper_sigma      <- xmean +    nsigma*xsd
    upper_pred_limit <- xmean + student_t*xsd*(1+1/length(x))^0.5
#    tol_out_normal <- normtol.int(x, alpha = alpha, P=pvalue, side=sided)     # dlh to fix
#    upper_tolerance_limit_norm <- tol_out$'1-sided.upper'                     # dlh to fix

    ## generate distribution
    xhist <- seq(min(x),max(x,upper_sigma,upper_pred_limit),by=0.001)          # dlh add upper tolerance limit to this when fix it
    xdensity <- dnorm(xhist,xmean,xsd)                                         # dlh add option to input distribution?
    xdensitymax <- max(xdensity)

    ## make histogram
    if (is.null(breaks)) {
        hist(x,                freq=FALSE)
    } else {
        hist(x, breaks=breaks, freq=FALSE)
    }
    ## add distribution
    lines(xhist, xdensity, col='red', lty=1)
    ## add lines for mean and limits
    abline(v=xmean,col="red")
    abline(v=upper_conf_limit          ,col="blue")
    abline(v=upper_sigma               ,col="darkviolet", lty=2)
    abline(v=upper_pred_limit          ,col="darkviolet")
#    abline(v=upper_tolerance_limit_norm,col="black")                          # dlh to fix
    cat("mean                                      =", xmean, "\n")
    cat(sided,  "sided upper", confidence, "confidence limit       =", upper_conf_limit,"\n")
    cat(nsigma, "sigma upper limit                       =", upper_sigma     , "\n")
    cat(sided,  "sided", confidence,"upper prediction limit       =", upper_pred_limit, "\n")
    upper_tolerance_limit <- NULL                                              # dlh to fix
    cat(sided,  "sided", confidence,"/",pvalue,"upper tolerance limit =", upper_tolerance_limit, "\n")
}

