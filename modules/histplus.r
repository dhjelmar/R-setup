histplus <- function(x, breaks=NULL, distribution='normal', sided=2, alpha=0.05, pvalue=0.95, nsigma=2) {
    ## plot histogram and normal distribution
    ## distribution = normal or Weibull
    ## sided  = -1 (lower), 1 (upper), or 2 sided interval
    ## alpha  = 1 - confidence (used in confidence, prediction, and tolerance intervals)
    confidence <- 1-alpha
    ## pvalue = % of data to be covered (used in tolerance interval)
    ## nsigma = number of standard deviations for mean +/- nsigma * standard deviation

    #### debug
    ## x <- mtcars$mpg
    ## distribution <- 'normal'
    ## sided  <- 2
    ## alpha  <- 0.05
    ## confidence <- 1-alpha
    ## pvalue <- 0.95
    ## nsigma <- 2
    ## breaks <- NULL

    ## calculate limits
    if (distribution=='Weibull') {
        ## not sure how to figure out other intervals for Weibull since mean and sd do not make sense
        tol_out <-  exttol.int(x, alpha = alpha, P=pvalue, side=abs(sided))
        shape <- tol_out$'shape.1'
        scale <- tol_out$'shape.2'
        if (sided==-1) {
            lower_conf_limit      <- NA   # not sure how to do this
            lower_sigma           <- NA   # N/A                    
            lower_pred_limit      <- NA   # not sure how to do this
            lower_tolerance_limit <- tol_out$'1-sided.lower'
            upper_conf_limit      <- NA   # N/A                    
            upper_sigma           <- NA   # N/A                    
            upper_pred_limit      <- NA   # N/A                    
            upper_tolerance_limit <- NA   # N/A                    
        } else if (sided==1) {
            lower_conf_limit      <- NA   # N/A                    
            lower_sigma           <- NA   # N/A                    
            lower_pred_limit      <- NA   # N/A                    
            lower_tolerance_limit <- NA   # N/A                    
            upper_conf_limit      <- NA   # not sure how to do this
            upper_sigma           <- NA   # N/A                    
            upper_pred_limit      <- NA   # not sure how to do this
            upper_tolerance_limit <- tol_out$'1-sided.upper'
        } else {                         
            lower_conf_limit      <- NA   # not sure how to do this
            lower_sigma           <- NA   # N/A                    
            lower_pred_limit      <- NA   # not sure how to do this
            lower_tolerance_limit <- tol_out$'2-sided.lower'
            upper_conf_limit      <- NA   # not sure how to do this
            upper_sigma           <- NA   # N/A                    
            upper_pred_limit      <- NA   # not sure how to do this
            upper_tolerance_limit <- tol_out$'2-sided.upper'
        }                                
    } else {  # normal distribution (default)
        xmean <- mean(x)
        xsd   <- sd(x)
        tol_out <- normtol.int(x, alpha = alpha, P=pvalue, side=abs(sided))
        if (sided==-1) {
            student_t <- abs(qt(alpha, length(x)-1))
            lower_conf_limit      <- xmean - student_t*xsd*(  1/length(x))^0.5
            lower_sigma           <- xmean -    nsigma*xsd
            lower_pred_limit      <- xmean - student_t*xsd*(1+1/length(x))^0.5
            lower_tolerance_limit <- tol_out$'1-sided.lower'
            upper_conf_limit      <- NA
            upper_sigma           <- NA
            upper_pred_limit      <- NA
            upper_tolerance_limit <- NA
        } else if (sided==1) {
            student_t <- abs(qt(alpha, length(x)-1))
            lower_conf_limit      <- NA
            lower_sigma           <- NA
            lower_pred_limit      <- NA
            lower_tolerance_limit <- NA
            upper_conf_limit      <- xmean + student_t*xsd*(  1/length(x))^0.5
            upper_sigma           <- xmean +    nsigma*xsd
            upper_pred_limit      <- xmean + student_t*xsd*(1+1/length(x))^0.5
            upper_tolerance_limit <- tol_out$'1-sided.upper'
        } else {
            student_t <- abs(qt(alpha/2, length(x)-1))
            lower_conf_limit      <- xmean - student_t*xsd*(  1/length(x))^0.5
            lower_sigma           <- xmean -    nsigma*xsd
            lower_pred_limit      <- xmean - student_t*xsd*(1+1/length(x))^0.5
            lower_tolerance_limit <- tol_out$'2-sided.lower'
            upper_conf_limit      <- xmean + student_t*xsd*(  1/length(x))^0.5
            upper_sigma           <- xmean +    nsigma*xsd
            upper_pred_limit      <- xmean + student_t*xsd*(1+1/length(x))^0.5
            upper_tolerance_limit <- tol_out$'2-sided.upper'
        }
    }        
    
    ## generate distribution
    xmin <- min(x,lower_sigma,lower_pred_limit,lower_tolerance_limit, na.rm=TRUE)
    xmax <- max(x,upper_sigma,upper_pred_limit,upper_tolerance_limit, na.rm=TRUE)
    chuncks <- (xmax-xmin)/1000
    xhist   <- seq(xmin,xmax,by=chuncks)
    if (distribution=='Weibull') {
        xdensity <- dweibull(xhist,shape=shape,scale=scale)
    } else {
        xdensity <- dnorm(xhist,xmean,xsd)
    }
    xdensitymax <- max(xdensity)

    ## make histogram
    if (is.null(breaks)) {
        ##  hist(x, xlim=c(xmin,xmax), freq=FALSE)  # specifying xlim messes up axes
        hist(x,                freq=FALSE)
    } else {
        hist(x, breaks=breaks, freq=FALSE)
    }
    ## add distribution
    lines(x=xhist, y=xdensity, col='red', lty=1)
    ## add lines for mean and limits
    if (distribution=='normal') abline(v=xmean,col="red")
    abline(v=lower_conf_limit          ,col="blue")
    abline(v=lower_sigma               ,col="darkviolet", lty=2)
    abline(v=lower_pred_limit          ,col="darkviolet")
    abline(v=lower_tolerance_limit     ,col="black")
    abline(v=upper_conf_limit          ,col="blue")
    abline(v=upper_sigma               ,col="darkviolet", lty=2)
    abline(v=upper_pred_limit          ,col="darkviolet")
    abline(v=upper_tolerance_limit     ,col="black")
    cat(sided,"sided", distribution,"distribution intervals\n")
    cat("nsigma                 =",nsigma,"\n")
    cat("confidence = 1 - alpha =",confidence,"\n")
    cat("coverage   = pvalue    =",pvalue,"\n")
    if (distribution=='Weibull') {
        cat("shape =",shape,"\n")
        cat("scale =",scale,"\n")
        interval <- c("confidence", "prediction", "tolerance")
        df <- data.frame(interval)
        df$lower <- c(lower_conf_limit, lower_pred_limit, lower_tolerance_limit)
        df$upper <- c(upper_conf_limit, upper_pred_limit, upper_tolerance_limit)
    } else {     # assume normal
        cat("mean                   =",xmean,"\n")
        cat("standard deviation     =",xsd,"\n")
        interval <- c("confidence", "mean +/- nsigma * sd", "prediction", "tolerance")
        df <- data.frame(interval)
        df$lower <- c(lower_conf_limit, lower_sigma, lower_pred_limit, lower_tolerance_limit)
        df$upper <- c(upper_conf_limit, upper_sigma, upper_pred_limit, upper_tolerance_limit)
        ##    cat("mean                                      =", xmean, "\n")
        ##    cat(sided,  "sided upper", confidence, "confidence limit       =", upper_conf_limit,"\n")
        ##    cat(nsigma, "sigma upper limit                       =", upper_sigma     , "\n")
        ##    cat(sided,  "sided", confidence,"upper prediction limit       =", upper_pred_limit, "\n")
        ##    cat(sided,  "sided", confidence,"/",pvalue,"upper tolerance limit =", upper_tolerance_limit, "\n")
    }
    print(df)
}
# histplus(mtcars$mpg)

#--------------------------------

