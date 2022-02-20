# optim not converging

## consider the following dataset
x <- iris$Sepal.Width
mean(x)  # 3.057333
sd(x)    # 0.4358663
alpha <- 0.01
P     <- 0.99  # proportion or coverage
sided <- 1
tol_out_norm <- tolerance::normtol.int(x, alpha = alpha, P=P, side=sided)
upper_tol <- tol_out_norm$'1-sided.upper'
upper_tol

mle <- function(data, param, fit='n', alpha=0.01, P=0.99, sided=1, plots=FALSE, debug=FALSE) {

    x <- data

    out <- NULL
    
    if (grepl('n', fit)) {

        ## determine best fit using nnl
        nll.normal <- function(data, param, debug=FALSE){
            ## calculate nll (negative log likelihhod) for normal distribution
            x       <- data
            xbar    <- param[[1]]
            sdev    <- param[[2]]
            z       <- (x - xbar) / sdev
            pdf     <- 1/(sqrt(2 * pi) * sdev) * exp(-z^2 / 2) 
            ## the following is equivalent
            ## pdf     <- dnorm(x, mean = xbar, sd = sdev, log = FALSE)
            nll     <- -sum(log(pdf))
            if (isTRUE(debug)) cat('xbar=', signif(xbar,11), 'sdev=', signif(sdev,11), 'nll=', signif(nll,11), "\n")
            return(nll)
        }        
        out.bestfit <- optim(par     = list(xbar=0, sdev=1), 
                             fn      = nll.normal, 
                             data    = x,
                             debug   = debug,
                             control = list(trace=TRUE),
                             method  = "BFGS")
        nll.max.bestfit <- out.bestfit$value
        xbar <- out.bestfit$par[[1]]
        sdev <- out.bestfit$par[[2]]
        params <- list(xbar=xbar, sdev=sdev)

        ## redefine nnl function to fit on desired quantile
        nll.normal.q <- function(data, param, P, debug=FALSE){
            ## calculate nll (negative log likelihhod) for normal distribution
            x       <- data
            quant   <- param[[1]]  # substituted for xbar
            sdev    <- param[[2]]
            ## the following does not work well
            ## xbar    <- quant - sqrt( 2*sdev^2 * (-log(qnorm(P)) - log(sqrt(2*pi)*sdev)) )      
            ## repalced it using an alternate form of the normal PDF found here:
            ## F(x) = 0.5 *( 1 + erf( (x-xbar)/(sdev * sqrt(2))))
            xbar <- quant - sqrt(2) * sdev * pracma::erfinv(2*P-1)
            z       <- (x - xbar) / sdev
            pdf     <- 1/(sqrt(2 * pi) * sdev) * exp(-z^2 / 2) 
            nll     <- -sum(log(pdf))
            if (isTRUE(debug)) cat('quant=', signif(quant,11), 'sdev=', signif(sdev,11), 'nll=', signif(nll,11), "\n")
            return(nll)
        }

        ##---------------------- aside
        ## f <- function(x, xbar, sdev) {
        ##   1/(sqrt(2*pi)*sdev) * exp(-(x-xbar)^2/(2*sdev^2))
        ## }
        ## f(xbar, xbar, sdev)
        ## par(mfrow=c(2,2))
        ## curve(f(x, xbar, sdev), xbar-3*sdev, xbar+3*sdev, xlab='x=observation', ylab='PDF')
        ## curve(dnorm(x, xbar, sdev), xbar-3*sdev, xbar+3*sdev, xlab='x=observation') # dnorm(x, mean, sd)
        ## curve(pnorm(x, xbar, sdev), 1, 5, xlab='x=quantile')                # pnorm(q, mean, sd)
        ## curve(qnorm(x, xbar, sdev), 0, 1, xlab='x=probability')                     # qnorm(p, mean, sd)
        ##---------------------- aside
        
        
        ## find confidence limit at level alpha for requested coverage, P
        ## confidence limit: P=0.5
        ## tolerance limit:  P=0.99
        quant.coverage <- qnorm(P, xbar, sdev)   # for p=0.5, this should be xbar
        quant.param <- c(quant=quant.coverage, sdev=sdev)

        if (isTRUE(plots)) {
            ## plot the likelihood as a function of the quantile
            quant.min <- xbar-4*sdev
            quant.max <- xbar+4*sdev
            quant     <- seq(quant.min, quant.max, length.out=100)
            ll <- NA
            for (i in 1:100) {
                ## vary quant parameter to see impact on likelihood
                param.vary.q <- list(quant  = quant[i],
                                     sdev   = sdev)
                ll[i] <- -nll.normal.q(x, param.vary.q, P=P, debug=FALSE)
            }
            plot(quant, ll, xlab='quantile', ylab='Log Likelihood')
            abline(v=xbar)
        }

        out.bestfit.q <- optim(par     = quant.param, 
                               fn      = nll.normal.q, 
                               data    = x,
                               P       = P,
                               debug   = debug,
                               control = list(trace=TRUE),
                               method  = "BFGS")
        nll.max.bestfit.q <- out.bestfit.q$value
        quant.conf <- out.bestfit.q$par[[1]]
        sdev.conf  <- out.bestfit.q$par[[2]]
        
        
        ## calculate confidence limits using LR (Likelihood Ratio)
        ## confidene limit is defined at likelihood that is lower than max by chi-squared
        ll.max.P <- -nll.max.bestfit
        ll.tol <-  ll.max.P - qchisq(1 - alpha/sided, 1)   # qchisq(1-0.01/1, 1) = 6.634897 
        if (isTRUE(plots)) {
            abline(h=ll.tol)
        }


        ## conf <- data.frame(alpha, P, sided, conf.lower, xbar, conf.upper)
        conf <- NA
        
    }

    return(list(out.bestit=out.bestfit, conf=conf))
}
out <- mle(x, c(mean=0, sd=1), fit='n', alpha=0.01, P=0.5 , sided=1, plots=TRUE, debug=FALSE)
out <- mle(x, c(mean=0, sd=1), fit='n', alpha=0.01, P=0.99, sided=1, plots=TRUE, debug=FALSE)
