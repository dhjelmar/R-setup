## https://www.r-bloggers.com/2019/08/maximum-likelihood-estimation-from-scratch/
data = data.frame(t = c(0, 16, 22, 29, 36, 58),
                  G = c(0, 0.12, 0.32, 0.6, 0.79, 1))
plot(data)

G = function(pars, t) {
  # Extract parameters of the model
  Gmax = pars[1]
  k = pars[2]
  th = pars[3]
  # Prediction of the model
  DG = Gmax/(1 - 1/(1 + exp(k*th)))
  Go = DG/(1 + exp(k*th))
  DG/(1 + exp(-k*(t - th))) - Go
}

## add function G to the curve for initial guesses of parameters
pars <- c(Gmax = 1, k = 0.15, th = 30)
curve(G(pars, x), 0, 60, add = TRUE)

## to instead estimate pars, use MLE

## first define NNL, negative log-likelihood
NLL <- function(pars, data) {
  # Values predicted by the model
  Gpred = G(pars, data$t)
  # Negative log-likelihood 
  -sum(dnorm(x = data$G, mean = Gpred, sd = pars[4], log = TRUE))
}

## guess parameters including a guess for standard deviation
par0 = c(pars, sd = 0.01)
## optimize model
fit = optim(par = par0, fn = NLL, data = data, control = list(parscale = abs(par0)), 
            hessian = TRUE)
fit$par
curve(G(fit$par[1:3], x), 0, 60, col='blue', add = TRUE)

##-----------------------------------------------------------------------------
## https://www.r-bloggers.com/2019/08/maximum-likelihood-estimation-from-scratch/
source('setup.r')

##----------------------
## example using normal distribution to fit
x <- rnorm(100)
x <- mtcars$mpg
norm.fit <- function(mu, sigma) {
    -sum(dnorm(x, mu, sigma, log=TRUE))
}
## x is the name of the variable containing the data to be fitted
mle.results <- bbmle::mle2(norm.fit, start=list(mu=1,sigma=1), data=list(x))  # fails
mle.results <- bbmle::mle2(norm.fit, start=list(mu=18,sigma=5), data=list(x)) # works
summary(mle.results)@coef

## cannot put x in the function call if using mle2
## the following fails
norm.fit.x <- function(x, mu, sigma) {
    if (missing(x)) 
    -sum(dnorm(x, mu, sigma, log=TRUE))
}
bbmle::mle2(norm.fit.x, start=list(mu=18,sigma=5), data=list(x)) # fails "x is missing"

## also cannot put x in the function call if using stats4::mle or maxLik::maxLik

## the following works (convergence = 0) and allows data in function call
## There are some warnings online about optim, but the above mle and mle2 packages
## seem to be wrappers on optim.
norm.fit <- function(data, fit) {
  x     <- data
  mu    <- fit[[1]]
  sigma <- fit[[2]]
  -sum(dnorm(x, mu, sigma, log=TRUE))
}
fit <- optim(par  = c(mu=18, sigma=5),  # parameters to be optimized and initial guess
                  fn   = norm.fit,           # NNL function
                  data = x,                  # data
                  method = 'L-BFGS-B',
                  control = list(parscale = c(mu=18, sigma=5)), # not clear this was needed
                  hessian = TRUE)
fit
fit$par
standard.error <- sqrt(diag(solve(fit$hessian)))   # matches bbmle::mle2 output
standard.error


##-------------------------------------------------------------------
##-------------------------------------------------------------------
## now try with JohnsonSU
plotspace(2,2)

## first get JohnsonSU parameters from R package fits and plot histogram
x <- mtcars$mpg
out <- ExtDist::eJohnsonSU(x)
jparms.ExtDist <- list(gamma   = out$gamma,
                       delta   = out$delta,
                       xi      = out$xi,
                       lambda  = out$lambda,
                       type    = 'SU')
## out <- hist_nwj(x, type='j', jfit=jparms.ExtDist)  # failure

jparms.SuppDists  <- SuppDists::JohnsonFit(x)
out <- hist_nwj(x, type='j', jfit=jparms.SuppDists)  # SU fit

## use SuppDists parameters for initial guesses below
jparms <- jparms.SuppDists


## try fitting JohnsonSU with bbmle::mle2()
nll.johnsonsu.mle2 <- function(gamma, delta, xi, lambda) {
    ## x is the name of the variable containing the data to be fitted
    ## x is defined in the bbmle::mle2 function
    # PDF for Johnson SU
    pdf <- delta /( lambda * sqrt(2 * pi)   ) *
        1 / sqrt(1 +            ( (x-xi)/lambda )^2)  *
        exp( -0.5*(gamma + delta * asinh( (x-xi)/lambda ))^2 )
    -sum(log(pdf))
    ## the following is equivalent to the above as long as it picks Johnson SU
    ## type <- 'SU'
    ## -sum(SuppDists::dJohnson(x, parms=c(gamma, delta, xi, lambda, type), log=TRUE))
}
out.mle2 <- bbmle::mle2(nll.johnsonsu.mle2, start=jparms[1:4], data=list(x))
out.mle2@coef
summary(out.mle2)@coef
jparms.mle2 <- list(gamma  = summary(out.mle2)@coef[1],
                   delta  = summary(out.mle2)@coef[2],
                   xi     = summary(out.mle2)@coef[3],
                   lambda = summary(out.mle2)@coef[4],
                   type   = 'SU')
curve(SuppDists::dJohnson(x, jparms.mle2), min(x), max(x, 70), col='red', add=TRUE)


##----------------------
## try fitting JohnsonSU  with optim()
nll.johnsonsu <- function(x, jparms, debug=FALSE) {
    ## x is the name of the variable containing the data to be fitted
    ## if use optim, cannot access as jparms$gamma because jparms is passed in as atomic vector
    gamma <- jparms[[1]]  
    delta <- jparms[[2]]
    xi    <- jparms[[3]]
    lambda <- jparms[[4]]
    ## PDF for Johnson SU
    pdf <- delta /( lambda * sqrt(2 * pi)   ) *
        1 / sqrt(1 +            ( (x-xi)/lambda )^2)  *
        exp( -0.5*(gamma + delta * asinh( (x-xi)/lambda ))^2 )
    -sum(log(pdf))
    ## the following is equivalent to the above as long as it picks Johnson SU
    ## type <- 'SU'
    ## -sum(SuppDists::dJohnson(x, parms=c(gamma, delta, xi, lambda, type), log=TRUE))
    if (isFALSE(debug)) {
        return( -sum(log(pdf)) )
    } else {
        return(list(gamma=gamma, delta=delta, xi=xi, lambda=lambda, pdf=pdf))
    }
}
out.optim <- optim(par = jparms[1:4], # parameters to be optimized and initial guess
                   fn  = nll.johnsonsu,     # NNL function
                   x   = x,                 # data
                   control = list(parscale = jparms[1:4]) # not clear this was needed
                   )
jparms.optim <- as.list(out.optim$par)
jparms.optim$type <- 'SU'
curve(SuppDists::dJohnson(x, jparms.optim), min(x), max(x, 70), col='blue', lty=2, add=TRUE)

## compare parameters
unlist(jparms[1:4])
unlist(jparms.mle2[1:4])
unlist(jparms.optim[1:4])


## compare qqplots
qqplot_nwj(x, type='j', jfit=jparms      , main='SuppDists auto fit')
qqplot_nwj(x, type='j', jfit=jparms.mle2 , main='MLE2 forced JohnsonSU fit')
qqplot_nwj(x, type='j', jfit=jparms.optim, main='optim forced JohnsonSU fit')

##------------------------------
## turn around johnson.fit for tolerance limit calculation
nll.johnsonsu.quant <- function(x, P, jparms, debug=FALSE) {
    ## x is the name of the variable containing the data to be fitted
    ## P is the percent coverage
    ## if use optim, cannot access as jparms$gamma because jparms is passed in as atomic vector
    ## replacing gamma with quant
    quant  <- jparms[[1]]  
    delta  <- jparms[[2]]
    xi     <- jparms[[3]]
    lambda <- jparms[[4]]
    ## write gamma as a function of quant, delta, xi and lambda
    gamma <- qnorm(P) - delta * asinh( (quant-xi)/lambda )
    ## PDF for Johnson SU
    pdf <- delta /( lambda * sqrt(2 * pi)   ) *
           1 / sqrt(1 +            ( (x-xi)/lambda )^2)  *
           exp( -0.5*(gamma + delta * asinh( (x-xi)/lambda ))^2 )
    if (isFALSE(debug)) {
        return( -sum(log(pdf)) )
    } else {
        return(list(quant=quant, gamma=gamma, delta=delta, xi=xi, lambda=lambda, pdf=pdf))
    }
}
quant.P <- ExtDist::qJohnsonSU(0.99, params = jparms.optim)
quant.parms <- c(quant=quant.P, jparms.optim[2:4])
out.optim.quant <- optim(par = quant.parms, # parameters to be optimized and initial guess
                         fn  = nll.johnsonsu.quant,     # NNL function
                         x   = x,                 # data
                         P   = 0.99,              # percent coverage
                         control = list(parscale = quant.parms), # not clear this was needed
                         hessian = TRUE
                         )
estimates      <- as.numeric(out.optim.quant$par)
standard.error <- as.numeric( sqrt(diag(solve(out.optim.quant$hessian))) )
coef           <- data.frame(estimates, standard.error)
rownames(coef) <- names(out.optim.quant$par)

jparms.optim.quant <- as.list(out.optim.quant$par)
jparms.optim.quant$type <- 'SU'


## check the optimzied PDF from the new parameters
plotspace(1,1)
## first plot histogram with R package fit
out <- hist_nwj(x, type='j', jfit=jparms)
unlist(jparms[1:4])
## add PDF from original optimization
curve(SuppDists::dJohnson(x, jparms.optim), min(x), max(x, 70), col='blue', lty=2, add=TRUE)
unlist(jparms.optim[1:4])
## now add PDF from optimization when reparameterized
out <- nll.johnsonsu.quant(x, P=0.99, jparms=jparms.optim.quant[1:4], debug = TRUE)
unlist(jparms.optim.quant[1:4])
df <- data.frame(x=x, pdf=out$pdf)
df <- df[sort(df$x),]
points(df$x, df$pdf, col='red')
## now add same PDF but using the calculated gamma sicne very different
xord <- sort(x)
jparms.debug <- c(gamma=out$gamma, delta=out$delta, xi=out$xi, lambda=out$lambda)
out.debug <- nll.johnsonsu(xord, jparms=jparms.debug, debug=TRUE)
points(xord, out.debug$pdf, col='red', pch=2)

jparms.debug.su <- as.list(jparms.debug)
jparms.debug.su$type <- 'SU'
out <- hist_nwj(x, type='j', jfit=jparms.debug.su)

## compare parameters
allfits <- rbind(as.numeric(jparms.ExtDist[1:4]),
                 as.numeric(jparms.SuppDists[1:4]),
                 as.numeric(jparms.mle2[1:4]),
                 as.numeric(jparms.optim[1:4]),
                 as.numeric(jparms.debug))
rownames(allfits) <- c('ExtDist::JohnsonSU',
                       'SuppDists::eJohnson',
                       'bbmle::mle2',
                       'optim',
                       'optim quant fit')
allfits

##-------------------------------------------------------------------
set.seed(123);
x=rnorm(1000, mean = 2, sd = 5)

nll.calc <- function(param,data){
  mu=param[1]
  sdev=param[2]
  loglik=dnorm(data, mean = mu, sd = sdev, log = TRUE)
  # cat(mu,sdev,loglik,"\n")
  return(-sum(loglik))
}

theta.start = c(2,4)
ans = optim(par=theta.start, 
            fn=nll.calc, 
            data=x,
            control=list(trace=TRUE),
            method="BFGS")
ans$par
mean(x)
sd(x)
