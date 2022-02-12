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

x <- rnorm(100)
norm.fit <- function(mu, sigma) {
    -sum(dnorm(x, mu, sigma, log=TRUE))
}
## x is the name of the variable containing the data to be fitted
mle.results <- bbmle::mle2(norm.fit, start=list(mu=1,sigma=1), data=list(x))

plotspace(1,3)
x <- mtcars$mpg
nll.johnson.fit <- function(gamma, delta, xi, lambda) {
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
out <- hist_nwj(x, type='j')
jparms <- out$jparms
mle.results <- bbmle::mle2(nll.johnson.fit, start=jparms[1:4], data=list(x))
mle.results@coef
unlist(out$jparms[1:4])
summary(mle.results)@coef

mle.jparms <- list(gamma  = summary(mle.results)@coef[1],
                   delta  = summary(mle.results)@coef[2],
                   xi     = summary(mle.results)@coef[3],
                   lambda = summary(mle.results)@coef[4],
                   type   = 'SU')
curve(SuppDists::dJohnson(x, mle.jparms), min(x), max(x, 70), col='red', add=TRUE)

## compare qqplots
qqplot_nwj(x, type='j', jfit=    jparms, main='SuppDists auto fit')
qqplot_nwj(x, type='j', jfit=mle.jparms, main='MLE forced JohnsonSU fit')

##------------------------------
## try similar fit with optim
nll.johnson.fit <- function(x, pars) {
  ## x is the name of the variable containing the data to be fitted
  ## x is defined in the bbmle::mle2 function
  gamma <- pars[1]  # optim does not pass the parameters as nicely
  delta <- pars[2]
  xi    <- pars[3]
  lambda <- pars[4]
  # PDF for Johnson SU
  pdf <- delta /( lambda * sqrt(2 * pi)   ) *
    1 / sqrt(1 +            ( (x-xi)/lambda )^2)  *
    exp( -0.5*(gamma + delta * asinh( (x-xi)/lambda ))^2 )
  -sum(log(pdf))
  ## the following is equivalent to the above as long as it picks Johnson SU
  ## type <- 'SU'
  ## -sum(SuppDists::dJohnson(x, parms=c(gamma, delta, xi, lambda, type), log=TRUE))
}
mle.optim = optim(par = jparms[1:4],      # parameters to be optimized and initial guess
            fn  = nll.johnson.fit,  # NNL function
            x   = x,                # data
            control = list(parscale = jparms[1:4]) # not clear this was needed
)
unlist(jparms[1:4])
unlist(mle.jparms[1:4])
mle.optim$par

##------------------------------
## turn around johnson.fit for tolerance limit calculation
