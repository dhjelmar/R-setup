source('~/Programs/GitHub_home/R-setup/setup.r')

# create dataframe of weibull data
num   <- 1000
shape <- 2
scale <- 1
df <- rweibull(num, shape, scale)
df

# create 2x2 space for 4 plots
plotspace(2,2)

# create histogram with normal, Weibull, and Johnson SU fits and tolerance limits
out_hist <- hist_nwj(df, type='nwj', side='upper', sided=1, P=0.99, conf=0.99)

# create qq plots
out_n <- qqplot_nwj(df, type='n')
out_w <- qqplot_nwj(df, type='w')
out_j <- qqplot_nwj(df, type='j')

# distribution free statistics
out_non <- nonparametric.tol(num, conf=0.99, P=0.99)
out_non$message
out_non$df

highpoints <- 2
out_non <- nonparametric.tol(num, conf=0.99, tol.index=num-highpoints)
out_non$message
out_non$df

################################################################################
# 1-sided vs. 2-sided tolerance bounds

set.seed(100)
x <- rnorm(1000,mean=0,sd=1)

###########################

# using tolerance factors (k)
# tolerance factor table
# https://ntrs.nasa.gov/api/citations/19670023646/downloads/19670023646.pdf
k <- 4.383    # for n=25, conf=0.95, and P=0.99
k <- 2.719    # n=1000, conf=0.99, P (reliabiilty) = 0.99
k <- 2.036    # n=1000, conf=0.95, P (reliabiilty) = 0.95
lower_tol_limit = mean(x) - k * sd(x) 
upper_tol_limit = mean(x) + k * sd(x)
lower_tol_limit                                 # -2.081
upper_tol_limit                                 #  2.115

# my 2-sided MLE function gives a very 2-sided bounds to the above
mle.normal.tol(x, sided=2, conf=0.95, P=0.95)   # -2.115, 2.148
out <- hist_nwj(x, type='n', sided=2, conf=0.95, P=0.95) 
out$tol_out_norm$tolerance                      # -2.115, 2.148

# using my MLE functions for 1-sided with same conf and P not as wide as 2-sided bounds
mle.normal.tol(x, side='lower', sided=1, conf=0.95, P=0.05)  # -1.762  not as low
mle.normal.tol(x, side='upper', sided=1, conf=0.95, P=0.95)  #  1.780  not as high

# what needs to be done to use 1-sided bounds to match 2-sided bounds?
# 2-sided alpha above was 0.05
#     conf= 1-alpha/sided = 1 - 0.05/2 = 0.975
# 2-sided coverage (P) was 95%, so missing 0.025 on the low side and 0.025 on the high side
#     1-sided, lower P = 0.025 
#     1-sided, upper P = 1 - 0.025 = 0.975
mle.normal.tol(x, side='lower',  sided=1, conf=0.975, P=0.025)  # -2.115 <- close to 2-sided lower from NTRS
mle.normal.tol(x, side='upper',  sided=1, conf=0.975, P=0.975)  #  2.148 <- close to 2-sided upper from NTRS
out <- hist_nwj(x, type='n', side='lower', sided=1, conf=0.975, P=0.025) 
out$tol_out_norm$tolerance                                     #  -2.115
out <- hist_nwj(x, type='n', side='upper', sided=1, conf=0.975, P=0.975) 
out$tol_out_norm$tolerance                                     #  2.148

# using standard R functions
conf <- 0.95
tolerance::normtol.int(x, side=2, alpha=  1 - conf        , P=0.95)   # -2.082, 2.115 <- almost identical to NTRS
tolerance::normtol.int(x, side=1, alpha=  1 - conf        , P=0.95)   # -1.763, 1.797 <- not as wide as 2-sided bounds (same as found above)
# what needs to be done to use 1-sided bounds to match 2-sided bounds?
sided <- 2
alpha <- (1-conf)*sided   # 0.1
tolerance::normtol.int(x, side=1, alpha= (1 - conf) * sided, P=0.975) # -2.077, 2.110 <- close to above 2-sided

