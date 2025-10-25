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
# per looking online

# 1-sided
alpha <- 0.1
conf <- 1 - alpha      # 0.90
P <- 0.95
tolerance::normtol.int(x, side=1, alpha=0.1, P=0.95)           # -1.744, 1.778

# 2-sided
# conf = 1 - alpha_1_sided / 2
alpha <- (1 - conf)*2  # 0.2
P <- (1+P)/2           # 0.975
tolerance::normtol.int(x, side=2, alpha=0.2, P=0.975)          # -2.34, 2.37


###########################
# upper, 1-sided
out1U <- hist_nwj(x, type='n', side='upper', sided=1, conf=0.90, P=0.95)     #  1.776
mle.normal.tol(   x,           side='upper', sided=1, conf=0.90, P=0.95)     

# lower, 1-sided
out1L <- hist_nwj(x, type='n', side='lower', sided=1, conf=0.90, P=0.05)     # -1.743
mle.normal.tol(   x,           side='lower', sided=1, conf=0.90, P=0.05)

# 2-sided
# P = 1 - (1-Pupper) - Plower = 0.90 rather since missing 0.05 on each side
# But what to do with confidence (or alpha)?
#
# the following is needed to match the two 1-sided calculations
out2  <- hist_nwj(x, type='n', side='both' , sided=2, conf=0.80, P=0.90)     # -1.743, 1.776
mle.normal.tol(   x,                         sided=2, conf=0.80, P=0.90)
mle.normal.tol(   x,                         sided=2, alpha=0.1*2, P=0.90)
#
# but should it be the following instead since now half the conf is on the low side and half is on the upper side?
out2  <- hist_nwj(x, type='n', side='both' , sided=2, conf=0.95, P=0.90)     # -1.778, 1.812
mle.normal.tol(   x,                         sided=2, alpha=0.1/2, P=0.90)
mle.normal.tol(   x,                         sided=2, conf=0.95, P=0.90)

