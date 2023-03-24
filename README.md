# R-setup

## FUNCTIONS
addfit                = adds data and a fit to an existing plot
cd                    = alias to setwd
closestdata           = 
df.duplicates         = identifies duplicates in a dataframe
df.init               = initialize a dataframe
droplevels.all        = drops all levels for a factored parameter
estbound              = prints useful regression estimates to a table
fastmerge(df1,df2)    = merges 2 dataframes with potentially different columns	 
ggcorplot	      = better correlation plot than pairs(df)			 
ggplotRegression(fit) = plot fit and supporting data				 
hist_nwj	      = histogram with normal, Weibull, and Johnson distributions
holes_bruteforce      = 
holes_liu_dlh         = 
holes_liu_functions   = functions used in other holes functions
holes_liu             = 
interpolate           = interpolation function
inv                   = inverse of a matrix
is.nothing            = test for any form of empty
john_z                = transforms data set x to z assuming Johnson distribution
ls                    = alias to list.files(directory)
mindistance           = returns the minimum distance between a point and the nearest point in a dataframe
mle.johnsonsu         = MLE (maximum likelihood estimate) of Johnson SU distribution
                        handles known and censored data
mle.johnsonsu.tol     = tolerance limit for Johnson SU distribution (uses MLE)
                        handles known and censored data
mle.normal            = MLE (maximum likelihood estimate) of normal distribution
                        handles known and censored data
mle.normal.tol	      = tolerance limit for normal distribution (uses MLE)      
                        handles known and censored data
mle.weibull           = MLE (maximum likelihood estimate) of Weibull distribution
                        handles known and censored data
mle.weibull.tol	      = tolerance limit for Weibull distribution (uses MLE)      
                        handles known and censored data
moi                   = calculates 1st moment of intertia
newton.raphson        = Newton Raphson optimization routine
newton.raaphson.target= wrapper on pracma::newtonRaphson() to allow non-zero target value
nonparametric.tol     = non-parameteric (a.k.a. distribution free) tolerance limit
operators             = similar to %in% but preserves matching order
pairsdf               = variation on standard pairs function
plotdf                = create plot from dataframe (have not used this in a while; not sure there is any advantage over plotfit)
plotfit               = x-y plot with labeling and fitting options
plot_interactive      = x-y plot with ability to use mouse to identify points
plotres               = plots residuals (THIS IS AN ABANDONED WORK IN PROGRESS)
plotspace(r,c)        = define plots per window by defining rows and columns
predictfit            = combine dataframe with regression output
printdf               = print specific # of rows and variables from dataframe
pwd()                 = unix pwd command instead of getwd
qqplot_nwj            = wrapper on qqplot_nwj_only with planned but not coded functionality for censored data
qqplot_nwj_xonly      = qqplots for normal, Weibull, and Johnson distributions and upper tolerance limits
readall               = reads Excel, csv, and other format input
scaledf               = scale (and center) a dataframe
set_plot_dimensions   = sets plot dimensions (set_plot_dimensions(width,height))
shinyplot             = creates an interactive x-y plot
stopexecution         = stops execution
unscaledf             = unscale (and uncenter) a dataframe
warning.reset         = clears warnings in memory


#############################################################################

## USEFUL FUNCTIONS THAT ARE PART OF R
library(help="package")  = useful for seeing info on a package

names(df)     = lists all variables in dataframe

nrow(df)      = returns number of rows in dataframe

ncol(df)      = returns number of columns in dataframe

df$name       = extracts vector name from dataframe df

v1 <- c(1,2,3)

df <- data.frame(v1,v2)

merge(df1,df2,by='v1')

left_join(df1,df2,by='v1')

read_excel("data.xlsx", col_names=TRUE)

df <- subset(df1, (!is.na(v1) & !is.na(v2)))   # create df from df1 but skip if v1 and v2 missing

pairs(df)     = makes cross-plot of all variables in df

interactive bar chart (example: plot_ly(diamonds, x = ~cut, color = ~clarity, colors = "Accent"))

?command     = internet search for help (or can try, for example, help.search("trimmed mean"))

setwd("C:/temp") = set the working folder (make sure use / or \\ since \ is a special character)

getwd()

list.files(".")
