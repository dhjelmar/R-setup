# R-setup

# MY FUNCTIONS PLUS A FEW NOTES

## MY FUNCTIONS
cd(directory) = unix cd  command instead of setwd

estbound      = prints useful regression estimates into a table

fastmerge(d1,d2) = merges 2 dataframes with potentially different columns

ggcorplot(df) = better correlation plot than pairs(df)

ggplotRegression(fit) = plot fit and supporting data

histplus  = plot histogram with normal mean and limits <-- delete this?

histplus2 = plot histogram with mean and limits for normal or Weibull

hist_nw   = histogram with normal and Weibull distributions

hist_nwj  = histogram with normal, Weibull, and Johnson distributions

inv       = inverse of a matrix

last(vector)  = returns last entry in vector (short for x(length(x))

ls(directory) = unix ls  command instead of list.files(directory)

pairsdf(df,var) = for dataframe df, plot all variables in vector var against eachother

plot_interactive(df,'x','y') = interactive plot and table

plotfit       = creates scatter plot with linear regression and 95% confidence bounds

plotfitcol    = same as plotfit    but colors points based on a continuous variable

plotfitcold   = same as plotfitcol but colors points based on a discrete variable

plotspace(r,c)= define plots per window by defining rows and columns

predictfit    = combine dataframe with regression output

printdf       = print specific # of rows and variables from dataframe

pwd()         = unix pwd command instead of getwd

qqplot_dlh    = attempt to get qqplot for Johnson with confidence limits (not currently successful)

qqplot_nwj    = side by side qqplots for normal, Weibull, and Johnson distributions and upper tolerance limits

scaledf       = scale (and center) a dataframe

set_plot_dimensions(width,height) = sets plot dimensions

tolerance     = not needed if can install R tolerance package
                (copied from tolerance package)

unscaledf     = unscale (and uncenter) a dataframe


## WORKS IN PROGRESS
plotres

plotxy


## OTHER USEFUL FUNCTIONS THAT ARE PART OF R
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
