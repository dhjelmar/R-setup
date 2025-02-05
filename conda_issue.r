## If launch from Conda, may not be able to read some *.r files.
## To fix, set the encoding by default to UTF-8 as follows:
##  
## 1. click on Tools (2nd option starting by the left on the top menu).
## 2. choose Global Options.
## 3. An option box will appear. There choose "code" on the left menu.
## 4. Under code, choose Savings and change the Default Text Encoding to UTF-8.
##
## Now the files that were open as blank should appear with code.


source("F:\\Documents\\01_Dave's Stuff\\Programs\\GitHub_home\\R-setup\\setup.r")
path <- "F:\\Documents\\01_Dave's Stuff\\Programs\\GitHub_home\\R-setup"
setwd(path)

x     <- iris$Sepal.Width
P     <- 0.99  # proportion or coverage
alpha <- 0.01
sided <- 1

plotspace(2,2)
out.fit <- mle.johnsonsu(x, plots=TRUE)
out.tol <- mle.johnsonsu.tol(x, plots=TRUE)

plotspace(1,2)
out.hist <- hist_nwj(x, type='j')
out.qq   <- qqplot_nwj(x, type='j')
