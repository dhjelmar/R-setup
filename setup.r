### INSTRUCTIONS
##
##
## Run this script with
##     On Windows
##          source('~/Engineering/GitHub_repos/R-setup/setup.r')
##     On Chromebook
##          source('~/GitHub_repos/R-setup/setup.r')
##     On Colab
##          Not sure yet. Likely need to run under Python to access Google Drive.
##          GitHub repo with setup.r file is here:
##             /content/gdrive/MyDrive/Colab Notebooks/github_dhjelmar/R-setup/modules
##
## What it does
##    
##     Load libraries, but first need to install packages (more details below)
##         In terminal, install packages with following if using conda
##             conda install r-essentials
##             conda install r-plotly
##        
##     Source local .r and .R scripts
##        Some of these are R packages that wold not install otherwise
##        Some of these are R functions I wrote

##----------------------------------------------------------------------------

## identify operating system
os <- .Platform$OS.type

### INSTALL PACKAGES
##  Do not use Anaconda for R packages. Packages are often not available and it complicates things.
##  For info, it is faster to install many packages in R than rstudio.

##  If using a project:
##        Create .Rprofile file in top level folder of project containing single line: .libPaths("./libs")
##        Create libs folder
##        Install packages to libs folder

install <- FALSE    # instead of changeing this to TRUE, it may be better to Ctrl-enter each install.packages command
if (install) {
    ## ## install.packages("essentials")   # not sure if this is needed
    install.packages("ggplot2")
    # plotly failed in LMDE installation. Need to first install: sudo apt install libcurl4-openssl-dev
    install.packages("plotly")      
    install.packages("RColorBrewer")
    install.packages("dplyr") 
    install.packages("reticulate")
    install.packages("matlib")   # did not install on Debian 12k
    install.packages("rgl")
    install.packages("IAPWS95")
    install.packages("DT")
    install.packages('tolerance')
    install.packages('DescTools')       # need for Closest() used in plotfit; requires R >= 4.0.0
    install.packages('superml')
    install.packages('xts')
    install.packages('quantmod')
    install.packages('pracma')          # need for pracma::newtonRaphson()
    # ExtDist failed. Need to insatll cmake first: sudo apt install cmake
    install.packages('ExtDist')         # need for Johnson SU distribution
    install.packages('SuppDists')       # need for Johnson SU distribution
    install.packages('maxLik')          # need for MLE (Maximum Likelihood Estimate) fits
    install.packages('expandFunctions') # need for warning.reset()
    install.packages('purrr')           # need for pmap
    install.packages('rpy2')
    ## install.packages('gtools')          # need for mixedorder() sort function
    if (os == 'windows') {
        install.packages("installr")
    }
}

## ## qualityTools no longer supported; had a good qqplot() function
## ## can get an old version if needed
## remotes::install_version('qualityTools', '1.55')  # google search shows v1.55 is latest
## library(qualityTools)  # need for qqplots <-- REMOVED FROM CRAN BECAUSE NOT BEING MAINTAINED


##----------------------------------------------------------------------------

#### LOAD LIBRARIES
## load packages
## following packages do not load on conda
##    tolerance
##    qualityTools
##    plotly
##    IAPWS95
##    rgl

#library(matlib)       # did not install on Debian 12
library(tibble)
library(readxl)
library(ggplot2)       # alternative to base plots (used by plotly)
library(plotly)        # interactive plots
library(RColorBrewer)  # need if use brewer.pal color pallettes
library(dplyr)         # need for left_join(df1,df2,by="byvar")
##                       alternate to be more directly like match/index in Excel would be to use:
##                       parent_dataset[sub_dataset$ID,'column_name_of_value_you_want']
library(reticulate)    # interface with Python (see https://datascienceplus.com/how-to-make-seaborn-pairplot-and-heatmap-in-r-write-python-in-r/)
library(superml)
## library(FrF2)          # fractional factorial design of experiments; will not install on chrome
library(DT)
library(stringi)       # need for stri_split_fixed function
library(stringr)       # need for str_extract function
library(tolerance)
library(purrr)

## ## fonts
## install.packages('extrafont')       # need to use font families
## library(extrafont)
## font_import()                       # only need to run once then comment out
## fonts()                             # use to list available fonts
## usage: loadfonts(device='win', quiet=TRUE)

## print version of R
R.Version()$version.string

loaded <- function() {
    ## identify loaded packages
    package <- data.frame(Package=(.packages()))
    ## (.packages())          # shows packages that are loaded
    ## search()               # little different from above but not sure how (includes more)
    ## sessionInfo()          # info on R-version loaded packages
    ## identify all installed packages and version numbers
    ip <- as.data.frame(installed.packages()[,c(1,3)])
    ## merge to only list loaded packages and version numbers
    package <- merge(package, ip, by = 'Package')   # all=TRUE would keep all lines in both dataframes
    ## convert to matrix to change from factors to character then back to dataframe
    package <- as.data.frame(as.matrix(package))
    return(package)
}
## loaded()

###----------------------------------------------------------------------------
### Load local .R and .r files
##     Source R packages that would not install with conda
##         Instead, source for each of following was copied to ~/ProgramFiles/R_packages and using source to pull in .r and .R files
##         library(qualityTools)  # rnorm, dnorm, rweibull, dweibull, qqPlot (different from qqplot)
##         library(rgl)           # plot3d       
##         library(tolerance)     # need for fitdistr
##         library(IAPWS95)       # water properties: hfT(T)=saturated liq enthalpy,h, in kJ/kg for temperature,T, in K
##                                  hgT(T)=saturated steam enthalpy
##                                  hTp(T,P) = enthalpy for given T and pressure, P
##                                  TSatp(P)
##                                  Tph(P,h) = Temperature for given P and h
##         Note: Only qualityTools would install with intall.packages("qualityTools") but kept having to reinstall because library kept failing
##
##     Source my own .r files in ~/RSTUDIO/modules

## source all files in specified folders herein
if (os == "windows") {
    setup.path <- c("~/Programs/GitHub_home/R-setup/modules")
} else if (os == 'unix') {
    setup.path <- c(##"~/ProgramFiles/R_packages/tolerance/R", # now available on ChromeOS and Debian
        ##"~/ProgramFiles/R_packages/rgl/R",       # available on Debian
        "~/Documents/GitHub/R-setup/modules")          # these are my modules
} else {
    ## assume Colab (.Platform returns NULL)
    setup.path <- c("/content/gdrive/MyDrive/Colab Notebooks/github_dhjelmar/R-setup/modules")
}
r_files <- list.files(setup.path, pattern="*.[rR]$", full.names=TRUE)
for (f in r_files) {
    ## cat("f =",f,"\n")
    source(f)
}

#################################################
## TO UPDATE R AND RSTUDIO
update_r = FALSE
if ((os == 'windows') & (update_r)) {
    ## following is only for Windows and should update R, Rstudio, and installed packages
    installr::updateR()
    ## for Linux, need to do this manually
    ## - install new version of R
    ## - install new version of Rstudio
    ## - find location of new base R packages
    ## - copy non-base R packages from old R installation to new version installation
    ## - run: update.packages()
    ## - check that you got what you think you did
    ##   - version
    ##   - packageStatus()
}
