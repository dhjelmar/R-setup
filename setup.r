## comment added to test SSH

### INSTRUCTIONS
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
##        Some of these are R packages that wold not install with conda
##        Some of these are R functions I wrote

#----------------------------------------------------------------------------

## identify operating system
os <- .Platform$OS.type

### INSTALL PACKAGES
##  Two options
##      Install from terminal if using conda (this made installing rstudio easy but all packages not supported)
##      Install inside R using (may be faster in R than rstudio), e.g.:
##             install.packages("r-plotly")

### INSTALL PACKAGES WITH CONDA
## input following into terminal (or Anaconda Prompt)
## conda install r-essentials  # not sure if this is needed
## conda install r-readxl
## conda install r-ggplot2     
## conda install r-plotly      
## conda install r-RColorBrewer
## conda install r-dplyr 
## conda install r-reticulate
## conda install r-matlib      <-- worked on chromebook but windows required conda install -c conda-forge r-matlib
## conda install r-qualityTools     <-- may have needed conda install -c conda-forge r-qualityTools
##
## Following list would not install so will source r functions instead
## conda install r-rgl
## conda install r-tolerance
## conda install r-IAPWS95
##
## Following not tried with Conda but installed inside JL
## install.packages("DT")
## install.packages('tolerance')   # worked on Windows in Rstudio

##----------------------------------------------------------------------------

#### LOAD LIBRARIES
## load packages
## following packages do not load on conda
##    tolerance
##    qualityTools
##    plotly
##    IAPWS95
##    rgl
##
## of the above, only would load inside R
##    install.packages("qualityTools")
##    may also have installed plotly inside R
library(qualityTools)
library(matlib)
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
if (os != 'unix') library(tolerance)

(.packages())          # shows packages that are loaded
##search()              # little different from above but not sure how (includes more)

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
    path <- c("~/Programs/GitHub_home/R-setup/modules")
} else if (os == 'unix') {
    path <- c("~/GitHub_repos/R-setup/modules",
              "~/ProgramFiles/R_packages/tolerance/R",
              "~/ProgramFiles/R_packages/rgl/R"      )
} else {
    # assume Colab (.Platform returns NULL)
    path <- c("/content/gdrive/MyDrive/Colab Notebooks/github_dhjelmar/R-setup/modules")
}
r_files <- list.files(path, pattern="*.[rR]$", full.names=TRUE)
for (f in r_files) {
  ## cat("f =",f,"\n")
  source(f)
}
