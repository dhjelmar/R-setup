readall <- function(datasource,
                    header=TRUE,      # row where header starts
                    nheader=1,        # number of header rows (all but first will be skipped)
                    datastartrow=2,   # row where data starts
                    skip=0,           # number of lines to skip before reading data
                    row.names=NULL,   # vector, column number, or column header
                    sheet=NULL,       # Excel only: sheet name in quotes; defaults to 1st sheet
                    range=NULL,       # Excel only: cell range in quotes including header; overrides headerrow and datastart row
                    rename=TRUE,      # If true, change special characters and spaces to "_"
                    suppress=FALSE    # default prints output of variable names; set to true to suppress output
                    ) {
    ## e.g.,: readall("mydata.xlsx", sheet="Test Geometry", range="b3:x99")

    ## datasource <- "test.xlsx"

    ## determine if datasource is actually a data table rather than a file
    lastone <- substr(datasource, nchar(datasource), nchar(datasource))
            
    ## find extension of the datasource
    ## find a dot that is followed by anything but a dot until the string ends
    lastdot  <- regexpr("\\.[^\\.]*$", datasource)[1]
    ext <- substr(datasource, lastdot+1, nchar(datasource))

    if (lastone == '\n') {
        df <- read.table(text=datasource, header=TRUE)
        datasourcename <- 'internal table'
    } else if (ext == 'csv') {
        if (nheader == 1) {
            df <- read.table(datasource, header=TRUE, sep=",", skip=skip)
        } else {
            hd <- read.table(datasource, header=TRUE,  sep=",", skip=skip,         nrows=nheader)
            df <- read.table(datasource, header=FALSE, sep=",", skip=skip+nheader)
            names(df) <- colnames(hd)
        }
        datasourcename <- datasource
    } else if (ext == 'xls' | ext == 'xlsx') {
        if (nheader == 1) {
            ## df <- read_excel(datasource, col_names=TRUE,       skip=skip,         row.names=row.names)
            if (skip > 0) {
                df <- read_excel(datasource, sheet=sheet, col_names=TRUE,       skip=skip)
            } else if (!is.null(range)) {
                df <- read_excel(datasource, sheet=sheet, col_names=TRUE,       range=range)
            } else {
                df <- read_excel(datasource, sheet=sheet, col_names=TRUE)
            }                
        } else {
            hd <- read_excel(datasource, sheet=sheet, col_names=TRUE,       skip=skip,         n_max=nheader, range=range)
            df <- read_excel(datasource, sheet=sheet, col_names=FALSE,      skip=skip+nheader               , range=range)
            names(df) <- colnames(hd)
        }
        datasourcename <- datasource
    } else {
        datasourcename <- datasource
        cat('FATAL ERROR: unknown datasource type\n')
        stop()
    }
    
    if (suppress == FALSE) {
        cat('\nFrom datasource:\n')
        print(datasourcename)
        cat('As read header:\n')
        print(names(df))
    }

    ## use regular expression to replace space, parenthesis, /, +, -, ., and comma with _
    if (rename == TRUE) names(df) <- gsub(" |\\(|\\)|\\/|\\+|\\-|\\.|\\,", "_", names(df))

    if (suppress == FALSE) {
        cat('Output dataframe names:\n')
        print(names(df))
        ## print(head(df))
    }
    
    return(df)
}

## following is read in as a string where the last two characters are \n
## use above info to identify this as data table rather than external file
## mydata <- "
## num     x1      x2      x3      x4      x5      x6
## 1 	55 	10 	5 	1580 	17.5 	20 
## 2 	75 	10 	5 	1580 	29 	26 
## 3 	55 	25 	5 	1580 	29 	20 
## "
## df <- read.table(text=data, header=TRUE)
