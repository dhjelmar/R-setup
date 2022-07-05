readall <- function(datasource,
                    header.row=1,         # row number of the header; header.row=0 means no header
                    data.start.row=NULL,  # row where data starts; default is header.row + 1
                    data.end.row=NULL,    # Excel only: row where data ends; default is end of file
                    sheet=NULL,           # Excel only: sheet name in quotes; defaults to 1st sheet
                    range=NULL,           # Excel only: cell range in quotes including 1 line header; overrides header.row and data.start.row
                    col.names=NULL,       # Used to replace column names
                    rename=TRUE,          # If true, change special characters and spaces in read names to "_";  not used if col.names is specified
                    na=c('NA', 'na', 'NaN', 'N/A', 'n/a', ''), # strings to read as NA rather than character
                    suppress=FALSE        # default prints output of variable names; set to true to suppress output
                    ) {

    ## e.g.,: readall("mydata.xlsx", sheet="Test Geometry", range="b3:x99")
    ## e.g.,: readall("mydata.xlsx", sheet="Test Geometry", header.row=4, data.start.row=6)

    if (is.null(data.start.row)) data.start.row <- header.row + 1
    
    ## datasource <- "test.xlsx"

    ## determine if datasource is actually a data table rather than a file
    lastone <- substr(datasource, nchar(datasource), nchar(datasource))
            
    ## find extension of the datasource
    ## find a dot that is followed by anything but a dot until the string ends
    lastdot  <- regexpr("\\.[^\\.]*$", datasource)[1]
    ext <- substr(datasource, lastdot+1, nchar(datasource))

    if (lastone == '\n') {
        df <- read.table(text=datasource, header=TRUE, stringsAsFactors=FALSE, na.strings=na)
        datasourcename <- 'internal table'
    } else if (ext == 'csv' | ext == 'txt') {
        if (header.row == 0) {
            ## read data; no header
            df <- read.table(datasource, header=FALSE, sep=",", skip=data.start.row-1, stringsAsFactors=FALSE, na.strings=na)
        } else {
            ## read header
            hd <- read.table(datasource, header=TRUE,  sep=",", skip=header.row-1, stringsAsFactors=FALSE, na.strings=na)
            ## read data
            df <- read.table(datasource, header=FALSE, sep=",", skip=data.start.row-1, stringsAsFactors=FALSE, na.strings=na)
            ## assign header names to columns
            names(df) <- colnames(hd)
        }
        datasourcename <- datasource

    } else if (ext == 'xls' | ext == 'xlsx') {
        if (!is.null(range)) {
            ## range is defined so this takes precedence
            if (header.row !=0) {
                ## simple structure with a single header row followed by data in the provided range
                df <- read_excel(datasource, sheet=sheet, col_names=TRUE,  range=range, na=na)
            } else {
                ## no header provided
                df <- read_excel(datasource, sheet=sheet, col_names=FALSE, range=range, na=na)
            }

        } else {
            if (!is.null(data.end.row)) {
                ## set max number of rows to read
                n_max <- data.end.row - data.start.row + 1
            } else {
                ## use default to read through the end of the file
                n_max <- NULL
            }
            if (header.row == 0) {
                ## no header row supplied
                if (is.null(n_max)) {
                    ## read to end of file
                    df <- read_excel(datasource, sheet=sheet, col_names=FALSE, skip=data.start.row-1, na=na)
                } else {
                    ## read specified number of rows
                    df <- read_excel(datasource, sheet=sheet, col_names=FALSE, skip=data.start.row-1, n_max=n_max, na=na)
                }                    
            } else if (data.start.row == header.row + 1) {
                ## data immediately follows single row header so read both together
                if (is.null(n_max)) {
                    ## read to end of file
                    df <- read_excel(datasource, sheet=sheet, col_names=TRUE, skip=header.row-1, na=na)
                } else {
                    ## add 1 to n_max to read header plus specified number of data rows
                    df <- read_excel(datasource, sheet=sheet, col_names=TRUE, skip=header.row-1, n_max=n_max+1, na=na)
                }
            } else {
                ## lines between header row and data
                ## read header
                hd <- read_excel(datasource, sheet=sheet, col_names=TRUE,  skip=header.row-1, n_max=1, na=na)
                if (is.null(n_max)) {
                    ## read data to end of file
                    df <- read_excel(datasource, sheet=sheet, col_names=FALSE, skip=data.start.row-1, na=na)
                } else {
                    ## read specified number of data rows
                    df <- read_excel(datasource, sheet=sheet, col_names=FALSE, skip=data.start.row-1, n_max=n_max, na=na)
                }                    
                ## assign header names to columns (added [] to following to strip off blank columns with header names)
                names(df) <- colnames(hd)[1:length(names(df))]
            }                
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

    ## use regular expression to replace space, parenthesis, /, +, -, and comma with _
    if (rename == TRUE) names(df) <- gsub(" |\\(|\\)|\\/|\\+|\\-|\\,", "_", names(df))

    if (!is.null(col.names)) names(df) <- col.names
    
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
## df <- read.table(text=mydata, header=TRUE)
## df <- readall(mydata, suppress=TRUE)
## df <- readall('test.xlsx', header.row=2, data.end.row=5)
## df <- readall('test.xlsx', range='b2:d11')
## readall('readall.xlsx', range='b2:e5', col.names=c('fred', 'Ethyl Mertz', 'lucy', 'desi'))
