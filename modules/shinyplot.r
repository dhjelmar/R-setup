shinyplot <- function(df, xx, yy, xline=NULL, yline=NULL) {
    ## use mouse to drag rectangle over points to see table of the points
    ## need to close the plot window to give control back to command line

    ## extract names of variables for plot
    xlabel <- xx
    ylabel <- yy

    ## add xx and yy columns to dataframe
    xcol <- match(xx, names(df))
    ycol <- match(yy, names(df))
    df$xx   <- df[, xcol]
    df$yy   <- df[, ycol]

    ## setup shiny user interface
    ui <- shiny::fluidPage(
        shiny::plotOutput("plot", brush = "plot_brush"),
        shiny::tableOutput("data")
    )

    if (is.null(xline) | is.null(yline)) {
        ## define shiny plot function
        server <- function(input, output, session) {
            output$plot <- shiny::renderPlot({
                ggplot(df, aes(xx, yy)) + geom_point() + xlab(xlabel) + ylab(ylabel)
                ## plot(df$xx, df$yy, xlab=xlabel, ylab=ylabel)      # does not work
                ## with(df, plot(xx, yy, xlab=xlabel, ylab=ylabel))  # does not work
            }, res = 96)  
            output$data <- renderTable({
                ##        shiny::brushedPoints(df, input$plot_brush)
                result <- shiny::brushedPoints(df, input$plot_brush)
                ## identify the numeric columns of result
                result.num <- as.numeric( which(unlist(lapply(df, is.numeric))) )
                ## format numbers with 4 significant figures
                for (i in result.num) {
                    result[[i]] <- format( signif(result[[i]], 4) )
                }
                result
            })
        }
        
    } else {
        ## add line to shiny plot

        ## first create 2nd dataframe with line info
        linedf <- data.frame(xx=xline, yy=yline)
        
        server <- function(input, output, session) {
            output$plot <- shiny::renderPlot({
                ggplot(NULL, aes(xx, yy)) +  xlab(xlabel) + ylab(ylabel) +
                    ## add data points for scatter plot
                    geom_point(data = df) +
                    ## add separately specified line with points
                    geom_line(data = linedf) + geom_point(data = linedf, shape = 3)  
            }, res = 96)  
            output$data <- renderTable({
                ##        shiny::brushedPoints(df, input$plot_brush)
                result <- shiny::brushedPoints(df, input$plot_brush)
                ## identify the numeric columns of result
                result.num <- as.numeric( which(unlist(lapply(df, is.numeric))) )
                ## format numbers with 4 significant figures
                for (i in result.num) {
                    result[[i]] <- format( signif(result[[i]], 4) )
                }
                result
            })
        }
    }        

    ## call shiny plotting app
    shiny::shinyApp(ui = ui, server = server)
}
## shinyplot(mtcars, 'wt', 'mpg')
