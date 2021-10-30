shinyplot <- function(df, xx, yy) {
    ## use mouse to drag rectangle over points to see table of the points
    ## need to close the plot window to give control back to command line

    ## extract names of variables for plot
    xlabel <- xx
    ylabel <- yy

    ## add xx and yy columns to dataframe
    xcol <- which(grepl(xx, names(df)))
    ycol <- which(grepl(yy, names(df)))
    df$xx   <- df[, xcol]
    df$yy   <- df[, ycol]

    ## setup shiny user interface
    ui <- fluidPage(
        plotOutput("plot", brush = "plot_brush"),
        tableOutput("data")
    )

    ## define shiny plot function
    server <- function(input, output, session) {
        output$plot <- renderPlot({
            ggplot(df, aes(xx, yy)) + geom_point() + xlab(xlabel) + ylab(ylabel)
            ## plot(df$xx, df$yy, xlab=xlabel, ylab=ylabel)      # does not work
            ## with(df, plot(xx, yy, xlab=xlabel, ylab=ylabel))  # does not work
        }, res = 96)  
        output$data <- renderTable({
            brushedPoints(df, input$plot_brush)
        })
    }

    ## call shiny plotting app
    shinyApp(ui = ui, server = server)
}
## shinyplot(mtcars, 'wt', 'mpg')
