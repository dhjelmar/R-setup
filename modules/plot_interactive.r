plot_interactive <- function(df,xx,yy) {
    ## interactive plot and table to find identified points by line # in df
    ## Usage: plot_interactive(df, "xvar", "yvar")
    ##
    ## function does not work as well as this example
    ## not sure how to replace aes(displ,hwy) with aes(xx,yy) where xx='displ'and yy='hwy'
    ## ideas here: https://stackoverflow.com/questions/2641653/pass-a-data-frame-column-name-to-a-function
    ##     df <- mpg
    ##     m <- highlight_key(df)
    ##     p <- ggplot(m, aes(displ, hwy)) + geom_point()
    ##     gg <- highlight(ggplotly(p), "plotly_selected")
    ##     crosstalk::bscols(gg, DT::datatable(m))
    ##
    m <- highlight_key(df)
    p <- ggplot(m, aes(df[[xx]],df[[yy]])) + geom_point() +
        xlab(xx) +  ylab(yy)
    gg <- highlight(ggplotly(p), "plotly_selected")
    ## DT package
    crosstalk::bscols(gg, DT::datatable(m))
}
##plot_interactive(mpg,'displ','hwy')

#--------------------------------

