## calculate moi (1st moment of inertia) for piecewise linear function y = f(x)
moi <- function(x, y) {

    moi <- 0
    for (i in 2:length(x)) {
        ## for each segment

        if (x[i] != x[i-1]) {
            ## not a step change so segment will add to moi for segment of length base
            base      <-  x[i] - x[i-1]

            ## calculate moi for rectangular part
            rheight   <- min(y[i], y[i-1])
            rarea     <- base * rheight
            rcentroid <- (x[i] + x[i-1]) / 2
            moi      <- moi + rarea * rcentroid

            ## calculate additional moi for triangular part, if any
            if (y[i-1] != y[i]) {
                ## there is additional area to consider
                theight <- max(y[i], y[i-1]) - rheight
                tarea   <- base * theight / 2
                if (y[i-1] < y[i]) {
                    ## triangle height increasing with distance
                    tcentroid <- x[i-1] + base * 2 / 3
                } else {
                    ## triangle height decreases with distance
                    tcentroid <- x[i-1] + base * 1 / 3
                }
                moi <- moi + tarea * tcentroid
            }
        }
    }
    return(moi)
}    

## x <- c(0, 10, 10, 15, 20)
## y <- c(5,  5, 10, 10,  5)
## moi(x,y)
