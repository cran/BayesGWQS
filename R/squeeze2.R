squeeze2 <- function(x, bounds = c(min(x), max(x))){

    unders <- length(x[x < bounds[1]])
    overs <- length(x[x > bounds[2]])
    x[x < bounds[1]] <- stats::runif(unders, min = bounds[1], max = bounds[2])
    x[x > bounds[2]] <- stats::runif(overs, min = bounds[1], max = bounds[2])
    return(x)
}
