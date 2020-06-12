quantile.fn <- function(data, n.quantiles){
    q <- matrix(0, dim(data)[1], dim(data)[2])
    I <- dim(data)[2]
    for(i in 1:I){
        a <- rank(data[,i], ties.method = "first")
        q[,i] <- cut(a, stats::quantile(a, probs = c(0:n.quantiles/n.quantiles)), include.lowest=TRUE)}
    q <- q-1
    colnames(q) <- colnames(data)
    return(q)
}
