zscore.fn <- function(data){
    data <- as.matrix(data)
    z <- matrix(0, dim(data)[1], dim(data)[2])
    C <- dim(data)[2]
    N <- dim(data)[1]
    for(i in 1:C){
        chem_mean <- mean(data[,i])
        chem_sd <- stats::sd(data[,i])
        for (j in 1:N){
            z[j,i] <- (data[j,i] - chem_mean)/chem_sd
        }
    }
return(z)
}
