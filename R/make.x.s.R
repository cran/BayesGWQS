#' Forms component group ID vector of X
#'
#' This function returns a vector which lets WQS.fit know the size and order of groups in X
#'
#' @param df A dataframe containing named component variables
#' @param num.groups An integer representing the number of component groups desired
#' @param groups A list, each item in the list being a string vector of variable names for one component group
#'
#' @return A vector of integers, each integer relating how many columns are in each group
#'
#' @examples
#' data("simdata")
#' group_list <- list(c("pcb_118", "pcb_138", "pcb_153", "pcb_180", "pcb_192"),
#'                    c("as", "cu", "pb", "sn"),
#'                    c("carbaryl", "propoxur", "methoxychlor", "diazinon", "chlorpyrifos"))
#' x.s <- make.x.s(simdata, 3, group_list)
#' x.s
#' @export
make.x.s <- function(df, num.groups, groups){

    x.s <- numeric()
    for (i in 1:num.groups){
        x.s[i] <- length(groups[[i]])
    }
    return(x.s)
}


