#' Generates Plots of weights by group
#'
#' This function takes the object created by the bgwqs.fit function and a vector of group names and generates a
#' random forest variable importance plot for each group. The weights in each group are listed in descending order.
#'
#' @param fit.object The object that is returned by the bgwqs.fit function
#' @param group.names A string vector containing the name of each group included in the BGWQS regression. Will be used for plot titles.
#' @param group.list A list, each item in the list being a string vector of variable names for one component group.
#' @param x.s A vector of the number of components in each index.
#'
#' @return A plot for each group of the BGWQS regression
#'
#'
#' @export
weight.plot <- function(fit.object, group.names, group.list, x.s){
    fit_sum <- summary(fit.object[[1]])[[1]]
    subset_row <- nrow(fit_sum) - sum(x.s) + 1         # row index of first weight row
    fit_sum2 <- fit_sum[c(subset_row:nrow(fit_sum)),]
    rownames(fit_sum2) <- unlist(group.list)
    plot_num <- length(x.s)
    c <- 1 #counter
    for (i in 1:plot_num){
        weight_vec <- fit_sum2[c:(x.s[i]+c-1),1]
        c <- c + x.s[i]
        sorted_vec <- sort(weight_vec)
        graphics::dotchart(sorted_vec, labels = names(sorted_vec), xlab = "Weight", pch = 19)
        graphics::title(main = paste("Weights for", group.names[i]), cex.main = 1.00)
    }
}

