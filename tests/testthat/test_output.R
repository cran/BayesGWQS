library(testthat)
library(BayesGWQS)

context("bgwqs.fit() tests")

data("simdata")
group_list <- list(c("pcb_118", "pcb_138"),
                   c("as", "cu"),
                   c("carbaryl", "propoxur"))
x.s <- make.x.s(simdata, 3, group_list)
X <- make.X(simdata, 3, group_list)
Y <- simdata$Y
work_dir <- tempdir()
results1 <- bgwqs.fit(y = Y, x = X, x.s = x.s, n.quantiles=4, working.dir = work_dir,
                     n.iter = 2, n.burnin = 0, n.thin = 0)

test_that("bgwqs.fit output 1",{
    expect_that(str(results1), prints_text("List of 3"))
    expect_that(str(results1[[1]]$sims.list), prints_text("List of 11"))
})

