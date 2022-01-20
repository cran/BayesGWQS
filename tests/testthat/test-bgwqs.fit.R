library(testthat)
library(BayesGWQS)

context("bgwqs.fit() tests")

data("simdata")
group_list <- list(c("pcb_118", "pcb_138"),
                   c("as"),
                   c("carbaryl", "propoxur", "cu"))
x.s <- make.x.s(simdata, 3, group_list)
X <- make.X(simdata, 3, group_list)
Y <- simdata$Y
work_dir <- tempdir()

# if (Sys.which("JAGS")==""){ # Test for when user doesn't have JAGS installed
#
#     test_that("bgwqs.fit no-JAGS error",{
#         expect_error(bgwqs.fit(y = Y, x = X, x.s = x.s, n.quantiles=4, working.dir = work_dir,
#                                n.iter = 2, n.burnin = 1, n.thin = 1, n.adapt = 100), "JAGS must be installed for this function to work.", fixed = TRUE)
#     })
#
# } else { # Tests for when user has JAGS installed (most of them)

    test_that("(1) binary, quantile, nomiss, chain=1",{
    expect_error(bgwqs.fit(y = Y, x = X, x.s = x.s, n.quantiles=4, working.dir = work_dir, n.chains = 1,
                           n.iter = 2, n.burnin = 1, n.thin = 1, n.adapt = 100), NA)
    })
#}
