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

if (Sys.which("JAGS")==""){

    test_that("bgwqs.fit no-JAGS error",{
        expect_error(bgwqs.fit(y = Y, x = X, x.s = x.s, n.quantiles=4, working.dir = work_dir, mcmc = "jags",
                               n.iter = 2, n.burnin = 1, n.thin = 1, n.adapt = 100), "JAGS must be installed for this function to work.", fixed = TRUE)
    })
} else {

    expect_error(bgwqs.fit(y = Y, x = X, x.s = x.s, n.quantiles=4, working.dir = work_dir, mcmc = "jags",
                           n.iter = 2, n.burnin = 1, n.thin = 1, n.adapt = 100), NA)

}
