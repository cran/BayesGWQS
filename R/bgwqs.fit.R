#' Bayesian Grouped WQS Regression
#'
#' This function fits a Bayesian grouped weighted quantile sum (BGWQS) regression model.
#'
#' @param y A vector containing outcomes.
#' @param x A matrix of component data.
#' @param z A vector or matrix of covariates.
#' @param x.s A vector of the number of components in each index.
#' @param n.quantiles The number of quantiles to apply to the component data.
#' @param working.dir A file path to the directory.
#' @param n.iter The number of total iterations per chain, including burn in.
#' @param n.burnin The number of iterations to discard at the beginning.
#' @param n.thin The thinning rate, which must be a positive integer.
#' @param debug False by default, when true OpenBUGS remains open for further investigation.
#'
#' @return A list of 3 which includes the BUGS output, coda output, and convergence test results.
#'
#' @examples
#' \donttest{
#' data("simdata")
#' group_list <- list(c("pcb_118", "pcb_138", "pcb_153", "pcb_180", "pcb_192"),
#'                    c("as", "cu", "pb", "sn"),
#'                    c("carbaryl", "propoxur", "methoxychlor", "diazinon", "chlorpyrifos"))
#' x.s <- make.x.s(simdata, 3, group_list)
#' X <- make.X(simdata, 3, group_list)
#' Y <- simdata$Y
#' work_dir <- tempdir()
#' results <- bgwqs.fit(y = Y, x = X, x.s = x.s, n.quantiles=4, working.dir = work_dir,
#'                     n.iter = 10000, n.burnin = 5000, n.thin = 1)
#'}
#'
#' @export
bgwqs.fit <- function(y, x, z = NULL, x.s, n.quantiles=4, working.dir,
                     n.iter = 10000, n.burnin = 5000, n.thin = 1, debug=FALSE){


    ### Error Checks ###
    if(n.iter <= n.burnin)stop("n.iter must be greater than n.burnin")
    if(n.iter <= n.thin)stop("n.iter must be greater than n.thin")
    if(missing(working.dir))stop("Working directory not specified")
    if(!requireNamespace("R2OpenBUGS", quietly = TRUE)) {
            stop("Package R2OpenBUGS needed for this function to work. Please install it.",
                 call. = FALSE)
        }
    ###############

    orig_user_wd <- getwd() # saving user's original working dir
    on.exit(setwd(orig_user_wd)) # reset working dir back to original on function exit

    setwd(working.dir)
    K <- length(x.s)
    C <- dim(x)[2] # number of total components
    q <- quantile.fn(x, n.quantiles)

    # Create Model Variables and Setting Initial Values

    N <- length(y)

    C_list <- as.character()

    x_list <- as.character()

    delta_list <- as.character()

    sigmad_list <- as.character()

    beta_list <- as.character()
    beta_list[1] <- "beta0"
    assign(beta_list[1], .5)

    sigma_list <- as.character()
    sigma_list[1] <- "sigma0"
    assign(sigma_list[1], 1)

    beta_list <- as.character()
    beta_list[1] <- "beta0"

    w_vec <- as.character()

    # Initialize model file Objects
    index_vec <- as.character()

    counter <- 1

    for (i in 1:K){
        # Model Variables
        C_list[i] <- paste0("C", i)
        assign(C_list[i], x.s[i])

        w_vec[i] <- paste0("w", i)

        x_list[i] <- paste0("x", i)
        endcol <- counter+x.s[i]-1
        temp_df <- quantile.fn(x[,counter:endcol], n.quantiles)
        assign(x_list[i], temp_df)
        counter <- x.s[i]+counter

        # Initial Values
        delta_list[i] <- paste0("delta", i)
        assign(delta_list[i], rep(0, get(C_list[i])))
        sigmad_list[i] <- paste0("sigmad", i)
        assign(sigmad_list[i], 1)
        beta_list[i+1] <- paste0("beta", i)
        assign(beta_list[i+1], .1)
        sigma_list[i+1] <- paste0("sigma", i)
        assign(sigma_list[i+1], 1)

        # model file objects
        index_vec[i] <- paste0("index", i, "[i]")
    }

    d <- list()
    d_comps <- c("N", C_list, "y", x_list)
    for (i in 1:length(d_comps)){
        d[[i]] <- d_comps[i]
    }

    init1 <- list()
    init1_comps <- c(delta_list, sigmad_list, beta_list, sigma_list)
    for (i in 1:length(init1_comps)){
        init1[[i]] <- get(init1_comps[i])
    }
    names(init1) <- init1_comps

    # Model

    sig_params <- sigma_list[-1]
    parameters <- c(beta_list, w_vec, sig_params)
    model_path <- model_gen(K, beta_list, index_vec) # Generating model file and returning the path to it

    # MCMC

    MCMCres <- R2OpenBUGS::bugs(data = d, inits = list(init1), parameters.to.save = parameters, model.file = model_path,
                                n.chains = 1, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                                working.directory = working.dir,
                                DIC = TRUE, debug=debug)

    coda_output <- coda::read.coda("CODAchain1.txt", "CODAindex.txt", quiet = TRUE)
    convergence <- coda::geweke.diag(coda_output)

    final_out <- list(MCMCres, coda_output, convergence)
    names(final_out) <- c("BUGS Results", "Coda Output", "Convergence")



    return(final_out)
}





















