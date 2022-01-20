#' Bayesian Grouped WQS Regression
#'
#' This function fits a Bayesian grouped weighted quantile sum (BGWQS) regression model.
#'
#' @param y A vector containing outcomes.
#' @param x A matrix of component data.
#' @param z A vector or matrix of controlling covariates.
#' @param x.s A vector of the number of components in each index.
#' @param n.quantiles The number of quantiles to apply to the component data.
#' @param working.dir A file path to the directory.
#' @param n.chains The number of Markov chains; must be a positive integer.
#' @param n.iter The number of total iterations per chain, including burn in.
#' @param n.burnin The number of iterations to discard at the beginning.
#' @param n.thin The thinning rate; must be a positive integer.
#' @param n.adapt The number of adaption iterations.
#' @param DIC Logical; whether or not the user desires the function to return DIC.
#'
#' @return A list which includes BUGS output, sample chains post-burnin, and convergence test results.
#'
#' @examples
#' \dontrun{
#' data("simdata")
#' group_list <- list(c("pcb_118", "pcb_138", "pcb_153", "pcb_180", "pcb_192"),
#'                    c("as", "cu", "pb", "sn"),
#'                    c("carbaryl", "propoxur", "methoxychlor", "diazinon", "chlorpyrifos"))
#' x.s <- make.x.s(simdata, 3, group_list)
#' X <- make.X(simdata, 3, group_list)
#' Y <- simdata$Y
#' work_dir <- tempdir()
#' results <- bgwqs.fit(y = Y, x = X, x.s = x.s, n.quantiles=4,
#'                      working.dir = work_dir,
#'                      n.chains = 1, n.iter = 10000, n.burnin = 5000, n.thin = 1, n.adapt = 500)
#'
#'}
#'
#' @importFrom stats update
#'
#' @export
bgwqs.fit <- function(y, x, z, x.s, n.quantiles=4, working.dir, n.chains = 1,
                      n.iter = 10000, n.burnin = 5000, n.thin = 1, n.adapt = 500, DIC = FALSE){


    ### Error Checks ###
    if(n.iter <= n.burnin)stop("n.iter must be greater than n.burnin")
    if(n.iter <= n.thin)stop("n.iter must be greater than n.thin")
    if(missing(working.dir))stop("Working directory not specified")
    # if(Sys.which("JAGS")=="")stop("JAGS must be installed for this function to work.")
    if(DIC==TRUE && n.chains==1)stop("DIC calculation requires two or more parallel chains.")
    ####################

    orig_user_wd <- getwd() # saving user's original working dir
    on.exit(setwd(orig_user_wd)) # reset working dir back to original on function exit

    ###########################
    #### No Covariates(z) #####
    ###########################

    if (missing(z)) {

        setwd(working.dir)
        K <- length(x.s)
        C <- dim(x)[2] # number of total components


        # Create Model Variables and Setting Initial Values

        N <- length(y)

        C_list <- as.character()

        x_list <- as.character()

        delta_inits <- as.character()


        sigmad_inits <- as.character()


        beta_inits <- as.character()
        beta_inits[1] <- '"beta0"=stats::rnorm(1), '

        sigma_list <- as.character()
        sigma_list[1] <- "sigma0"

        sigma_inits <- as.character()
        sigma_inits[1] <- '"sigma0"=stats::runif(1), '

        beta_list <- as.character()
        beta_list[1] <- "beta0"

        w_vec <- as.character()

        # Initialize model file Objects
        index_vec <- as.character()
        colnum_vec <- as.character()

        counter <- 1

        for (i in 1:K){
            # Model Variables
            C_list[i] <- paste0("C", i)
            assign(C_list[i], x.s[i])

            w_vec[i] <- paste0("w", i)
            endcol <- counter+x.s[i]-1

            if (anyNA(x) == FALSE){

                x_list[i] <- paste0("x", i)
                temp_df <- quantile.fn(x[,counter:endcol], n.quantiles)
                assign(x_list[i], temp_df)
            }

            # Initial Values
            beta_list[i+1] <- paste0("beta", i)
            sigma_list[i+1] <- paste0("sigma", i)

            delta_inits[i] <- paste0('"delta', i, '"=stats::runif(', x.s[i], '), ')
            sigmad_inits[i] <- paste0('"sigmad', i, '"=stats::runif(1), ')
            beta_inits[i+1] <- paste0('"beta', i, '"=stats::rnorm(1), ')
            sigma_inits[i+1] <- paste0('"sigma', i, '"=stats::runif(1), ')

            # model file objects
            index_vec[i] <- paste0("index", i, "[i]")
            colnum_vec[i] <- paste0(counter, ":", endcol)

            counter <- x.s[i]+counter
        }
        sigma_inits[K+1] <- paste0('"sigma', K, '"=stats::runif(1)')

        sig_params <- sigma_list[-1]
        parameters <- c(beta_list, sig_params, w_vec)



        all_inits <- c('list(', delta_inits, sigmad_inits, beta_inits, sigma_inits, ')')
        init_string <- paste(all_inits, sep="", collapse="")

        if (n.chains == 1){
            init1 <- eval(parse(text=init_string))
        } else {
            init_chain_list <- list()
            for (i in 1:n.chains){
                init_chain_list[[i]] <- eval(parse(text=init_string))
            }
            init1 <- init_chain_list
        }

        model_path <- model_gen_nomiss_quant_bin(K, beta_list, index_vec) # Generating model file and returning the path to it

        d <- list()
        d_comps <- c("N", C_list, "y", x_list)
        for (i in 1:length(d_comps)){
            d[[i]] <- get(d_comps[i])
        }
        names(d) <- d_comps



        ############################
        #### With Covariates(z) ####
        ############################

    } else {

        setwd(working.dir)
        K <- length(x.s)
        C <- dim(x)[2] # number of total components
        Zcol <- dim(z)[2] # number of covariate columns

        # Create Model Variables and Setting Initial Values
        N <- length(y)

        C_list <- as.character()
        x_list <- as.character()
        delta_inits <- as.character()
        sigmad_inits <- as.character()
        beta_inits <- as.character()
        beta_inits[1] <- '"beta0"=stats::rnorm(1), '
        sigma_list <- as.character()
        sigma_list[1] <- "sigma0"
        sigma_inits <- as.character()
        sigma_inits[1] <- '"sigma0"=stats::runif(1), '
        beta_list <- as.character()
        beta_list[1] <- "beta0"
        w_vec <- as.character()
        # Covariate Parameter Initilizations
        phi_inits <- as.character()
        sigmaz_inits <- as.character()

        # Initialize model file Objects
        index_vec <- as.character()
        colnum_vec <- as.character()

        counter <- 1

        for (i in 1:K){
            # Model Variables
            C_list[i] <- paste0("C", i)
            assign(C_list[i], x.s[i])

            w_vec[i] <- paste0("w", i)
            endcol <- counter+x.s[i]-1

            x_list[i] <- paste0("x", i)
            temp_df <- quantile.fn(x[,counter:endcol], n.quantiles)
            assign(x_list[i], temp_df)

            # Initial Values

            beta_list[i+1] <- paste0("beta", i)
            sigma_list[i+1] <- paste0("sigma", i)

            delta_inits[i] <- paste0('"delta', i, '"=stats::runif(', x.s[i], '), ')
            sigmad_inits[i] <- paste0('"sigmad', i, '"=stats::runif(1), ')
            beta_inits[i+1] <- paste0('"beta', i, '"=stats::rnorm(1), ')
            sigma_inits[i+1] <- paste0('"sigma', i, '"=stats::runif(1), ')

            # model file objects
            index_vec[i] <- paste0("index", i, "[i]")
            colnum_vec[i] <- paste0(counter, ":", endcol)

            counter <- x.s[i]+counter
        }
        sigma_inits[K+1] <- paste0('"sigma', K, '"=stats::runif(1)') # This initial value must be at the end!
        # Covariate Initial Values
        phi_inits <- paste0('"phi"=stats::rnorm(', Zcol, '),')
        sigmaz_inits <- paste0('"sigmaz"=stats::runif(', Zcol, '),')


        sig_params <- sigma_list[-1]
        parameters <- c(beta_list, sig_params, w_vec, "phi")


        all_inits <- c('list(', delta_inits, sigmad_inits, beta_inits, phi_inits, sigmaz_inits, sigma_inits, ')') # note that sigma_inits is placed at the end
        init_string <- paste(all_inits, sep="", collapse="")

        if (n.chains == 1){
            init1 <- eval(parse(text=init_string))
        } else {
            init_chain_list <- list()
            for (i in 1:n.chains){
                init_chain_list[[i]] <- eval(parse(text=init_string))
            }
            init1 <- init_chain_list
        }

        model_path <- model_gen_nomiss_quant_binz(K, beta_list, index_vec) # Generating model file and returning the path to it

        d <- list()
        d_comps <- c("N", C_list, "y", x_list, "z", "Zcol")
        for (i in 1:length(d_comps)){
            d[[i]] <- get(d_comps[i])
        }
        names(d) <- d_comps
    }


    # MCMC and Convergence Stats

    jags_mod <- rjags::jags.model(file = model_path, data = d, inits = init1, n.chains=n.chains, n.adapt = n.adapt)
    update(object = jags_mod, n.burnin)
    Samples <- rjags::coda.samples(jags_mod, variable.names = parameters, n.iter=n.iter, thin = n.thin)


    if (n.chains > 1){

        if (DIC == TRUE){

            dic_list <- list()
            dic_samps <- rjags::dic.samples(jags_mod, 1000, "pD") # Deviance Information Criterion
            dic_string <- utils::capture.output(dic_samps) # returns a 3-element string with Mean Deviance, penalty, and Penalized Deviance
            dic_list[[1]] <- as.numeric(stringr::str_remove(dic_string[3], stringr::coll("Penalized deviance: "))) # Pulls out Penalized Dev, turns to a numeric
            dic_list[[2]] <- as.numeric(stringr::str_remove(dic_string[2], stringr::coll("penalty ")))
            names(dic_list) <- c("DIC", "pD")

            convergence <- coda::gelman.diag(Samples, multivariate = FALSE)

            final_out <- list(Samples, convergence, dic_list)
            names(final_out) <- c("Samples","Convergence","DIC")

        } else {

            convergence <- coda::gelman.diag(Samples, multivariate = FALSE)

            final_out <- list(Samples, convergence)
            names(final_out) <- c("Samples","Convergence")

        }

    } else { # One chain convergence stats

        Convergence <- coda::geweke.diag(Samples[[1]])
        final_out <- list(Samples, Convergence)
        names(final_out) <- c("Samples", "Convergence")
    }
    return(final_out)
}
