#####################################################
# No Missing Data, Binary Outcome, No Covariates
#####################################################

model_gen_nomiss_quant_bin <- function(K, beta_list, index_vec){

# Chunks of Model Text

beta_line <- c("logit(p[i]) <- beta0")
inprod_vec <- as.character()
weights_list <- list()
index_eff_list <- list()
for (j in 1:K){
    term <- paste0(" + ", beta_list[j+1], " * ", index_vec[j])
    beta_line[j+1] <- term
    inprod_vec[j] <- paste0(index_vec[j], " <- inprod(w", j, "[], x", j, "[i,])")
    weights_list[[j]] <- c(paste0("for (c in 1:C", j, "){"),
                           paste0("w", j, "[c] <- edelta", j, "[c] / sum(edelta", j, "[])"),
                           paste0("edelta", j, "[c] <- exp(delta", j, "[c])"),
                           paste0("delta", j, "[c] ~ dnorm(0, taud", j, ")"),
                           "}",
                           paste0("taud", j, " <- 1/sqrt(sigmad", j, ")"),
                           paste0("sigmad", j, " ~ dunif(0,10)"))
    index_eff_list[[j]] <- c(paste0("beta", j, " ~ dnorm(0, tau", j, ")"),
                             paste0("tau", j, " <- 1/pow(sigma", j, ",2)"),
                             paste0("sigma", j, " ~ dunif(0,100)"))
}
beta_line2 <- paste0(beta_line, collapse="")
weights_chunk <- unlist(weights_list)
intercept_chunk <- c("# Intercept",
                     "beta0 ~ dnorm(0, tau0)",
                     "tau0 <- 1/sqrt(sigma0)",
                     "sigma0 ~ dunif(0,100)")
index_eff_chunk <- unlist(index_eff_list)

# Combining Chunks
model_text <- c("model", "{", "for (i in 1:N)", "{", "y[i]  ~ dbern(p[i])", beta_line2,
                inprod_vec, "}", "meanp <- mean(p[])", "# Index weights", weights_chunk, intercept_chunk,
                "# Index Effects", index_eff_chunk, "}")

# Writing model file and returning path to file

conn <- file("bayesWQS_model.txt")
writeLines(model_text, con = conn) # writing the model file
model_path <- "bayesWQS_model.txt"

close(conn)

return(model_path)
}

####################################################################
# No Missing Data, Binary Outcome, With Covariates
####################################################################

model_gen_nomiss_quant_binz <- function(K, beta_list, index_vec){

    # Chunks of Model Text

    beta_line <- c("logit(p[i]) <- beta0")
    inprod_vec <- as.character()
    weights_list <- list()
    index_eff_list <- list()
    for (j in 1:K){
        term <- paste0(" + ", beta_list[j+1], " * ", index_vec[j])
        beta_line[j+1] <- term
        inprod_vec[j] <- paste0(index_vec[j], " <- inprod(w", j, "[], x", j, "[i,])")
        weights_list[[j]] <- c(paste0("for (c in 1:C", j, "){"),
                               paste0("w", j, "[c] <- edelta", j, "[c] / sum(edelta", j, "[])"),
                               paste0("edelta", j, "[c] <- exp(delta", j, "[c])"),
                               paste0("delta", j, "[c] ~ dnorm(0, taud", j, ")"),
                               "}",
                               paste0("taud", j, " <- 1/sqrt(sigmad", j, ")"),
                               paste0("sigmad", j, " ~ dunif(0,10)"))
        index_eff_list[[j]] <- c(paste0("beta", j, " ~ dnorm(0, tau", j, ")"),
                                 paste0("tau", j, " <- 1/pow(sigma", j, ",2)"),
                                 paste0("sigma", j, " ~ dunif(0,100)"))
    }
    beta_line2 <- paste0(beta_line, collapse="")
    beta_line2 <- paste0(beta_line2, " + phi_z[i]")
    phi_z_inprod <- "phi_z[i] <- inprod(phi[], z[i,])"
    weights_chunk <- unlist(weights_list)
    covariate_chunk <- c("# Covariate Effects",
                         "for (j in 1:Zcol){",
                         "phi[j] ~ dnorm(0, tauz[j])",
                         "tauz[j] <- 1/pow(sigmaz[j],2)",
                         "sigmaz[j] ~ dunif(0,100)",
                         "}")
    intercept_chunk <- c("# Intercept",
                         "beta0 ~ dnorm(0, tau0)",
                         "tau0 <- 1/sqrt(sigma0)",
                         "sigma0 ~ dunif(0,100)")
    index_eff_chunk <- unlist(index_eff_list)

    # Combining Chunks
    model_text <- c("model", "{", "for (i in 1:N)", "{", "y[i]  ~ dbern(p[i])", beta_line2,
                    inprod_vec, phi_z_inprod, "}", "meanp <- mean(p[])", "# Index weights", weights_chunk, covariate_chunk, intercept_chunk,
                    "# Index Effects", index_eff_chunk, "}")

    # Writing model file and returning path to file

    conn <- file("bayesWQS_model.txt")
    writeLines(model_text, con = conn) # writing the model file
    model_path <- "bayesWQS_model.txt"

    close(conn)

    return(model_path)
}




