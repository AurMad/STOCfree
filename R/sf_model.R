#' Draws samples from the posterior distributions of model parameters
#'
#' @param STOCfree_data a STOC free data object
#' @param n_chains number of MCMC chains
#' @param n_burnin number of burnin iterations
#' @param n_iter number of iterations to monitor
#' @param n_thin thinning interval for monitors
#' @param method method to be passed to the runjags::run.jags() function. Default is parallel.
#' @param save_model if TRUE, the JAGS model code is saved to out_path
#' @param save_output if TRUE, the JAGS model output is saved to out_path in a tidy format using the tidybayes package
#' @param out_path folder where model code and output are saved. By default, a STOCfree_files folder is created in the working directory
#' @param ...
#'
#' @export
STOCfree_model <- function(STOCfree_data,
                                n_chains = 4,
                                n_burnin = 5000,
                                n_iter = 5000,
                                n_thin = 20,
                                method = "parallel",
                                save_model = TRUE,
                                save_output = TRUE,
                                out_path = "STOCfree_files",
                                ...){

  ## folder in which the different files are saved
  STOCfree_path <- STOCfree_files(out_path)

  ## running the JAGS model
  JAGS_samples <- runjags::run.jags(
    model = write_JAGS_model(STOCfree_data),
    monitor = JAGS_monitor(STOCfree_data),
    data = STOCfree_JAGS_data(STOCfree_data),
    n.chains = n_chains,
    burnin = n_burnin,
    sample = n_iter,
    thin = n_thin,
    method = method,
    ...)

  ## model file saved as a text file
  if(save_model == TRUE){

    cat(write_JAGS_model(STOCfree_data),
        file = paste0(STOCfree_path ,"/model.txt"))

  }

  ## model results saved in tidy format
  if(save_output == TRUE){

    ## making the output tidy with tidybayes
    tidy_output <- STOCfree_tidy_output(JAGS_samples)

    ## saving parameter values
    write.csv(tidy_output$parameters, file = paste0(STOCfree_path, "/parameters.csv"),
              row.names = FALSE)

    ## saving predicted probabilities of latent status
    write.csv(tidy_output$predictions, file = paste0(STOCfree_path , "/predictions.csv"),
              row.names = FALSE)

  }

  return(JAGS_samples)

  }


## Function to create the list of variables to monitor in JAGS
JAGS_monitor <- function(STOCfree_data){

  ## variables to monitor in JAGS
  parameters_to_monitor  <- c("Se", "Sp", "tau2")
  predictions_to_monitor <- c("predicted_proba")

  if(attributes(sfd)$`number of risk factors` == 0) parameters_to_monitor <- c(parameters_to_monitor, "tau1")
  if(attributes(sfd)$`number of risk factors` > 0)  parameters_to_monitor <- c(parameters_to_monitor, "theta")

  if(attributes(sfd)$level == "animal") parameters_to_monitor <- c(parameters_to_monitor, "pi_within")


  monitor <- c(parameters_to_monitor, predictions_to_monitor)

  return(monitor)

}

