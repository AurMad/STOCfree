#' JAGS implementation of the STOC free model - deprecated
#'
#' @param STOCfree_data
#' @param n_chains
#' @param n_burnin
#' @param n_iter
#' @param n_thin
#' @param method
#' @param save_model
#' @param save_data
#' @param save_output
#' @param out_path
#' @param ...
#'
#' @return
#' @export
STOCfree_model <- function(STOCfree_data,
                           n_chains = 4,
                           n_burnin = 5000,
                           n_iter = 5000,
                           n_thin = 20,
                           method = "parallel",
                           save_model = TRUE,
                           save_data = FALSE,
                           save_output = TRUE,
                           out_path = "STOCfree_files",
                           ...){


  message("The STOCfree_model() function will soon be removed. Please use STOCfree_JAGS() instead.")

  STOCfree_JAGS(STOCfree_data = STOCfree_data,
                n_chains = n_chains,
                n_burnin = n_burnin,
                n_iter = n_iter,
                n_thin = n_thin,
                method = method,
                save_model = save_model,
                save_data = save_data,
                save_output = save_output,
                out_path = out_path,
                ...)

  }


#' JAGS implementation of the STOC free model
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
STOCfree_JAGS <- function(STOCfree_data,
                                n_chains = 4,
                                n_burnin = 5000,
                                n_iter = 5000,
                                n_thin = 20,
                                method = "parallel",
                                save_model = TRUE,
                                save_data = FALSE,
                                save_output = TRUE,
                                out_path = "STOCfree_files",
                                ...){

  ## folder in which the different files are saved
  STOCfree_path <- STOCfree_files(out_path)

  ## model file saved as a text file
  if(save_model == TRUE){

    cat(write_JAGS_model(STOCfree_data),
        file = paste0(STOCfree_path ,"/model.txt"))

  }

  ## JAGS dataset created and saved if save_model = TRUE
  STOCfree_JAGS_data <- STOCfree_JAGS_data(STOCfree_data)
  if(save_data == TRUE){

  save(STOCfree_JAGS_data,
      file = paste0(STOCfree_path ,"/data.RData"))

  }

  ## running the JAGS model
  JAGS_samples <- runjags::run.jags(
    model = write_JAGS_model(STOCfree_data),
    monitor = JAGS_monitor(STOCfree_data),
    data = STOCfree_JAGS_data,
    n.chains = n_chains,
    burnin = n_burnin,
    sample = n_iter / n_thin,
    thin = n_thin,
    inits = STOCfree_model_inits(STOCfree_data, n_chains, engine = "JAGS"),
    method = method,
    ...)

  JAGS_samples$herd_id_corresp <- STOCfree_data$herd_id_corresp

  ## model results saved in tidy format
  if(save_output == TRUE){

    ## saving parameter values
    write.csv(extract_STOCfree_param(JAGS_samples), file = paste0(STOCfree_path, "/parameters.csv"),
              row.names = FALSE)

    ## saving predicted probabilities of latent status
    write.csv(extract_STOCfree_pred(JAGS_samples), file = paste0(STOCfree_path , "/predictions.csv"),
              row.names = FALSE)

    ## saving monthly prevalences
    write.csv(extract_STOCfree_month_prev(JAGS_samples, STOCfree_data), file = paste0(STOCfree_path , "/month_prev.csv"),
              row.names = FALSE)

  }

  return(JAGS_samples)

  }


## Function to create the list of variables to monitor in JAGS
JAGS_monitor <- function(STOCfree_data){

  ## variables to monitor in JAGS
  # final list of parameters to monitor
  parameters_to_monitor  <- c("Se", "Sp", "tau2")

  if(attributes(STOCfree_data)$`number of risk factors` == 0) parameters_to_monitor <- c(parameters_to_monitor, "tau1")
  if(attributes(STOCfree_data)$`number of risk factors` > 0)  parameters_to_monitor <- c(parameters_to_monitor, "theta")

  if(attributes(STOCfree_data)$level == "animal") parameters_to_monitor <- c(parameters_to_monitor, "pi_within")


  monitor <- c(parameters_to_monitor, "predicted_proba", "month_prev")

  return(monitor)

}


#' Stan implementation of the STOC free model - deprecated
#'
#' @param STOCfree_data
#' @param n_chains
#' @param n_iter
#' @param n_thin thinning interval for monitors
#' @param save_output
#' @param out_path
#'
#' @return
#' @export
#'
#' @examples
STOCfree_model_Stan <- function(STOCfree_data,
                                n_chains = 4,
                                n_iter = 1000,
                                n_thin = 1,
                                save_model = TRUE,
                                save_data = FALSE,
                                save_output = TRUE,
                                out_path = "STOCfree_files"){


  message("The STOCfree_model_Stan() function will soon be removed. Please use STOCfree_Stan() instead.")

  STOCfree_Stan(STOCfree_data = STOCfree_data,
                n_chains = n_chains,
                n_iter = n_iter,
                n_thin = n_thin,
                save_output = save_output,
                out_path = out_path)

}

#' Stan implementation of the STOC free model
#'
#' @param STOCfree_data a STOC free data object
#' @param save_output if TRUE, the JAGS model output is saved to out_path in a tidy format using the tidybayes package
#' @param out_path folder where model code and output are saved. By default, a STOCfree_files folder is created in the working directory
#' @param n_chains number of MCMC chains
#' @param n_iter number of iterations to monitor
#' @param n_thin thinning interval for monitors
#' @param n_warmup
#' @param save_model if TRUE, the Stan model code is saved to out_path
#' @param save_data
#'
#' @details The code used is an adaptation of Damiano et al. (2017): https://github.com/luisdamiano/stancon18
#'
#' @return
#' @export
STOCfree_Stan <- function(STOCfree_data,
                                n_chains = 4,
                                n_iter = 1000,
                                n_thin = 1,
                                n_warmup = NULL,
                                save_model = TRUE,
                                save_data = FALSE,
                                save_output = TRUE,
                                out_path = "STOCfree_files"){

  ## folder in which the different files are saved
  STOCfree_path <- STOCfree_files(out_path)

  ## creating data object for Stan
  sf_Stan_data <- STOCfree_Stan_data(STOCfree_data)

  if(save_data == TRUE){

    save(sf_Stan_data,
         file = paste0(STOCfree_path ,"/data.RData"))

  }

  ## write the model
  Stan_model <- cmdstanr::write_stan_file(write_Stan_model(STOCfree_data))
  sf_Stan    <- cmdstanr::cmdstan_model(Stan_model)

  ## model file saved as a text file
  if(save_model == TRUE){

    sink(paste0(STOCfree_path ,"/model.txt"))
    sf_Stan$print()
    sink()

  }

  ## sample
  Stan_fit <- sf_Stan$sample(
    data = sf_Stan_data,
    chains = n_chains,
    iter_warmup = n_warmup,
    iter_sampling = n_iter,
    thin = n_thin,
    init = STOCfree_model_inits(STOCfree_data, n_chains, engine = "Stan")
    )

  ## model results saved in tidy format
  if(save_output == TRUE){

    ## folder in which the different files are saved
    STOCfree_path <- STOCfree_files(out_path)

    ## saving parameter values
    write.csv(extract_STOCfree_param(x = Stan_fit), file = paste0(STOCfree_path, "/parameters.csv"),
              row.names = FALSE)

    ## saving predicted probabilities of latent status
    write.csv(extract_STOCfree_pred(x = Stan_fit, STOCfree_data),
              file = paste0(STOCfree_path, "/predictions.csv"),
              row.names = FALSE)

  }

 return(Stan_fit)

}





## this function creates initial values for different model parameters
## using prior distributions
STOCfree_model_inits <- function(STOCfree_data, n_chains, engine){

  n_tests <- attr(STOCfree_data, "number of tests")
  status_dynamics_scale <- attr(STOCfree_data, "status dynamics scale")
  n_risk_factors <- attr(STOCfree_data, "number of risk factors")

  list_inits <- list()

  for(n in 1:n_chains){

    ## initial values for test characteristics
    inits <- list(Se = rep(NA, n_tests), Sp = rep(NA, n_tests))

    for(i in 1:n_tests){

      inits$Se[i] <- rbeta(1, STOCfree_data$test_perf_prior$Se_a[i], STOCfree_data$test_perf_prior$Se_b[i])
      inits$Sp[i] <- rbeta(1, STOCfree_data$test_perf_prior$Sp_a[i], STOCfree_data$test_perf_prior$Sp_b[i])

    }

    ## dynamics on the probability scale
    if(status_dynamics_scale == "proba"){


      if(engine == "JAGS"){

        inits <- c(inits, tau2 = NA)

        inits$tau2 <- rbeta(1, STOCfree_data$inf_dyn_priors["tau2_a"], STOCfree_data$inf_dyn_priors["tau2_b"])

      }


      if(engine == "Stan"){

        inits <- c(inits, pi1 = NA, tau2 = NA)

        inits$pi1  <- rbeta(1, STOCfree_data$inf_dyn_priors["pi1_a"], STOCfree_data$inf_dyn_priors["pi1_b"])
        inits$tau2 <- rbeta(1, STOCfree_data$inf_dyn_priors["tau2_a"], STOCfree_data$inf_dyn_priors["tau2_b"])

      }

      ## no risk factor, dynamics on the probability scale
      if(n_risk_factors == 0){

        inits <- c(inits, tau1 = NA)
        inits$tau1 <- rbeta(1, STOCfree_data$inf_dyn_priors["tau1_a"], STOCfree_data$inf_dyn_priors["tau1_b"])

      } else {

        inits <- c(inits, list(theta = rep(NA, nrow(STOCfree_data$risk_factors))))

        for(i in 1:nrow(STOCfree_data$risk_factors)){

          inits$theta[i] <- rnorm(1, STOCfree_data$risk_factors$mean_prior[i], STOCfree_data$risk_factors$sd_prior[i])

        }

      }
    } ## dynamics on the probability scale


    ## dynamics on the logit scale
    if(status_dynamics_scale == "logit"){

      if(engine == "JAGS"){

        inits <- c(inits, logit_tau2 = NA)
        inits$logit_tau2 <- rnorm(1, STOCfree_data$inf_dyn_priors["logit_tau2_mean"], STOCfree_data$inf_dyn_priors["logit_tau2_sd"])

      }

      if(engine == "Stan"){

        inits <- c(inits, pi1 = NA, tau2 = NA)

        inits$pi1  <- invlogit(rnorm(1, STOCfree_data$inf_dyn_priors["logit_pi1_mean"], STOCfree_data$inf_dyn_priors["logit_pi1_sd"]))
        inits$tau2 <- invlogit(rnorm(1, STOCfree_data$inf_dyn_priors["logit_tau2_mean"], STOCfree_data$inf_dyn_priors["logit_tau2_sd"]))

      }

      ## no risk factor, dynamics on the probability scale
      if(n_risk_factors == 0){

        if(engine == "JAGS"){

          inits <- c(inits, logit_tau1 = NA)
          inits$logit_tau1 <- rnorm(1, STOCfree_data$inf_dyn_priors["logit_tau1_mean"], STOCfree_data$inf_dyn_priors["logit_tau1_sd"])

        }

        if(engine == "Stan"){

          inits <- c(inits, tau1 = NA)
          inits$tau1 <- invlogit(rnorm(1, STOCfree_data$inf_dyn_priors["logit_tau1_mean"], STOCfree_data$inf_dyn_priors["logit_tau1_sd"]))

        }


      } else {

        inits <- c(inits, list(theta = rep(NA, nrow(STOCfree_data$risk_factors))))

        for(i in 1:nrow(STOCfree_data$risk_factors)){

          inits$theta[i] <- rnorm(1, STOCfree_data$risk_factors$mean_prior[i], STOCfree_data$risk_factors$sd_prior[i])

        }

      }
    } ## dynamics on the probability scale

    list_inits <- c(list_inits, list(inits))

  }

  list_inits
}
