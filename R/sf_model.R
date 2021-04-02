#' JAGS implementation of the STOC free model
#'
#' @param STOCfree_data a STOC free data object
#' @param n_chains number of MCMC chains
#' @param n_burnin number of burnin iterations
#' @param n_iter number of iterations to monitor
#' @param n_thin thinning interval for monitors
#' @param method method to be passed to the runjags::run.jags() function. Default is parallel.
#' @param save_data if TRUE, the data are saved to 'out_path' as a .RData file
#' @param save_inits if TRUE initial values are written to 'out_path' as a text file
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
                                save_inits = TRUE,
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

  ## initial values
  inits <- STOCfree_model_inits(STOCfree_data, n_chains, engine = "JAGS")
  if(save_inits == TRUE){

    sink(paste0(STOCfree_path ,"/inits.txt"))
    print(inits)
    sink()

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


#' Stan implementation of the STOC free model
#'
#' @param STOCfree_data a STOC free data object
#' @param n_chains number of MCMC chains
#' @param n_iter number of iterations to monitor
#' @param n_thin thinning interval for monitors
#' @param n_warmup number of warmup iterations
#' @param out_path folder where model code and output are saved. By default, a STOCfree_files folder is created in the working directory
#' @param save_output if TRUE, the JAGS model output is saved to out_path in a tidy format using the tidybayes package
#' @param save_model if TRUE, the Stan model code is saved to 'out_path'
#' @param save_data if TRUE, the data are saved to 'out_path' as a .RData file
#' @param save_inits if TRUE initial values are written to 'out_path' as a text file
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
                                save_inits = TRUE,
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

  ## initial values
  inits <- STOCfree_model_inits(STOCfree_data, n_chains, engine = "Stan")
  if(save_inits == TRUE){

    sink(paste0(STOCfree_path ,"/inits.txt"))
    print(inits)
    sink()

  }

  ## sample
  Stan_fit <- sf_Stan$sample(
    data = sf_Stan_data,
    chains = n_chains,
    iter_warmup = n_warmup,
    iter_sampling = n_iter,
    thin = n_thin,
    init = inits
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
## pi1, the prevalence of status positives on the first tests is estimated
## from the proportion of test positives on the first test,
## the initial values for sensitivity and specificity
## this was made this way because I got errors when pi1 was sampled
## randomly from the priors
## I figured that pi1 was constrained by these 3 parameters
## this seems to have resolved the problem
STOCfree_model_inits <- function(STOCfree_data, n_chains, engine){

  n_tests <- attr(STOCfree_data, "number of tests")
  status_dynamics_scale <- attr(STOCfree_data, "status dynamics scale")
  n_risk_factors <- attr(STOCfree_data, "number of risk factors")

  ## the initial value for pi1 is determined from other values
  ## which test is the most frequent on the first test occasion?
  test_m1 <- STOCfree_data$test_data[STOCfree_data$test_data$month_id == 1, ]
  tab_test_id <- table(test_m1$test_id)

  test1 <- as.integer(names(tab_test_id)[which(tab_test_id == max(tab_test_id))])

  ## proportion of test positives with the most used test on the first test occasion
  p_tpos_m1 <- nrow(test_m1[test_m1$test_id == test1 & test_m1$test_res == 1,]) / nrow(test_m1[test_m1$test_id == test1,])

  ## function to estimate pi1 from proportion of test positives on first test occasion, Se and Sp
  pi1_est <- function(p_tpos, Se, Sp){

    pi1 <- (p_tpos - 1 + Sp) / (Se + Sp - 1)

    return(pi1)

  }

  ## empty list to store initial values
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

        inits$pi1  <- pi1_est(p_tpos_m1, inits$Se[test1], inits$Se[test1])
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

        inits$pi1  <- pi1_est(p_tpos_m1, inits$Se[test1], inits$Se[test1])
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
