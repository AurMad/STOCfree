#' Draws samples from the posterior distributions of model parameters
#'
#' @param compiled_model a model compiled by compile_JAGS()
#' @param n_burnin number of burnin iterations
#' @param n_iter number of iterations to monitor
#' @param n_thin thinning interval for monitors
#'
#' @return a list with 2 variables
#' @export
#'
#' @examples
sample_model <- function(compiled_model,
                         n_burnin,
                         n_iter,
                         n_thin){

  UseMethod("sample_model")

}



#' @export
sample_model.default <- function(compiled_model,
                                 n_burnin,
                                 n_iter,
                                 n_thin){

  message("No method defined for this type of data")

}


#' @export
sample_model.herd <- function(compiled_model,
                              n_burnin,
                              n_iter,
                              n_thin){

  class(compiled_model) <- "jags"

  ## iterations for burnin
  message(paste0("Burnin (", n_burnin, " iterations)"))
  update(compiled_model, n.iter = n_burnin)

  ## samples from posterior distributions
  variables_to_save <- c("Se", "Sp", "tau1", "tau2",
                         "predicted_proba", "predicted_status")

  message(paste0("Sampling (", n_iter, " iterations)"))
  samples <- rjags::coda.samples(compiled_model,
                                 variable.names = variables_to_save,
                                 n.iter = n_iter,
                                 thin = n_thin)

  ## Gelman diag from the coda package to assess convergence
  r_hat <- coda::gelman.diag(samples, autoburnin = FALSE)$psrf
  r_hat <- r_hat[-grep("predicted", rownames(r_hat)),]

  ## tidying the results for predcited probabilities
  predictions <- tidybayes::spread_draws(samples,
                                         predicted_proba[herd_id],
                                         predicted_status[herd_id])

  ## tidying the results for model parameters
  n_tests <- length(grep("Se", colnames(samples[[1]])))

  if(n_tests == 1){

    parameters <- tidybayes::spread_draws(samples, Se, Sp, tau1, tau2)

  } else {

    parameters <- tidybayes::spread_draws(samples, Se[..], Sp[..], tau1, tau2)

  }

  list(
    parameters = parameters,
    proba_inf = predictions,
    Gelman_diag = r_hat
  )

}


#' @export
sample_model.herd_rf <- function(compiled_model,
                                 n_burnin,
                                 n_iter,
                                 n_thin){

  class(compiled_model) <- "jags"

  ## iterations for burnin
  message(paste0("Burnin (", n_burnin, " iterations)"))
  update(compiled_model, n.iter = n_burnin)

  ## samples from posterior distributions
  variables_to_save <- c("Se", "Sp", "theta", "tau2",
                         "predicted_proba", "predicted_status")

  message(paste0("Sampling (", n_iter, " iterations)"))
  samples <- rjags::coda.samples(compiled_model,
                                 variable.names = variables_to_save,
                                 n.iter = n_iter,
                                 thin = n_thin)

  ## Gelman diag from the coda package to assess convergence
  r_hat <- coda::gelman.diag(samples, autoburnin = FALSE)$psrf
  r_hat <- r_hat[-grep("predicted", rownames(r_hat)),]

  ## tidying the results for predicted probabilities
  predictions <- tidybayes::spread_draws(samples,
                                         predicted_proba[herd_id],
                                         predicted_status[herd_id])

  ## tidying the results for model parameters
  n_tests <- length(grep("Se", colnames(samples[[1]])))

  if(n_tests == 1){

    parameters <- tidybayes::spread_draws(samples, Se, Sp, tau2, theta[..])

  } else {

    parameters <- tidybayes::spread_draws(samples, Se[..], Sp[..], tau2, theta[..])

  }

  list(
    parameters = parameters,
    proba_inf = predictions,
    Gelman_diag = r_hat

  )

}

#' @export
sample_model.animal <- function(compiled_model,
                                n_burnin,
                                n_iter,
                                n_thin){

  class(compiled_model) <- "jags"

  ## iterations for burnin
  message(paste0("Burnin (", n_burnin, " iterations)"))
  update(compiled_model, n.iter = n_burnin)

  ## samples from posterior distributions
  variables_to_save <- c("Se", "Sp", "pi_within", "tau1", "tau2",
                         "predicted_proba", "predicted_status")

  message(paste0("Sampling (", n_iter, " iterations)"))
  samples <- rjags::coda.samples(compiled_model,
                                 variable.names = variables_to_save,
                                 n.iter = n_iter,
                                 thin = n_thin)

  ## Gelman diag from the coda package to assess convergence
  r_hat <- coda::gelman.diag(samples, autoburnin = FALSE)$psrf
  r_hat <- r_hat[-grep("predicted", rownames(r_hat)),]

  ## tidying the results for predicted probabilities
  predictions <- tidybayes::spread_draws(samples,
                                         predicted_proba[herd_id],
                                         predicted_status[herd_id])

  ## tidying the results for model parameters
  n_tests <- length(grep("Se", colnames(samples[[1]])))

  if(n_tests == 1){

    parameters <- tidybayes::spread_draws(samples, Se, Sp, pi_within, tau1, tau2)

  } else {

    parameters <- tidybayes::spread_draws(samples, Se[..], Sp[..], pi_within, tau1, tau2)

  }

  list(
    parameters = parameters,
    proba_inf = predictions,
    Gelman_diag = r_hat

  )

}


#' @export
sample_model.animal_rf <- function(compiled_model,
                                   n_burnin,
                                   n_iter,
                                   n_thin){

  class(compiled_model) <- "jags"

  ## iterations for burnin
  message(paste0("Burnin (", n_burnin, " iterations)"))
  update(compiled_model, n.iter = n_burnin)

  ## samples from posterior distributions
  variables_to_save <- c("Se", "Sp", "pi_within", "theta", "tau2",
                         "predicted_proba", "predicted_status")

  message(paste0("Sampling (", n_iter, " iterations)"))
  samples <- rjags::coda.samples(compiled_model,
                                 variable.names = variables_to_save,
                                 n.iter = n_iter,
                                 thin = n_thin)

  ## Gelman diag from the coda package to assess convergence
  r_hat <- coda::gelman.diag(samples, autoburnin = FALSE)$psrf
  r_hat <- r_hat[-grep("predicted", rownames(r_hat)),]

  ## tidying the results for predicted probabilities
  predictions <- tidybayes::spread_draws(samples,
                                         predicted_proba[herd_id],
                                         predicted_status[herd_id])

  ## tidying the results for model parameters
  n_tests <- length(grep("Se", colnames(samples[[1]])))

  if(n_tests == 1){

    parameters <- tidybayes::spread_draws(samples, Se, Sp, pi_within, tau2, theta[..])

  } else {

    parameters <- tidybayes::spread_draws(samples, Se[..], Sp[..], pi_within, tau2, theta[..])

  }

  list(
    parameters = parameters,
    proba_inf = predictions,
    Gelman_diag = r_hat

  )

}
