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
sample_model.herd_1test <- function(compiled_model,
                                    n_burnin,
                                    n_iter,
                                    n_thin){

  class(compiled_model) <- "jags"

  ## iterations for burnin
  message(paste0("Burnin (", n_burnin, " iterations)"))
  update(compiled_model, n.iter = n_iter)

  ## samples from posterior distributions
  variables_to_save <- c("Se", "Sp", "tau1", "tau2",
                         "predicted_proba", "predicted_status")

  message(paste0("Sampling (", n_iter, " iterations)"))
  samples <- rjags::coda.samples(compiled_model,
                                 variable.names = variables_to_save,
                                 n.iter = n_iter,
                                 thin = n_thin)

  ## tidying the results for predcited probabilities
  predictions <- tidybayes::spread_draws(samples,
                                         predicted_proba[herd_id],
                                         predicted_status[herd_id])

  ## tidying the results for model parameters
  parameters <- tidybayes::spread_draws(samples, Se, Sp, tau1, tau2)

  ## Percentiles of predicted probabilities of infection according to predicted status
  choice_cutoff <- data.frame(
    cutoff = seq(0, 1, by = .01),
    pred_Se = rep(NA),
    pred_Sp = rep(NA),
    pred_ppv = rep(NA),
    pred_npv = rep(NA)
  )

  for(i in 1:nrow(choice_cutoff)){

    cutoff <- choice_cutoff$cutoff[i]

    true_pos <- nrow(predictions[predictions$predicted_status == 1 &
                                   predictions$predicted_proba > cutoff,])
    true_neg <- nrow(predictions[predictions$predicted_status == 0 &
                                   predictions$predicted_proba <= cutoff,])
    false_pos <- nrow(predictions[predictions$predicted_status == 0 &
                                    predictions$predicted_proba > cutoff,])
    false_neg <- nrow(predictions[predictions$predicted_status == 1 &
                                    predictions$predicted_proba <= cutoff,])

    choice_cutoff$pred_Se[i] <- true_pos / (true_pos + false_neg)
    choice_cutoff$pred_Sp[i] <- true_neg / (true_neg + false_pos)
    choice_cutoff$pred_ppv[i] <- true_pos / (true_pos + false_pos)
    choice_cutoff$pred_npv[i] <- true_neg / (true_neg + false_neg)

  }

  choice_cutoff$pred_ppv[is.nan(choice_cutoff$pred_npv)] <- choice_cutoff$pred_Se[is.nan(choice_cutoff$pred_npv)]
  choice_cutoff$pred_ppv[is.nan(choice_cutoff$pred_ppv)] <- choice_cutoff$pred_Sp[is.nan(choice_cutoff$pred_ppv)]

  list(
    parameters = parameters,
    proba_inf = predictions,
    choice_cutoff = choice_cutoff
  )

}


#' @export
sample_model.herd_1test_rf <- function(compiled_model,
                                        n_burnin,
                                        n_iter,
                                        n_thin){

  class(compiled_model) <- "jags"

  ## iterations for burnin
  message(paste0("Burnin (", n_burnin, " iterations)"))
  update(compiled_model, n.iter = n_iter)

  ## samples from posterior distributions
  variables_to_save <- c("Se", "Sp", "theta", "tau2",
                         "predicted_proba", "predicted_status")

  message(paste0("Sampling (", n_iter, " iterations)"))
  samples <- rjags::coda.samples(compiled_model,
                                 variable.names = variables_to_save,
                                 n.iter = n_iter,
                                 thin = n_thin)

  ## tidying the results for predcited probabilities
  predictions <- tidybayes::spread_draws(samples,
                                         predicted_proba[herd_id],
                                         predicted_status[herd_id])

  ## tidying the results for model parameters
  parameters <- tidybayes::spread_draws(samples, Se, Sp, tau2, theta[..])

  ## Percentiles of predicted probabilities of infection according to predicted status
  choice_cutoff <- data.frame(
    cutoff = seq(0, 1, by = .01),
    pred_Se = rep(NA),
    pred_Sp = rep(NA),
    pred_ppv = rep(NA),
    pred_npv = rep(NA)
  )

  for(i in 1:nrow(choice_cutoff)){

    cutoff <- choice_cutoff$cutoff[i]

    true_pos <- nrow(predictions[predictions$predicted_status == 1 &
                                   predictions$predicted_proba > cutoff,])
    true_neg <- nrow(predictions[predictions$predicted_status == 0 &
                                   predictions$predicted_proba <= cutoff,])
    false_pos <- nrow(predictions[predictions$predicted_status == 0 &
                                    predictions$predicted_proba > cutoff,])
    false_neg <- nrow(predictions[predictions$predicted_status == 1 &
                                    predictions$predicted_proba <= cutoff,])

    choice_cutoff$pred_Se[i] <- true_pos / (true_pos + false_neg)
    choice_cutoff$pred_Sp[i] <- true_neg / (true_neg + false_pos)
    choice_cutoff$pred_ppv[i] <- true_pos / (true_pos + false_pos)
    choice_cutoff$pred_npv[i] <- true_neg / (true_neg + false_neg)

  }

  choice_cutoff$pred_ppv[is.nan(choice_cutoff$pred_npv)] <- choice_cutoff$pred_Se[is.nan(choice_cutoff$pred_npv)]
  choice_cutoff$pred_ppv[is.nan(choice_cutoff$pred_ppv)] <- choice_cutoff$pred_Sp[is.nan(choice_cutoff$pred_ppv)]

  list(
    parameters = parameters,
    proba_inf = predictions,
    choice_cutoff = choice_cutoff
  )

}


#' @export
sample_model.herd_ntests <- function(compiled_model,
                                     n_burnin,
                                     n_iter,
                                     n_thin){

  class(compiled_model) <- "jags"

  ## iterations for burnin
  message(paste0("Burnin (", n_burnin, " iterations)"))
  update(compiled_model, n.iter = n_iter)

  ## samples from posterior distributions
  variables_to_save <- c("Se", "Sp", "tau1", "tau2",
                         "predicted_proba", "predicted_status")

  message(paste0("Sampling (", n_iter, " iterations)"))
  samples <- rjags::coda.samples(compiled_model,
                                 variable.names = variables_to_save,
                                 n.iter = n_iter,
                                 thin = n_thin)

  ## tidying the results for predcited probabilities
  predictions <- tidybayes::spread_draws(samples,
                                         predicted_proba[herd_id],
                                         predicted_status[herd_id])

  ## tidying the results for model parameters
  parameters <- tidybayes::spread_draws(samples, Se[..], Sp[..], tau1, tau2)

  ## Percentiles of predicted probabilities of infection according to predicted status
  choice_cutoff <- data.frame(
    cutoff = seq(0, 1, by = .01),
    pred_Se = rep(NA),
    pred_Sp = rep(NA),
    pred_ppv = rep(NA),
    pred_npv = rep(NA)
  )

  for(i in 1:nrow(choice_cutoff)){

    cutoff <- choice_cutoff$cutoff[i]

    true_pos <- nrow(predictions[predictions$predicted_status == 1 &
                                   predictions$predicted_proba > cutoff,])
    true_neg <- nrow(predictions[predictions$predicted_status == 0 &
                                   predictions$predicted_proba <= cutoff,])
    false_pos <- nrow(predictions[predictions$predicted_status == 0 &
                                    predictions$predicted_proba > cutoff,])
    false_neg <- nrow(predictions[predictions$predicted_status == 1 &
                                    predictions$predicted_proba <= cutoff,])

    choice_cutoff$pred_Se[i] <- true_pos / (true_pos + false_neg)
    choice_cutoff$pred_Sp[i] <- true_neg / (true_neg + false_pos)
    choice_cutoff$pred_ppv[i] <- true_pos / (true_pos + false_pos)
    choice_cutoff$pred_npv[i] <- true_neg / (true_neg + false_neg)

  }

  choice_cutoff$pred_ppv[is.nan(choice_cutoff$pred_npv)] <- choice_cutoff$pred_Se[is.nan(choice_cutoff$pred_npv)]
  choice_cutoff$pred_ppv[is.nan(choice_cutoff$pred_ppv)] <- choice_cutoff$pred_Sp[is.nan(choice_cutoff$pred_ppv)]

  list(
    parameters = parameters,
    proba_inf = predictions,
    choice_cutoff = choice_cutoff
  )

}


#' @export
sample_model.herd_ntests_rf <- function(compiled_model,
                                        n_burnin,
                                        n_iter,
                                        n_thin){

  class(compiled_model) <- "jags"

  ## iterations for burnin
  message(paste0("Burnin (", n_burnin, " iterations)"))
  update(compiled_model, n.iter = n_iter)

  ## samples from posterior distributions
  variables_to_save <- c("Se", "Sp", "theta", "tau2",
                         "predicted_proba", "predicted_status")

  message(paste0("Sampling (", n_iter, " iterations)"))
  samples <- rjags::coda.samples(compiled_model,
                          variable.names = variables_to_save,
                          n.iter = n_iter,
                          thin = n_thin)

  ## tidying the results for predcited probabilities
  predictions <- tidybayes::spread_draws(samples,
                              predicted_proba[herd_id],
                              predicted_status[herd_id])

  ## tidying the results for model parameters
  parameters <- tidybayes::spread_draws(samples, Se[..], Sp[..], tau2, theta[..])

  ## Percentiles of predicted probabilities of infection according to predicted status
  choice_cutoff <- data.frame(
    cutoff = seq(0, 1, by = .01),
    pred_Se = rep(NA),
    pred_Sp = rep(NA),
    pred_ppv = rep(NA),
    pred_npv = rep(NA)
  )

  for(i in 1:nrow(choice_cutoff)){

    cutoff <- choice_cutoff$cutoff[i]

    true_pos <- nrow(predictions[predictions$predicted_status == 1 &
                                   predictions$predicted_proba > cutoff,])
    true_neg <- nrow(predictions[predictions$predicted_status == 0 &
                                   predictions$predicted_proba <= cutoff,])
    false_pos <- nrow(predictions[predictions$predicted_status == 0 &
                                    predictions$predicted_proba > cutoff,])
    false_neg <- nrow(predictions[predictions$predicted_status == 1 &
                                    predictions$predicted_proba <= cutoff,])

    choice_cutoff$pred_Se[i] <- true_pos / (true_pos + false_neg)
    choice_cutoff$pred_Sp[i] <- true_neg / (true_neg + false_pos)
    choice_cutoff$pred_ppv[i] <- true_pos / (true_pos + false_pos)
    choice_cutoff$pred_npv[i] <- true_neg / (true_neg + false_neg)

  }

  choice_cutoff$pred_ppv[is.nan(choice_cutoff$pred_npv)] <- choice_cutoff$pred_Se[is.nan(choice_cutoff$pred_npv)]
  choice_cutoff$pred_ppv[is.nan(choice_cutoff$pred_ppv)] <- choice_cutoff$pred_Sp[is.nan(choice_cutoff$pred_ppv)]

  list(
    parameters = parameters,
    proba_inf = predictions,
    choice_cutoff = choice_cutoff
  )

}



#' @export
sample_model.animal_1test <- function(compiled_model,
                                       n_burnin,
                                       n_iter,
                                       n_thin){

  class(compiled_model) <- "jags"

  ## iterations for burnin
  message(paste0("Burnin (", n_burnin, " iterations)"))
  update(compiled_model, n.iter = n_iter)

  ## samples from posterior distributions
  variables_to_save <- c("Se", "Sp", "pi_within", "tau1", "tau2",
                         "predicted_proba", "predicted_status")
  message(paste0("Sampling (", n_iter, " iterations)"))
  samples <- rjags::coda.samples(compiled_model,
                                 variable.names = variables_to_save,
                                 n.iter = n_iter,
                                 thin = n_thin)

  ## tidying the results for predcited probabilities
  predictions <- tidybayes::spread_draws(samples,
                                         predicted_proba[herd_id],
                                         predicted_status[herd_id])

  ## tidying the results for model parameters
  parameters <- tidybayes::spread_draws(samples, Se, Sp, pi_within, tau1, tau2)

  ## Percentiles of predicted probabilities of infection according to predicted status
  choice_cutoff <- data.frame(
    cutoff = seq(0, 1, by = .01),
    pred_Se = rep(NA),
    pred_Sp = rep(NA),
    pred_ppv = rep(NA),
    pred_npv = rep(NA)
  )

  for(i in 1:nrow(choice_cutoff)){

    cutoff <- choice_cutoff$cutoff[i]

    true_pos <- nrow(predictions[predictions$predicted_status == 1 &
                                   predictions$predicted_proba > cutoff,])
    true_neg <- nrow(predictions[predictions$predicted_status == 0 &
                                   predictions$predicted_proba <= cutoff,])
    false_pos <- nrow(predictions[predictions$predicted_status == 0 &
                                    predictions$predicted_proba > cutoff,])
    false_neg <- nrow(predictions[predictions$predicted_status == 1 &
                                    predictions$predicted_proba <= cutoff,])

    choice_cutoff$pred_Se[i] <- true_pos / (true_pos + false_neg)
    choice_cutoff$pred_Sp[i] <- true_neg / (true_neg + false_pos)
    choice_cutoff$pred_ppv[i] <- true_pos / (true_pos + false_pos)
    choice_cutoff$pred_npv[i] <- true_neg / (true_neg + false_neg)

  }

  choice_cutoff$pred_ppv[is.nan(choice_cutoff$pred_npv)] <- choice_cutoff$pred_Se[is.nan(choice_cutoff$pred_npv)]
  choice_cutoff$pred_ppv[is.nan(choice_cutoff$pred_ppv)] <- choice_cutoff$pred_Sp[is.nan(choice_cutoff$pred_ppv)]

  list(
    parameters = parameters,
    proba_inf = predictions,
    choice_cutoff = choice_cutoff
  )

  }



#' @export
sample_model.animal_1test_rf <- function(compiled_model,
                                      n_burnin,
                                      n_iter,
                                      n_thin){

  class(compiled_model) <- "jags"

  ## iterations for burnin
  message(paste0("Burnin (", n_burnin, " iterations)"))
  update(compiled_model, n.iter = n_iter)

  ## samples from posterior distributions
  variables_to_save <- c("Se", "Sp", "pi_within", "theta", "tau2",
                         "predicted_proba", "predicted_status")

  message(paste0("Sampling (", n_iter, " iterations)"))
  samples <- rjags::coda.samples(compiled_model,
                                 variable.names = variables_to_save,
                                 n.iter = n_iter,
                                 thin = n_thin)

  ## tidying the results for predcited probabilities
  predictions <- tidybayes::spread_draws(samples,
                                         predicted_proba[herd_id],
                                         predicted_status[herd_id])

  ## tidying the results for model parameters
  parameters <- tidybayes::spread_draws(samples, Se, Sp, pi_within, theta[..], tau2)

  ## Percentiles of predicted probabilities of infection according to predicted status
  choice_cutoff <- data.frame(
    cutoff = seq(0, 1, by = .01),
    pred_Se = rep(NA),
    pred_Sp = rep(NA),
    pred_ppv = rep(NA),
    pred_npv = rep(NA)
  )

  for(i in 1:nrow(choice_cutoff)){

    cutoff <- choice_cutoff$cutoff[i]

    true_pos <- nrow(predictions[predictions$predicted_status == 1 &
                                   predictions$predicted_proba > cutoff,])
    true_neg <- nrow(predictions[predictions$predicted_status == 0 &
                                   predictions$predicted_proba <= cutoff,])
    false_pos <- nrow(predictions[predictions$predicted_status == 0 &
                                    predictions$predicted_proba > cutoff,])
    false_neg <- nrow(predictions[predictions$predicted_status == 1 &
                                    predictions$predicted_proba <= cutoff,])

    choice_cutoff$pred_Se[i] <- true_pos / (true_pos + false_neg)
    choice_cutoff$pred_Sp[i] <- true_neg / (true_neg + false_pos)
    choice_cutoff$pred_ppv[i] <- true_pos / (true_pos + false_pos)
    choice_cutoff$pred_npv[i] <- true_neg / (true_neg + false_neg)

  }

  choice_cutoff$pred_ppv[is.nan(choice_cutoff$pred_npv)] <- choice_cutoff$pred_Se[is.nan(choice_cutoff$pred_npv)]
  choice_cutoff$pred_ppv[is.nan(choice_cutoff$pred_ppv)] <- choice_cutoff$pred_Sp[is.nan(choice_cutoff$pred_ppv)]

  list(
    parameters = parameters,
    proba_inf = predictions,
    choice_cutoff = choice_cutoff
  )

}
