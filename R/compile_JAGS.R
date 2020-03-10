#' Compilation of the JAGS model
#'
#' @param data an object of class STOCfree_data
#' @param n_chains number of MCMC chains
#' @param ...
#'
#' @return
#' a compiled JAGS model
#'
#' @export
#'
#' @examples
compile_JAGS <- function(data,
                         n_chains = 4, ...){

  UseMethod("compile_JAGS")

  }


#' @export
compile_JAGS.default = function(data, n_chains){

  print("No method defined for this type of data")

  }



#' @export
compile_JAGS.herd_1test = function(data, n_chains){

  n_herds <- attr(data, "number of herds")
  month_max <- max(data$test_data$month_id)

  ## test results
  test_res_col <- "test_res"
  test_res <- data$test_data[data$test_data$month_id < month_max, test_res_col]
  ## Indices when test performed in the expanded dataset
  ind_test <- data$test_data$status_id[data$test_data$month_id < month_max]

  ## Herds for which a test was performed on the month to predict
  ind_last_is_test <- which(data$herd_test_data$last_is_test == 1)
  ind_p_test <- data$herd_test_data$ind_p[ind_last_is_test]
  n_pred_test <- length(ind_last_is_test)

  ## Herds for which a test was NOT performed on the month to predict
  ind_last_is_not_test <- which(data$herd_test_data$last_is_test == 0)
  ind_p_no_test <- data$herd_test_data$ind_p[ind_last_is_not_test]
  n_pred_no_test <- length(ind_last_is_not_test)

  ## test results used for prediction
  test_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,test_res_col]), test_res_col]

  ## Priors for test characteristics
  test_char <- data$test_perf_prior

  ## Data used by JAGS
  JAGS_data <- list(
    n_herds = n_herds,
    ind_i = data$herd_test_data$ind_i,
    ind_j = data$herd_test_data$ind_j,
    ind_f = data$herd_test_data$ind_f,
    n_tests_perf = length(ind_test),
    ind_test = as.integer(ind_test),
    test_res = as.integer(test_res),
    test_for_pred = test_for_pred,
    n_pred_test = n_pred_test,
    ind_last_is_test = ind_last_is_test,
    n_pred_no_test = n_pred_no_test,
    ind_p_test = ind_p_test,
    ind_last_is_not_test = ind_last_is_not_test,
    ind_p_no_test = ind_p_no_test,
    pi1_beta_a = data$inf_dyn_priors["pi1_a"],
    pi1_beta_b = data$inf_dyn_priors["pi1_b"],
    tau1_beta_a = data$inf_dyn_priors["tau1_a"],
    tau1_beta_b = data$inf_dyn_priors["tau1_b"],
    tau2_beta_a = data$inf_dyn_priors["tau2_a"],
    tau2_beta_b = data$inf_dyn_priors["tau2_b"],
    Se_beta_a = test_char$Se_a,
    Se_beta_b = test_char$Se_b,
    Sp_beta_a = test_char$Sp_a,
    Sp_beta_b = test_char$Sp_b
  )

  write_JAGS_model(data)

  JAGS_model_compiled <- rjags::jags.model(
    file = "JAGS_model.txt",
    data = JAGS_data,
    n.chains = n_chains)


  if(file.exists("JAGS_model.txt")) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "herd_1test")

  JAGS_model_compiled

}


#' @export
compile_JAGS.herd_1test_rf = function(data, n_chains){

  n_herds <- attr(data, "number of herds")
  month_max <- max(data$test_data$month_id)

  ## test results
  test_res_col <- "test_res"
  test_res <- data$test_data[data$test_data$month_id < month_max, test_res_col]
  ## Indices when test performed in the expanded dataset
  ind_test <- data$test_data$status_id[data$test_data$month_id < month_max]

  ## Herds for which a test was performed on the month to predict
  ind_last_is_test <- which(data$herd_test_data$last_is_test == 1)
  ind_p_test <- data$herd_test_data$ind_p[ind_last_is_test]
  n_pred_test <- length(ind_last_is_test)

  ## Herds for which a test was NOT performed on the month to predict
  ind_last_is_not_test <- which(data$herd_test_data$last_is_test == 0)
  ind_p_no_test <- data$herd_test_data$ind_p[ind_last_is_not_test]
  n_pred_no_test <- length(ind_last_is_not_test)

  ## test results used for prediction
  test_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,test_res_col]), test_res_col]

  ## Priors for risk factors
  risk_factors <- data$risk_factors[data$risk_factors$ref == 0,]
  theta_norm_mean <- risk_factors$mean_prior
  theta_norm_prec <- 1 / risk_factors$sd_prior^2

  ## Priors for test characteristics
  test_char <- data$test_perf_prior

  ## Data used by JAGS
  JAGS_data <- list(
    n_herds = n_herds,
    ind_i = data$herd_test_data$ind_i,
    ind_j = data$herd_test_data$ind_j,
    ind_f = data$herd_test_data$ind_f,
    ind_p = data$herd_test_data$ind_p,
    n_tests_perf = length(ind_test),
    ind_test = as.integer(ind_test),
    test_res = as.integer(test_res),
    test_for_pred = test_for_pred,
    n_pred_test = n_pred_test,
    ind_last_is_test = ind_last_is_test,
    n_pred_no_test = n_pred_no_test,
    ind_p_test = ind_p_test,
    ind_last_is_not_test = ind_last_is_not_test,
    ind_p_no_test = ind_p_no_test,
    pi1_beta_a = data$inf_dyn_priors["pi1_a"],
    pi1_beta_b = data$inf_dyn_priors["pi1_b"],
    tau2_beta_a = data$inf_dyn_priors["tau2_a"],
    tau2_beta_b = data$inf_dyn_priors["tau2_b"],
    Se_beta_a = test_char$Se_a,
    Se_beta_b = test_char$Se_b,
    Sp_beta_a = test_char$Sp_a,
    Sp_beta_b = test_char$Sp_b,
    theta_norm_mean = theta_norm_mean,
    theta_norm_prec = theta_norm_prec,
    n_risk_factors = length(theta_norm_mean),
    risk_factors = as.matrix(data$risk_factor_data[, -(1:3)])
  )

  write_JAGS_model(data)

  JAGS_model_compiled <- rjags::jags.model(
    file = "JAGS_model.txt",
    data = JAGS_data,
    n.chains = n_chains)


  if(file.exists("JAGS_model.txt")) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "herd_1test_rf")

  JAGS_model_compiled

}

#' @export
compile_JAGS.herd_ntests = function(data, n_chains){

  n_herds <- attr(data, "number of herds")
  month_max <- max(data$test_data$month_id)

  ## test results
  test_res_col <- "test_res"
  test_res <- data$test_data[data$test_data$month_id < month_max, test_res_col]
  test_id  <- data$test_data$test_id[data$test_data$month_id < month_max]
  ## Indices when test performed in the expanded dataset
  ind_test <- data$test_data$status_id[data$test_data$month_id < month_max]

  ## Herds for which a test was performed on the month to predict
  ind_last_is_test <- which(data$herd_test_data$last_is_test == 1)
  ind_p_test <- data$herd_test_data$ind_p[ind_last_is_test]
  n_pred_test <- length(ind_last_is_test)

  ## Herds for which a test was NOT performed on the month to predict
  ind_last_is_not_test <- which(data$herd_test_data$last_is_test == 0)
  ind_p_no_test <- data$herd_test_data$ind_p[ind_last_is_not_test]
  n_pred_no_test <- length(ind_last_is_not_test)

  ## test results used for prediction
  test_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,test_res_col]), test_res_col]
  test_id_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,test_res_col]), "test_id"]

  ## Priors for test characteristics
  test_char <- data$test_perf_prior

  ## Data used by JAGS
  JAGS_data <- list(
    n_herds = n_herds,
    ind_i = data$herd_test_data$ind_i,
    ind_j = data$herd_test_data$ind_j,
    ind_f = data$herd_test_data$ind_f,
    n_tests_perf = length(ind_test),
    ind_test = as.integer(ind_test),
    test_res = as.integer(test_res),
    test_id = as.integer(test_id),
    test_for_pred = test_for_pred,
    test_id_for_pred = test_id_for_pred,
    n_pred_test = n_pred_test,
    ind_last_is_test = ind_last_is_test,
    n_pred_no_test = n_pred_no_test,
    ind_p_test = ind_p_test,
    ind_last_is_not_test = ind_last_is_not_test,
    ind_p_no_test = ind_p_no_test,
    pi1_beta_a = data$inf_dyn_priors["pi1_a"],
    pi1_beta_b = data$inf_dyn_priors["pi1_b"],
    tau1_beta_a = data$inf_dyn_priors["tau1_a"],
    tau1_beta_b = data$inf_dyn_priors["tau1_b"],
    tau2_beta_a = data$inf_dyn_priors["tau2_a"],
    tau2_beta_b = data$inf_dyn_priors["tau2_b"],
    n_tests = nrow(test_char),
    Se_beta_a = test_char$Se_a,
    Se_beta_b = test_char$Se_b,
    Sp_beta_a = test_char$Sp_a,
    Sp_beta_b = test_char$Sp_b
  )

  write_JAGS_model(data)

  JAGS_model_compiled <- rjags::jags.model(
    file = "JAGS_model.txt",
    data = JAGS_data,
    n.chains = n_chains)


  if (file.exists("JAGS_model.txt")) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "herd_ntests")

  JAGS_model_compiled

}

#' @export
compile_JAGS.herd_ntests_rf = function(data, n_chains){

  n_herds <- attr(data, "number of herds")
  month_max <- max(data$test_data$month_id)

  ## test results
  test_res_col <- "test_res"
  test_res <- data$test_data[data$test_data$month_id < month_max, test_res_col]
  test_id  <- data$test_data$test_id[data$test_data$month_id < month_max]
  ## Indices when test performed in the expanded dataset
  ind_test <- data$test_data$status_id[data$test_data$month_id < month_max]

  ## Herds for which a test was performed on the month to predict
  ind_last_is_test <- which(data$herd_test_data$last_is_test == 1)
  ind_p_test <- data$herd_test_data$ind_p[ind_last_is_test]
  n_pred_test <- length(ind_last_is_test)

  ## Herds for which a test was NOT performed on the month to predict
  ind_last_is_not_test <- which(data$herd_test_data$last_is_test == 0)
  ind_p_no_test <- data$herd_test_data$ind_p[ind_last_is_not_test]
  n_pred_no_test <- length(ind_last_is_not_test)

  ## test results used for prediction
  test_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,test_res_col]), test_res_col]
  test_id_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,test_res_col]), "test_id"]

  ## Priors for risk factors
  risk_factors <- data$risk_factors[data$risk_factors$ref == 0,]
  theta_norm_mean <- risk_factors$mean_prior
  theta_norm_prec <- 1 / risk_factors$sd_prior^2

  ## Priors for test characteristics
  test_char <- data$test_perf_prior

  ## Data used by JAGS
  JAGS_data <- list(
    n_herds = n_herds,
    ind_i = data$herd_test_data$ind_i,
    ind_j = data$herd_test_data$ind_j,
    ind_f = data$herd_test_data$ind_f,
    ind_p = data$herd_test_data$ind_p,
    n_tests_perf = length(ind_test),
    ind_test = as.integer(ind_test),
    test_res = as.integer(test_res),
    test_id = as.integer(test_id),
    test_for_pred = test_for_pred,
    test_id_for_pred = test_id_for_pred,
    n_pred_test = n_pred_test,
    ind_last_is_test = ind_last_is_test,
    n_pred_no_test = n_pred_no_test,
    ind_p_test = ind_p_test,
    ind_last_is_not_test = ind_last_is_not_test,
    ind_p_no_test = ind_p_no_test,
    pi1_beta_a = data$inf_dyn_priors["pi1_a"],
    pi1_beta_b = data$inf_dyn_priors["pi1_b"],
    tau2_beta_a = data$inf_dyn_priors["tau2_a"],
    tau2_beta_b = data$inf_dyn_priors["tau2_b"],
    n_tests = nrow(test_char),
    Se_beta_a = test_char$Se_a,
    Se_beta_b = test_char$Se_b,
    Sp_beta_a = test_char$Sp_a,
    Sp_beta_b = test_char$Sp_b,
    theta_norm_mean = theta_norm_mean,
    theta_norm_prec = theta_norm_prec,
    n_risk_factors = length(theta_norm_mean),
    risk_factors = as.matrix(data$risk_factor_data[, -(1:3)])
  )

  write_JAGS_model(data)

  JAGS_model_compiled <- rjags::jags.model(
    file = "JAGS_model.txt",
    data = JAGS_data,
    n.chains = n_chains)

  if (file.exists("JAGS_model.txt")) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "herd_ntests_rf")

  JAGS_model_compiled

}

#' @export
compile_JAGS.animal_1test = function(data, n_chains){

  n_herds <- attr(data, "number of herds")
  month_max <- max(data$test_data$month_id)

  ## test results
  n_pos    <- data$test_data[data$test_data$month_id < month_max, "n_pos"]
  n_tested <- data$test_data[data$test_data$month_id < month_max, "n_tested"]
  ## Indices when test performed in the expanded dataset
  ind_test <- data$test_data$status_id[data$test_data$month_id < month_max]

  ## Herds for which a test was performed on the month to predict
  ind_last_is_test <- which(data$herd_test_data$last_is_test == 1)
  ind_p_test <- data$herd_test_data$ind_p[ind_last_is_test]
  n_pred_test <- length(ind_last_is_test)

  ## Herds for which a test was NOT performed on the month to predict
  ind_last_is_not_test <- which(data$herd_test_data$last_is_test == 0)
  ind_p_no_test <- data$herd_test_data$ind_p[ind_last_is_not_test]
  n_pred_no_test <- length(ind_last_is_not_test)

  ## test results used for prediction
  n_pos_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,"n_tested"]), "n_pos"]
  n_tested_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,"n_tested"]), "n_tested"]
  test_id_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,"n_tested"]), "test_id"]

  ## Priors for test characteristics
  test_char <- data$test_perf_prior

  ## Data used by JAGS
  JAGS_data <- list(
    n_herds = n_herds,
    ind_i = data$herd_test_data$ind_i,
    ind_j = data$herd_test_data$ind_j,
    ind_f = data$herd_test_data$ind_f,
    n_tests_perf = length(ind_test),
    ind_test = as.integer(ind_test),
    n_pos = as.integer(n_pos),
    n_tested = as.integer(n_tested),
    n_pos_for_pred = n_pos_for_pred,
    n_tested_for_pred = n_tested_for_pred,
    n_pred_test = n_pred_test,
    ind_last_is_test = ind_last_is_test,
    n_pred_no_test = n_pred_no_test,
    ind_p_test = ind_p_test,
    ind_last_is_not_test = ind_last_is_not_test,
    ind_p_no_test = ind_p_no_test,
    pi1_beta_a = data$inf_dyn_priors["pi1_a"],
    pi1_beta_b = data$inf_dyn_priors["pi1_b"],
    pi_within_a = data$inf_dyn_priors["pi_within_a"],
    pi_within_b = data$inf_dyn_priors["pi_within_b"],
    tau1_beta_a = data$inf_dyn_priors["tau1_a"],
    tau1_beta_b = data$inf_dyn_priors["tau1_b"],
    tau2_beta_a = data$inf_dyn_priors["tau2_a"],
    tau2_beta_b = data$inf_dyn_priors["tau2_b"],
    Se_beta_a = test_char$Se_a,
    Se_beta_b = test_char$Se_b,
    Sp_beta_a = test_char$Sp_a,
    Sp_beta_b = test_char$Sp_b
  )

  write_JAGS_model(data)

  JAGS_model_compiled <- rjags::jags.model(
    file = "JAGS_model.txt",
    data = JAGS_data,
    n.chains = n_chains)

  if (file.exists("JAGS_model.txt")) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "animal_1test")

  JAGS_model_compiled

}


#' @export
compile_JAGS.animal_1test_rf = function(data, n_chains){

  n_herds <- attr(data, "number of herds")
  month_max <- max(data$test_data$month_id)

  ## test results
  n_pos    <- data$test_data[data$test_data$month_id < month_max, "n_pos"]
  n_tested <- data$test_data[data$test_data$month_id < month_max, "n_tested"]
  ## Indices when test performed in the expanded dataset
  ind_test <- data$test_data$status_id[data$test_data$month_id < month_max]

  ## Herds for which a test was performed on the month to predict
  ind_last_is_test <- which(data$herd_test_data$last_is_test == 1)
  ind_p_test <- data$herd_test_data$ind_p[ind_last_is_test]
  n_pred_test <- length(ind_last_is_test)

  ## Herds for which a test was NOT performed on the month to predict
  ind_last_is_not_test <- which(data$herd_test_data$last_is_test == 0)
  ind_p_no_test <- data$herd_test_data$ind_p[ind_last_is_not_test]
  n_pred_no_test <- length(ind_last_is_not_test)

  ## test results used for prediction
  n_pos_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,"n_tested"]), "n_pos"]
  n_tested_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,"n_tested"]), "n_tested"]

  ## Priors for risk factors
  risk_factors <- data$risk_factors[data$risk_factors$ref == 0,]
  theta_norm_mean <- risk_factors$mean_prior
  theta_norm_prec <- 1 / risk_factors$sd_prior^2

  ## Priors for test characteristics
  test_char <- data$test_perf_prior

  ## Data used by JAGS
  JAGS_data <- list(
    n_herds = n_herds,
    ind_i = data$herd_test_data$ind_i,
    ind_j = data$herd_test_data$ind_j,
    ind_f = data$herd_test_data$ind_f,
    ind_p = data$herd_test_data$ind_p,
    n_tests_perf = length(ind_test),
    ind_test = as.integer(ind_test),
    n_pos = as.integer(n_pos),
    n_tested = as.integer(n_tested),
    n_pos_for_pred = n_pos_for_pred,
    n_tested_for_pred = n_tested_for_pred,
    n_pred_test = n_pred_test,
    ind_last_is_test = ind_last_is_test,
    n_pred_no_test = n_pred_no_test,
    ind_p_test = ind_p_test,
    ind_last_is_not_test = ind_last_is_not_test,
    ind_p_no_test = ind_p_no_test,
    pi1_beta_a = data$inf_dyn_priors["pi1_a"],
    pi1_beta_b = data$inf_dyn_priors["pi1_b"],
    pi_within_a = data$inf_dyn_priors["pi_within_a"],
    pi_within_b = data$inf_dyn_priors["pi_within_b"],
    tau2_beta_a = data$inf_dyn_priors["tau2_a"],
    tau2_beta_b = data$inf_dyn_priors["tau2_b"],
    Se_beta_a = test_char$Se_a,
    Se_beta_b = test_char$Se_b,
    Sp_beta_a = test_char$Sp_a,
    Sp_beta_b = test_char$Sp_b,
    theta_norm_mean = theta_norm_mean,
    theta_norm_prec = theta_norm_prec,
    n_risk_factors = length(theta_norm_mean),
    risk_factors = as.matrix(data$risk_factor_data[, -(1:3)])
  )

  write_JAGS_model(data)

  JAGS_model_compiled <- rjags::jags.model(
    file = "JAGS_model.txt",
    data = JAGS_data,
    n.chains = n_chains)

  if (file.exists("JAGS_model.txt")) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "animal_1test_rf")

  JAGS_model_compiled

}

#' @export
compile_JAGS.animal_ntests_rf = function(data, n_chains){

  n_herds <- attr(data, "number of herds")
  month_max <- max(data$test_data$month_id)

  ## test results
  n_pos    <- data$test_data[data$test_data$month_id < month_max, "n_pos"]
  n_tested <- data$test_data[data$test_data$month_id < month_max, "n_tested"]
  test_id  <- data$test_data$test_id[data$test_data$month_id < month_max]
  ## Indices when test performed in the expanded dataset
  ind_test <- data$test_data$status_id[data$test_data$month_id < month_max]

  ## Herds for which a test was performed on the month to predict
  ind_last_is_test <- which(data$herd_test_data$last_is_test == 1)
  ind_p_test <- data$herd_test_data$ind_p[ind_last_is_test]
  n_pred_test <- length(ind_last_is_test)

  ## Herds for which a test was NOT performed on the month to predict
  ind_last_is_not_test <- which(data$herd_test_data$last_is_test == 0)
  ind_p_no_test <- data$herd_test_data$ind_p[ind_last_is_not_test]
  n_pred_no_test <- length(ind_last_is_not_test)

  ## test results used for prediction
  n_pos_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,"n_tested"]), "n_pos"]
  n_tested_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,"n_tested"]), "n_tested"]
  test_id_for_pred <- data$test_data[data$test_data$month_id == month_max & !is.na(data$test_data[,"n_tested"]), "test_id"]

  ## Priors for risk factors
  risk_factors <- data$risk_factors[data$risk_factors$ref == 0,]
  theta_norm_mean <- risk_factors$mean_prior
  theta_norm_prec <- 1 / risk_factors$sd_prior^2

  ## Priors for test characteristics
  test_char <- data$test_perf_prior

  ## Data used by JAGS
  JAGS_data <- list(
    n_herds = n_herds,
    ind_i = data$herd_test_data$ind_i,
    ind_j = data$herd_test_data$ind_j,
    ind_f = data$herd_test_data$ind_f,
    ind_p = data$herd_test_data$ind_p,
    n_tests_perf = length(ind_test),
    ind_test = as.integer(ind_test),
    n_pos = as.integer(n_pos),
    n_tested = as.integer(n_tested),
    test_id = as.integer(test_id),
    n_pos_for_pred = n_pos_for_pred,
    n_tested_for_pred = n_tested_for_pred,
    test_id_for_pred = as.integer(test_id_for_pred),
    n_pred_test = n_pred_test,
    ind_last_is_test = ind_last_is_test,
    n_pred_no_test = n_pred_no_test,
    ind_p_test = ind_p_test,
    ind_last_is_not_test = ind_last_is_not_test,
    ind_p_no_test = ind_p_no_test,
    pi1_beta_a = data$inf_dyn_priors["pi1_a"],
    pi1_beta_b = data$inf_dyn_priors["pi1_b"],
    pi_within_a = data$inf_dyn_priors["pi_within_a"],
    pi_within_b = data$inf_dyn_priors["pi_within_b"],
    tau2_beta_a = data$inf_dyn_priors["tau2_a"],
    tau2_beta_b = data$inf_dyn_priors["tau2_b"],
    Se_beta_a = test_char$Se_a,
    Se_beta_b = test_char$Se_b,
    Sp_beta_a = test_char$Sp_a,
    Sp_beta_b = test_char$Sp_b,
    theta_norm_mean = theta_norm_mean,
    theta_norm_prec = theta_norm_prec,
    n_risk_factors = length(theta_norm_mean),
    risk_factors = as.matrix(data$risk_factor_data[, -(1:3)])
  )

  write_JAGS_model(data)

  JAGS_model_compiled <- rjags::jags.model(
    file = "JAGS_model.txt",
    data = JAGS_data,
    n.chains = n_chains)

  if (file.exists("JAGS_model.txt")) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "animal_ntests_rf")

  JAGS_model_compiled

}
