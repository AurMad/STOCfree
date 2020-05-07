#' Compilation of the JAGS model
#'
#' @param data an object of class STOCfree_data
#' @param n_chains number of MCMC chains
#' @param keep_model should the JAGS model file be saved on the disk
#' @param ...
#'
#' @return
#' a compiled JAGS model
#'
#' @export
#'
#' @examples
compile_JAGS <- function(data,
                         n_chains = 4,
                         keep_model_file = FALSE,
                         ...){

  UseMethod("compile_JAGS")

  }


#' @export
compile_JAGS.default <- function(data, n_chains, keep_model_file){

  print("No method defined for this type of data")

  }

compile_JAGS.herd <- function(data, n_chains, keep_model_file){

  n_herds <- attr(data, "number of herds")
  month_max <- max(data$test_data$month_id)

  test_data <- data$test_data
  ## no test
  status_no_test <- data$risk_factor_data$status_id
  status_no_test <- sort(status_no_test[-which(status_no_test %in% test_data$status_id)])

  ## Priors for risk factors
  risk_factors <- data$risk_factors[data$risk_factors$ref == 0,]
  theta_norm_mean <- risk_factors$mean_prior
  theta_norm_prec <- 1 / risk_factors$sd_prior^2

  ## Priors for test characteristics
  test_char <- data$test_perf_prior

  ## Data used by JAGS
  JAGS_data <- list(
    ## First test in a herd - status type = 1
    n_status_typ1 = nrow(test_data[test_data$status_type == 1,]),
    status_typ1 = test_data$status_id[test_data$status_type == 1],
    test_res_typ1 = test_data$test_res[test_data$status_type == 1],
    test_id_typ1  = test_data$test_id[test_data$status_type == 1],
    ## 2: first test on a month which is not first test in herd
    n_status_typ2 = nrow(test_data[test_data$status_type == 2,]),
    status_typ2 = test_data$status_id[test_data$status_type == 2],
    test_res_typ2 = test_data$test_res[test_data$status_type == 2],
    test_id_typ2  = test_data$test_id[test_data$status_type == 2],
    ## 3: test > 1 on a month
    n_status_typ3 = nrow(test_data[test_data$status_type == 3,]),
    status_typ3 = test_data$status_id[test_data$status_type == 3],
    test_res_typ3 = test_data$test_res[test_data$status_type == 3],
    test_id_typ3 = test_data$test_id[test_data$status_type == 3],
    ## no test
    n_no_test = length(status_no_test),
    status_no_test = status_no_test,
    ## 4: status to predict without test result
    n_status_typ4 = nrow(test_data[test_data$status_type == 4,]),
    status_typ4 = test_data$status_id[test_data$status_type == 4],
    herd_id_pr4 = test_data$herd_id[test_data$status_type == 4],
    ## 5: status to predict with a single test performed
    n_status_typ5 = nrow(test_data[test_data$status_type == 5,]),
    status_typ5 = test_data$status_id[test_data$status_type == 5],
    test_res_typ5 = test_data$test_res[test_data$status_type == 5],
    test_id_typ5  = test_data$test_id[test_data$status_type == 5],
    herd_id_pr5 = test_data$herd_id[test_data$status_type == 5],
    ## 6: status to predict with several tests on this month
    n_status_typ6 = nrow(test_data[test_data$status_type == 6,]),
    status_typ6   = test_data$status_id[test_data$status_type == 6],
    test_res_typ6 = test_data$test_res[test_data$status_type == 6],
    test_id_typ6  = test_data$test_id[test_data$status_type == 6],
    herd_id_pr6 = test_data$herd_id[test_data$status_type == 6],
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

  if(file.exists("JAGS_model.txt") & isFALSE(keep_model_file)) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "herd_ntests")

  JAGS_model_compiled

}


#' @export
compile_JAGS.herd_rf <- function(data, n_chains, keep_model_file){

  n_herds <- attr(data, "number of herds")
  month_max <- max(data$test_data$month_id)

  test_data <- data$test_data
  ## no test
  status_no_test <- data$risk_factor_data$status_id
  status_no_test <- sort(status_no_test[-which(status_no_test %in% test_data$status_id)])

  ## several tests on the same month
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])

  ## Priors for risk factors
  risk_factors <- data$risk_factors[data$risk_factors$ref == 0,]
  theta_norm_mean <- risk_factors$mean_prior
  theta_norm_prec <- 1 / risk_factors$sd_prior^2

  ## Priors for test characteristics
  test_char <- data$test_perf_prior

  ## Data used by JAGS
  JAGS_data <- list(
    ## First test in a herd - status type = 1
    n_status_typ1 = nrow(test_data[test_data$status_type == 1,]),
    status_typ1 = test_data$status_id[test_data$status_type == 1],
    test_res_typ1 = test_data$test_res[test_data$status_type == 1],
    test_id_typ1  = test_data$test_id[test_data$status_type == 1],
    ## 2: first test on a month which is not first test in herd
    n_status_typ2 = nrow(test_data[test_data$status_type == 2,]),
    status_typ2 = test_data$status_id[test_data$status_type == 2],
    test_res_typ2 = test_data$test_res[test_data$status_type == 2],
    test_id_typ2  = test_data$test_id[test_data$status_type == 2],
    ## no test
    n_no_test = length(status_no_test),
    status_no_test = status_no_test,
    ## 4: status to predict without test result
    n_status_typ4 = nrow(test_data[test_data$status_type == 4,]),
    status_typ4 = test_data$status_id[test_data$status_type == 4],
    herd_id_pr4 = test_data$herd_id[test_data$status_type == 4],
    ## 5: status to predict with a single test performed
    n_status_typ5 = nrow(test_data[test_data$status_type == 5,]),
    status_typ5 = test_data$status_id[test_data$status_type == 5],
    test_res_typ5 = test_data$test_res[test_data$status_type == 5],
    test_id_typ5  = test_data$test_id[test_data$status_type == 5],
    herd_id_pr5 = test_data$herd_id[test_data$status_type == 5],
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
    risk_factors = as.matrix(data$risk_factor_data[, -(1:4)])
  )

  if(n_status_typ3 > 0){
  ## 3: test > 1 on a month
  JAGS_data <- c(JAGS_data,
                   n_status_typ3 = n_status_typ3,
                   status_typ3 = test_data$status_id[test_data$status_type == 3],
                   test_res_typ3 = test_data$test_res[test_data$status_type == 3],
                   test_id_typ3 = test_data$test_id[test_data$status_type == 3])
  }


  if(n_status_typ6 > 0){
    ## 6: status to predict with several tests on this month
  JAGS_data <- c(JAGS_data,
                 n_status_typ6 = n_status_typ6,
                 status_typ6   = test_data$status_id[test_data$status_type == 6],
                 test_res_typ6 = test_data$test_res[test_data$status_type == 6],
                 test_id_typ6  = test_data$test_id[test_data$status_type == 6],
                 herd_id_pr6 = test_data$herd_id[test_data$status_type == 6])
  }

  write_JAGS_model(data)

  JAGS_model_compiled <- rjags::jags.model(
    file = "JAGS_model.txt",
    data = JAGS_data,
    n.chains = n_chains)

  if(file.exists("JAGS_model.txt") & isFALSE(keep_model_file)) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "herd_ntests_rf")

  JAGS_model_compiled

}


#' @export
compile_JAGS.animal_rf = function(data, n_chains, keep_model_file){

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

  if(file.exists("JAGS_model.txt") & isFALSE(keep_model_file)) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "animal_ntests_rf")

  JAGS_model_compiled

}
