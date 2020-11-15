#' Compilation of the JAGS model
#'
#' @param data an object of class STOCfree_data
#' @param n_chains number of MCMC chains
#' @param keep_model_file if TRUE, the JAGS model is written as a text file in the working folder
#' @param status 'discrete' latent status (0/1) is modelled, 'proba' probability of latent status is modelled, 'predict' predictive simulation from the model prior and data
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
compile_JAGS.default <- function(data,
                                 n_chains,
                                 keep_model_file){

  print("No method defined for this type of data")

  }

#' @export
compile_JAGS.herd <- function(data,
                              n_chains = 4,
                              keep_model_file = FALSE){

  n_herds <- attr(data, "number of herds")
  month_max <- max(data$test_data$month_id)

  test_data <- data$test_data
  ## no test
  status_no_test <- data$risk_factor_data$status_id
  status_no_test <- sort(status_no_test[-which(status_no_test %in% test_data$status_id)])

  ## number of tests of optional types
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ4 <- nrow(test_data[test_data$status_type == 4,])
  n_status_typ5 <- nrow(test_data[test_data$status_type == 5,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])

  ## Priors for risk factors
  risk_factors <- data$risk_factors[data$risk_factors$ref == 0,]
  theta_norm_mean <- data$risk_factors$mean_prior
  theta_norm_prec <- 1 / data$risk_factors$sd_prior^2

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

    if(n_status_typ3 > 0){
    ## 3: test > 1 on a month
    JAGS_data$n_status_typ3 <- n_status_typ3
    JAGS_data$status_typ3 <- test_data$status_id[test_data$status_type == 3]
    JAGS_data$test_res_typ3 <- test_data$test_res[test_data$status_type == 3]
    JAGS_data$test_id_typ3 <- test_data$test_id[test_data$status_type == 3]
  }

  if(n_status_typ4 > 0){
    ## 4: status to predict without test result
    JAGS_data$n_status_typ4 <- n_status_typ4
    JAGS_data$status_typ4 <- test_data$status_id[test_data$status_type == 4]
    JAGS_data$herd_id_pr4 <- test_data$herd_id[test_data$status_type == 4]
  }

  if(n_status_typ5 > 0){
    ## 5: status to predict with a single test performed
    JAGS_data$n_status_typ5 <- n_status_typ5
    JAGS_data$status_typ5 <- test_data$status_id[test_data$status_type == 5]
    JAGS_data$test_res_typ5 <- test_data$test_res[test_data$status_type == 5]
    JAGS_data$test_id_typ5 <- test_data$test_id[test_data$status_type == 5]
    JAGS_data$herd_id_pr5 <- test_data$herd_id[test_data$status_type == 5]
  }

  if(n_status_typ6 > 0){
    ## 6: status to predict with several tests on this month
    JAGS_data$n_status_typ6 <- n_status_typ6
    JAGS_data$status_typ6   <- test_data$status_id[test_data$status_type == 6]
    JAGS_data$test_res_typ6 <- test_data$test_res[test_data$status_type == 6]
    JAGS_data$test_id_typ6  <- test_data$test_id[test_data$status_type == 6]
    JAGS_data$herd_id_pr6 <- test_data$herd_id[test_data$status_type == 6]
  }

  write_JAGS_model(data)

  JAGS_model_compiled <- rjags::jags.model(
    file = "JAGS_model.txt",
    data = JAGS_data,
    n.chains = n_chains)

  if(file.exists("JAGS_model.txt") & isFALSE(keep_model_file)) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "herd")

  JAGS_model_compiled

}


#' @export
compile_JAGS.herd_rf <- function(data,
                                 n_chains = 4,
                                 keep_model_file = FALSE,
                                 status = c("discrete", "proba", "predict")){

  n_herds <- attr(data, "number of herds")
  month_max <- max(data$test_data$month_id)

  test_data <- data$test_data
  ## no test
  status_no_test <- data$risk_factor_data$status_id
  status_no_test <- sort(status_no_test[-which(status_no_test %in% test_data$status_id)])

  ## number of tests of optional types
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ4 <- nrow(test_data[test_data$status_type == 4,])
  n_status_typ5 <- nrow(test_data[test_data$status_type == 5,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])

  ## Priors for risk factors
  risk_factors <- data$risk_factors[data$risk_factors$ref == 0,]
  theta_norm_mean <- data$risk_factors$mean_prior
  theta_norm_prec <- 1 / data$risk_factors$sd_prior^2

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

  if(n_status_typ3 > 0){
  ## 3: test > 1 on a month
    JAGS_data$n_status_typ3 <- n_status_typ3
    JAGS_data$status_typ3 <- test_data$status_id[test_data$status_type == 3]
    JAGS_data$test_res_typ3 <- test_data$test_res[test_data$status_type == 3]
    JAGS_data$test_id_typ3 <- test_data$test_id[test_data$status_type == 3]
  }

  if(n_status_typ4 > 0){
    ## 4: status to predict without test result
    JAGS_data$n_status_typ4 <- n_status_typ4
    JAGS_data$status_typ4 <- test_data$status_id[test_data$status_type == 4]
    JAGS_data$herd_id_pr4 <- test_data$herd_id[test_data$status_type == 4]
  }

  if(n_status_typ5 > 0){
    ## 5: status to predict with a single test performed
    JAGS_data$n_status_typ5 <- n_status_typ5
    JAGS_data$status_typ5 <- test_data$status_id[test_data$status_type == 5]
    JAGS_data$test_res_typ5 <- test_data$test_res[test_data$status_type == 5]
    JAGS_data$test_id_typ5 <- test_data$test_id[test_data$status_type == 5]
    JAGS_data$herd_id_pr5 <- test_data$herd_id[test_data$status_type == 5]
  }

  if(n_status_typ6 > 0){
    ## 6: status to predict with several tests on this month
    JAGS_data$n_status_typ6 <- n_status_typ6
    JAGS_data$status_typ6   <- test_data$status_id[test_data$status_type == 6]
    JAGS_data$test_res_typ6 <- test_data$test_res[test_data$status_type == 6]
    JAGS_data$test_id_typ6  <- test_data$test_id[test_data$status_type == 6]
    JAGS_data$herd_id_pr6 <- test_data$herd_id[test_data$status_type == 6]
  }

  write_JAGS_model(data)

  JAGS_model_compiled <- rjags::jags.model(
    file = "JAGS_model.txt",
    data = JAGS_data,
    n.chains = n_chains)

  if(file.exists("JAGS_model.txt") & isFALSE(keep_model_file)) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "herd_rf")

  JAGS_model_compiled

}


#' @export
compile_JAGS.animal <- function(data,
                               n_chains = 4,
                               keep_model_file = FALSE,
                               status = c("discrete", "proba", "predict")){

  n_herds <- attr(data, "number of herds")
  month_max <- max(data$test_data$month_id)

  test_data <- data$test_data
  ## no test
  status_no_test <- data$risk_factor_data$status_id
  status_no_test <- sort(status_no_test[-which(status_no_test %in% test_data$status_id)])

  ## number of tests of optional types
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ4 <- nrow(test_data[test_data$status_type == 4,])
  n_status_typ5 <- nrow(test_data[test_data$status_type == 5,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])

  ## Priors for test characteristics
  test_char <- data$test_perf_prior

  ## Data used by JAGS
  JAGS_data <- list(
    ## First test in a herd - status type = 1
    n_status_typ1 = nrow(test_data[test_data$status_type == 1,]),
    status_typ1 = test_data$status_id[test_data$status_type == 1],
    n_pos_typ1 = test_data$n_pos[test_data$status_type == 1],
    n_tested_typ1 = test_data$n_tested[test_data$status_type == 1],
    test_id_typ1  = test_data$test_id[test_data$status_type == 1],
    ## 2: first test on a month which is not first test in herd
    n_status_typ2 = nrow(test_data[test_data$status_type == 2,]),
    status_typ2 = test_data$status_id[test_data$status_type == 2],
    n_pos_typ2 = test_data$n_pos[test_data$status_type == 2],
    n_tested_typ2 = test_data$n_tested[test_data$status_type == 2],
    test_id_typ2  = test_data$test_id[test_data$status_type == 2],
    ## no test
    n_no_test = length(status_no_test),
    status_no_test = status_no_test,
    pi1_beta_a = data$inf_dyn_priors["pi1_a"],
    pi1_beta_b = data$inf_dyn_priors["pi1_b"],
    pi_within_a = data$inf_dyn_priors["pi_within_a"],
    pi_within_b = data$inf_dyn_priors["pi_within_b"],
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

  if(n_status_typ3 > 0){
    ## 3: test > 1 on a month
    JAGS_data$n_status_typ3 <- n_status_typ3
    JAGS_data$status_typ3 <- test_data$status_id[test_data$status_type == 3]
    JAGS_data$n_pos_typ3 <- test_data$n_pos[test_data$status_type == 3]
    JAGS_data$n_tested_typ3 <- test_data$n_tested[test_data$status_type == 3]
    JAGS_data$test_id_typ3 <- test_data$test_id[test_data$status_type == 3]
  }

  if(n_status_typ4 > 0){
    ## 4: status to predict without test result
    JAGS_data$n_status_typ4 <- n_status_typ4
    JAGS_data$status_typ4 <- test_data$status_id[test_data$status_type == 4]
    JAGS_data$herd_id_pr4 <- test_data$herd_id[test_data$status_type == 4]
  }

  if(n_status_typ5 > 0){
    ## 5: status to predict with a single test performed
    JAGS_data$n_status_typ5 <- n_status_typ5
    JAGS_data$status_typ5   <- test_data$status_id[test_data$status_type == 5]
    JAGS_data$n_pos_typ5    <- test_data$n_pos[test_data$status_type == 5]
    JAGS_data$n_tested_typ5 <- test_data$n_tested[test_data$status_type == 5]
    JAGS_data$test_id_typ5  <- test_data$test_id[test_data$status_type == 5]
    JAGS_data$herd_id_pr5   <- test_data$herd_id[test_data$status_type == 5]
  }

  if(n_status_typ6 > 0){
    ## 6: status to predict with several tests on this month
    JAGS_data$n_status_typ6 <- n_status_typ6
    JAGS_data$status_typ6 <- test_data$status_id[test_data$status_type == 6]
    JAGS_data$n_pos_typ6 <- test_data$n_pos[test_data$status_type == 6]
    JAGS_data$n_tested_typ6 <- test_data$n_tested[test_data$status_type == 6]
    JAGS_data$test_id_typ6 <- test_data$test_id[test_data$status_type == 6]
    JAGS_data$herd_id_pr6 <- test_data$herd_id[test_data$status_type == 6]
  }

  write_JAGS_model(data)

  JAGS_model_compiled <- rjags::jags.model(
    file = "JAGS_model.txt",
    data = JAGS_data,
    n.chains = n_chains)

  if(file.exists("JAGS_model.txt") & isFALSE(keep_model_file)) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "animal")

  JAGS_model_compiled

}


#' @export
compile_JAGS.animal_rf <- function(data,
                                  n_chains = 4,
                                  keep_model_file = FALSE,
                                  status = c("discrete", "proba", "predict")){

  n_herds <- attr(data, "number of herds")
  month_max <- max(data$test_data$month_id)

  test_data <- data$test_data
  ## no test
  status_no_test <- data$risk_factor_data$status_id
  status_no_test <- sort(status_no_test[-which(status_no_test %in% test_data$status_id)])

  ## number of tests of optional types
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ4 <- nrow(test_data[test_data$status_type == 4,])
  n_status_typ5 <- nrow(test_data[test_data$status_type == 5,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])

  ## Priors for risk factors
  risk_factors <- data$risk_factors[data$risk_factors$ref == 0,]
  theta_norm_mean <- data$risk_factors$mean_prior
  theta_norm_prec <- 1 / data$risk_factors$sd_prior^2

  ## Priors for test characteristics
  test_char <- data$test_perf_prior

  ## Data used by JAGS
  JAGS_data <- list(
    ## First test in a herd - status type = 1
    n_status_typ1 = nrow(test_data[test_data$status_type == 1,]),
    status_typ1 = test_data$status_id[test_data$status_type == 1],
    n_pos_typ1 = test_data$n_pos[test_data$status_type == 1],
    n_tested_typ1 = test_data$n_tested[test_data$status_type == 1],
    test_id_typ1  = test_data$test_id[test_data$status_type == 1],
    ## 2: first test on a month which is not first test in herd
    n_status_typ2 = nrow(test_data[test_data$status_type == 2,]),
    status_typ2 = test_data$status_id[test_data$status_type == 2],
    n_pos_typ2 = test_data$n_pos[test_data$status_type == 2],
    n_tested_typ2 = test_data$n_tested[test_data$status_type == 2],
    test_id_typ2  = test_data$test_id[test_data$status_type == 2],
    ## no test
    n_no_test = length(status_no_test),
    status_no_test = status_no_test,
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

  if(n_status_typ3 > 0){
    ## 3: test > 1 on a month
    JAGS_data$n_status_typ3 <- n_status_typ3
    JAGS_data$status_typ3 <- test_data$status_id[test_data$status_type == 3]
    JAGS_data$n_pos_typ3 <- test_data$n_pos[test_data$status_type == 3]
    JAGS_data$n_tested_typ3 <- test_data$n_tested[test_data$status_type == 3]
    JAGS_data$test_id_typ3 <- test_data$test_id[test_data$status_type == 3]
  }

  if(n_status_typ4 > 0){
    ## 4: status to predict without test result
    JAGS_data$n_status_typ4 <- n_status_typ4
    JAGS_data$status_typ4 <- test_data$status_id[test_data$status_type == 4]
    JAGS_data$herd_id_pr4 <- test_data$herd_id[test_data$status_type == 4]
  }

  if(n_status_typ5 > 0){
    ## 5: status to predict with a single test performed
    JAGS_data$n_status_typ5 <- n_status_typ5
    JAGS_data$status_typ5 <- test_data$status_id[test_data$status_type == 5]
    JAGS_data$n_pos_typ5 <- test_data$n_pos[test_data$status_type == 5]
    JAGS_data$n_tested_typ5 <- test_data$n_tested[test_data$status_type == 5]
    JAGS_data$test_id_typ5 <- test_data$test_id[test_data$status_type == 5]
    JAGS_data$herd_id_pr5 <- test_data$herd_id[test_data$status_type == 5]
  }

  if(n_status_typ6 > 0){
    ## 6: status to predict with several tests on this month
    JAGS_data$n_status_typ6 <- n_status_typ6
    JAGS_data$status_typ6   <- test_data$status_id[test_data$status_type == 6]
    JAGS_data$n_pos_typ6 <- test_data$n_pos[test_data$status_type == 6]
    JAGS_data$n_tested_typ6 <- test_data$n_tested[test_data$status_type == 6]
    JAGS_data$test_id_typ6  <- test_data$test_id[test_data$status_type == 6]
    JAGS_data$herd_id_pr6 <- test_data$herd_id[test_data$status_type == 6]
  }

  write_JAGS_model(data)

  JAGS_model_compiled <- rjags::jags.model(
    file = "JAGS_model.txt",
    data = JAGS_data,
    n.chains = n_chains)

  if(file.exists("JAGS_model.txt") & isFALSE(keep_model_file)) file.remove("JAGS_model.txt")

  class(JAGS_model_compiled) <- c("jags", "animal_rf")

  JAGS_model_compiled

}
