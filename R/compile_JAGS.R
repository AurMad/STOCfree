#' Compilation of the JAGS model
#'
#' @param test_data a tibble or a data.frame with one row per month
#' @param herd_id name of the column with herd identifier
#' @param row_id name of the column with the row IDs
#' @param month name of the colummn with month number
#' @param test_result name of the colummn with test results
#' @param risk_factors character string for the risk factors to include
#' @param test_priors list of parameters for the priors for test characteristics
#' @param infection_priors list of parameters for the infection priors
#' @param risk_factor_priors list of parameters for the prior for the risk factor associations
#' @param n_chains number of MCMC chains
#'
#' @return
#' a compiled JAGS model
#'
#' @export
#'
#' @examples
compile_JAGS <- function(test_data, herd_id, row_id, month,
                         risk_factors,
                         test_result,
                         test_priors, infection_priors, risk_factor_priors,
                         n_chains){

  herd   <- dplyr::enquo(herd_id)
  mnthid <- dplyr::enquo(month)
  rowid  <- dplyr::enquo(row_id)
  test   <- dplyr::enquo(test_result)

  ## Matrix of risk factors
  X <- matrix(rep(1, nrow(test_data)))


  for(i in 1:length(risk_factors)){

    cln <- match(risk_factors[i], colnames(test_data))
    Z <- test_data[, cln]
    Z[is.na(Z)] <- 0
    X <- cbind(X, Z)
  }



  ## Indices for herd level data
  herd <- group_by_herd(data = test_data,
                         herd_colname = !! herd,
                         month_colname = !! mnthid,
                         row_id_colname = !! rowid,
                         test_res_colname = !! test)

  ## Herds for which a test was performed on the month to predict
  ind_last_is_test <- which(herd$last_is_test == 1)
  n_pred_test <- length(ind_last_is_test)


  ## Herds for which a test was NOT performed on the month to predict
  ind_last_is_not_test <- which(herd$last_is_test == 0)
  n_pred_no_test <- length(ind_last_is_not_test)

  ## Indices when test performed in the expanded dataset
  ind_test <- dplyr::filter(test_data, !is.na(!! test) & !! mnthid < max(!! mnthid))
  ind_test <- unlist(dplyr::select(ind_test, !! rowid))

  pi1_beta_a <- 1
  pi1_beta_b <- 2

  tau2_beta_a <- 30
  tau2_beta_b <- 2

  Se_beta_a <- 12
  Se_beta_b <- 2

  Sp_beta_a <- 200
  Sp_beta_b <- 4

  theta1_norm_mean <- -3.8
  theta1_norm_sd <- .15

  theta2_norm_mean <- -.6
  theta2_norm_sd <- 1.2

  risk_factors <- matrix(
    rep(0, nrow(test_data) * 2),
    ncol = 2
  )

  ## Data used by JAGS
  JAGS_data <- list(
    n_herds = nrow(herd),
    ind_i = herd$ind_i,
    ind_j = herd$ind_j,
    ind_f = herd$ind_f,
    ind_p = herd$ind_p,
    n_tests = length(ind_test),
    ind_test = as.numeric(ind_test),
    test_res = as.numeric(unlist(dplyr::select(test_data, !! test))),
    n_pred_test = n_pred_test,
    ind_last_is_test = ind_last_is_test,
    n_pred_no_test = n_pred_no_test,
    ind_last_is_not_test = ind_last_is_not_test,
    pi1_beta_a = infection_priors$pi1_beta_a,
    pi1_beta_b = infection_priors$pi1_beta_b,
    tau2_beta_a = infection_priors$tau2_beta_a,
    tau2_beta_b = infection_priors$tau2_beta_b,
    Se_beta_a = test_priors$Se_beta_a,
    Se_beta_b = test_priors$Se_beta_b,
    Sp_beta_a = test_priors$Sp_beta_a,
    Sp_beta_b = test_priors$Sp_beta_b,
    theta_norm_mean = risk_factor_priors$theta_norm_mean,
    theta_norm_prec = 1 / risk_factor_priors$theta_norm_sd^2,
    n_risk_factors = dim(X)[2],
    risk_factors = X
  )

  write_JAGS_model()

  JAGS_model_compiled <- rjags::jags.model(
    file = "JAGS_model.txt",
    data = JAGS_data,
    n.chains = n_chains)

  if (file.exists("JAGS_model.txt")) file.remove("JAGS_model.txt")

  JAGS_model_compiled

}
