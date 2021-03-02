STOCfree_Stan_data <- function(STOCfree_data){

  ## test data
  test_data <- STOCfree_data$test_data[, c("status_id", "herd_id", "test_res")]
  rf_data   <- STOCfree_data$risk_factor_data[, c("status_id", "herd_id", "month_id")]

  Stan_test_data <- merge(rf_data, test_data, all.x = TRUE)

  ## repalcing missing test data with the value 3
  Stan_test_data$test_res[is.na(Stan_test_data$test_res)] <- 3

  Stan_test_data <- Stan_test_data[order(Stan_test_data$status_id),]

  ## indices for herd level tests
  herds_t1 <- with(Stan_test_data,
                   tapply(status_id, herd_id, min))

  herds_T <- with(Stan_test_data,
                  tapply(status_id, herd_id, max))

  ## number of tests
  N <- nrow(Stan_test_data)

  ## number of herds
  n_herds <- attr(STOCfree_data, "number of herds")

  Stan_data <- list(
    N = N,
    n_herds = n_herds,
    herds_t1 = herds_t1,
    herds_t2 = herds_t1 + 1,
    herds_T = herds_T,
    test_res = Stan_test_data$test_res,
    Se_beta_a = STOCfree_data$test_perf_prior$Se_a,
    Se_beta_b = STOCfree_data$test_perf_prior$Se_b,
    Sp_beta_a = STOCfree_data$test_perf_prior$Sp_a,
    Sp_beta_b = STOCfree_data$test_perf_prior$Sp_b,
    pi1_beta_a = STOCfree_data$inf_dyn_priors["pi1_a"],
    pi1_beta_b = STOCfree_data$inf_dyn_priors["pi1_b"],
    tau1_beta_a = STOCfree_data$inf_dyn_priors["tau1_a"],
    tau1_beta_b = STOCfree_data$inf_dyn_priors["tau1_b"],
    tau2_beta_a = STOCfree_data$inf_dyn_priors["tau2_a"],
    tau2_beta_b = STOCfree_data$inf_dyn_priors["tau2_b"]
  )

  return(Stan_data)

}
