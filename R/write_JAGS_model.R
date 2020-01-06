#' Writes a text file with the JAGS model. For internal use.
#'
#' This function is for internal use.
#'
#' @return a text file with the JAGS model
#' @export
#'
#' @examples
write_JAGS_model <- function(){

cat('model{

  ## Loop for infection dynamics up to test before last
  for(h in 1:n_herds){

    # First month
    pi[ind_i[h]] ~ dbeta(pi1_beta_a, pi1_beta_b)
    Status[ind_i[h]] ~ dbern(pi[ind_i[h]])

    # Months 2 to 1 minus final
    for(t in ind_j[h]:ind_f[h]){

      logit(tau1[t]) <- inprod(risk_factors[t,], theta[, 1])

      pi[t] <- tau1[t] * (1 - Status[t-1]) +
        tau2 * Status[t-1]

      Status[t] ~ dbern(pi[t])

    } #t

    ## tau1 for test to predict
    logit(tau1[ind_p[h]]) <- inprod(risk_factors[ind_p[h],], theta[, 1])

  }


  ## Loop for statuses to predict when test result is available
  for(j in 1:n_pred_test){

    # After having observed the risk factors and the test result
    pi[ind_p[ind_last_is_test[j]]] <- tau1[ind_p[ind_last_is_test[j]]] * (1 - Status[ind_f[ind_last_is_test[j]]]) +
      tau2 * Status[ind_f[ind_last_is_test[j]]]

    predicted_proba[ind_last_is_test[j]] <- test_res[ind_p[ind_last_is_test[j]]] *
      (Se * pi[ind_p[ind_last_is_test[j]]]) /
      (Se * pi[ind_p[ind_last_is_test[j]]] + (1 - Sp) * (1 - pi[ind_p[ind_last_is_test[j]]])) +
      (1 - test_res[ind_p[ind_last_is_test[j]]]) *
      (1 - Se) * pi[ind_p[ind_last_is_test[j]]] /
      ((1 - Se) * pi[ind_p[ind_last_is_test[j]]] + Sp * (1 - pi[ind_p[ind_last_is_test[j]]]))

    predicted_status[ind_last_is_test[j]] ~ dbern(predicted_proba[ind_last_is_test[j]])

  }


  ## Loop for statuses to predict when test result is not available
  for(k in 1:n_pred_no_test){

    predicted_proba[ind_last_is_not_test[k]] <- tau1[ind_p[ind_last_is_not_test[k]]] * (1 - Status[ind_f[ind_last_is_not_test[k]]]) +
      tau2 * Status[ind_f[ind_last_is_not_test[k]]]

    predicted_status[ind_last_is_not_test[k]] ~ dbern(predicted_proba[ind_last_is_not_test[k]])

  }

  ## Loop for test results
  for(i in 1:n_tests){

    pTestPos[i] <-  Se * Status[ind_test[i]] +
      (1 - Sp) * (1 - Status[ind_test[i]])

    test_res[ind_test[i]] ~ dbern(pTestPos[i])

  }


  # Priors

  ## Priors for sensitivities and specificities
  Se ~ dbeta(Se_beta_a, Se_beta_b)
  Sp ~ dbeta(Sp_beta_a, Sp_beta_b)

  ## Probability of not eliminating the infection
  tau2 ~ dbeta(tau2_beta_a, tau2_beta_b)

  ## Logistic regression coefficients
  for(i_rf in 1:n_risk_factors){

    theta[i_rf, 1] ~ dnorm(theta_norm_mean[i_rf], theta_norm_prec[i_rf])

  }

}',
    file =  "JAGS_model.txt")
}
