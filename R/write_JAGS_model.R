#' Writes a text file with the JAGS model. For internal use.
#'
#' This function is for internal use.
#'
#' @return a text file with the JAGS model
#' @export
#'
#' @examples
write_JAGS_model <- function(data){

  UseMethod("write_JAGS_model")

  }


#' @export
write_JAGS_model.default <- function(data){

  print("No method defined for this type of data")

}


#' @export
write_JAGS_model.herd_1test <- function(data){

  cat('model{

##############################################################################
###  Inference from historical data
##############################################################################

  ## Loop for infection dynamics up to test before last
  for(h in 1:n_herds){

    # First month
    pi[ind_i[h]] ~ dbeta(pi1_beta_a, pi1_beta_b)
    Status[ind_i[h]] ~ dbern(pi[ind_i[h]])

    # Months 2 to 1 minus final
    for(t in ind_j[h]:ind_f[h]){

      pi[t] <- tau1 * (1 - Status[t-1]) +
               tau2 * Status[t-1]

      Status[t] ~ dbern(pi[t])

    } #t


  }


  ## Loop for test results - historical data
  for(i in 1:n_tests_perf){

    pTestPos[i] <-  Se * Status[ind_test[i]] +
                    (1 - Sp) * (1 - Status[ind_test[i]])

    test_res[i] ~ dbern(pTestPos[i])

  }


##############################################################################
###  Prediction of probability of infection
##############################################################################

  ## Loop for statuses to predict when test result is available
  for(j in 1:n_pred_test){

    ## Indices in loops:
    # ind_p_test: relates to status index
    # ind_last_is_test relates to herd_id -> 1 to number of herds between this loop and the next
    # test_for_pred relates to the cases when test is available for prediction -> 1 to number of such cases

    # After having observed the risk factors and the test result
    pi[ind_p_test[j]] <- tau1 * (1 - Status[ind_p_test[j] - 1]) +
                         tau2 * Status[ind_p_test[j] - 1]

    predicted_proba[ind_last_is_test[j]] <-
    test_for_pred[j] *
      (Se * pi[ind_p_test[j]]) /
      (Se * pi[ind_p_test[j]] + (1 - Sp) * (1 - pi[ind_p_test[j]])) +
      (1 - test_for_pred[j]) *
      (1 - Se) * pi[ind_p_test[j]] /
      ((1 - Se) * pi[ind_p_test[j]] + Sp * (1 - pi[ind_p_test[j]]))

    predicted_status[ind_last_is_test[j]] ~ dbern(predicted_proba[ind_last_is_test[j]])

  }


  ## Loop for statuses to predict when test result is not available
  for(k in 1:n_pred_no_test){

    predicted_proba[ind_last_is_not_test[k]] <- tau1 * (1 - Status[ind_p_no_test[k] - 1]) +
                                                tau2 * Status[ind_f[ind_last_is_not_test[k]]]

    predicted_status[ind_last_is_not_test[k]] ~ dbern(predicted_proba[ind_last_is_not_test[k]])

  }


##############################################################################
###  Priors
##############################################################################

  ## Priors for sensitivities and specificities
  Se ~ dbeta(Se_beta_a, Se_beta_b)
  Sp ~ dbeta(Sp_beta_a, Sp_beta_b)

  ## Probability of not eliminating the infection
  tau1 ~ dbeta(tau1_beta_a, tau1_beta_b)
  tau2 ~ dbeta(tau2_beta_a, tau2_beta_b)

  ## pi1 <- tau1 / (1 + tau1 - tau2)


}',
      file =  "JAGS_model.txt")
}


#' @export
write_JAGS_model.herd_1test_rf <- function(data){

  cat('model{

##############################################################################
###  Inference from historical data
##############################################################################

  ## Loop for infection dynamics up to test before last
  for(h in 1:n_herds){

    # First month
    pi[ind_i[h]] ~ dbeta(pi1_beta_a, pi1_beta_b)
    Status[ind_i[h]] ~ dbern(pi[ind_i[h]])

    # Months 2 to 1 minus final
    for(t in ind_j[h]:ind_f[h]){

      logit(tau1[t]) <- inprod(risk_factors[t,], theta)

      pi[t] <- tau1[t] * (1 - Status[t-1]) +
               tau2 * Status[t-1]

      Status[t] ~ dbern(pi[t])

    } #t

    ## tau1 for test to predict
    logit(tau1[ind_p[h]]) <- inprod(risk_factors[ind_p[h],], theta)

  }


  ## Loop for test results - historical data
  for(i in 1:n_tests_perf){

    pTestPos[i] <-  Se * Status[ind_test[i]] +
                    (1 - Sp) * (1 - Status[ind_test[i]])

    test_res[i] ~ dbern(pTestPos[i])

  }


##############################################################################
###  Prediction of probability of infection
##############################################################################

  ## Loop for statuses to predict when test result is available
  for(j in 1:n_pred_test){

    ## Indices in loops:
    # ind_p_test: relates to status index
    # ind_last_is_test relates to herd_id -> 1 to number of herds between this loop and the next
    # test_for_pred relates to the cases when test is available for prediction -> 1 to number of such cases

    # After having observed the risk factors and the test result
    pi[ind_p_test[j]] <- tau1[ind_p_test[j]] * (1 - Status[ind_p_test[j] - 1]) +
                               tau2 * Status[ind_p_test[j] - 1]

    predicted_proba[ind_last_is_test[j]] <-
    test_for_pred[j] *
      (Se * pi[ind_p_test[j]]) /
      (Se * pi[ind_p_test[j]] + (1 - Sp) * (1 - pi[ind_p_test[j]])) +
      (1 - test_for_pred[j]) *
      (1 - Se) * pi[ind_p_test[j]] /
      ((1 - Se) * pi[ind_p_test[j]] + Sp * (1 - pi[ind_p_test[j]]))

    predicted_status[ind_last_is_test[j]] ~ dbern(predicted_proba[ind_last_is_test[j]])

  }


  ## Loop for statuses to predict when test result is not available
  for(k in 1:n_pred_no_test){

    predicted_proba[ind_last_is_not_test[k]] <- tau1[ind_p_no_test[k]] * (1 - Status[ind_p_no_test[k] - 1]) +
                                                tau2 * Status[ind_f[ind_last_is_not_test[k]]]

    predicted_status[ind_last_is_not_test[k]] ~ dbern(predicted_proba[ind_last_is_not_test[k]])

  }


##############################################################################
###  Priors
##############################################################################

  ## Priors for sensitivities and specificities
  Se ~ dbeta(Se_beta_a, Se_beta_b)
  Sp ~ dbeta(Sp_beta_a, Sp_beta_b)

  ## Probability of not eliminating the infection
  tau2 ~ dbeta(tau2_beta_a, tau2_beta_b)

  ## Logistic regression coefficients
  for(i_rf in 1:n_risk_factors){

    theta[i_rf] ~ dnorm(theta_norm_mean[i_rf], theta_norm_prec[i_rf])

  }

}',
      file =  "JAGS_model.txt")
}


#' @export
write_JAGS_model.herd_ntests <- function(data){

  cat('model{

##############################################################################
###  Inference from historical data
##############################################################################

  ## Loop for infection dynamics up to test before last
  for(h in 1:n_herds){

    # First month
    pi[ind_i[h]] ~ dbeta(pi1_beta_a, pi1_beta_b)
    Status[ind_i[h]] ~ dbern(pi[ind_i[h]])

    # Months 2 to 1 minus final
    for(t in ind_j[h]:ind_f[h]){

      pi[t] <- tau1 * (1 - Status[t-1]) +
               tau2 * Status[t-1]

      Status[t] ~ dbern(pi[t])

    } #t


  }


  ## Loop for test results - historical data
  for(i in 1:n_tests_perf){

    pTestPos[i] <-  Se[test_id[i]] * Status[ind_test[i]] +
                    (1 - Sp[test_id[i]]) * (1 - Status[ind_test[i]])

    test_res[i] ~ dbern(pTestPos[i])

  }


##############################################################################
###  Prediction of probability of infection
##############################################################################

  ## Loop for statuses to predict when test result is available
  for(j in 1:n_pred_test){

    ## Indices in loops:
    # ind_p_test: relates to status index
    # ind_last_is_test relates to herd_id -> 1 to number of herds between this loop and the next
    # test_for_pred relates to the cases when test is available for prediction -> 1 to number of such cases

    # After having observed the risk factors and the test result
    pi[ind_p_test[j]] <- tau1 * (1 - Status[ind_p_test[j] - 1]) +
                         tau2 * Status[ind_p_test[j] - 1]

    predicted_proba[ind_last_is_test[j]] <-
    test_for_pred[j] *
      (Se[test_id_for_pred[j]] * pi[ind_p_test[j]]) /
      (Se[test_id_for_pred[j]] * pi[ind_p_test[j]] + (1 - Sp[test_id_for_pred[j]]) * (1 - pi[ind_p_test[j]])) +
      (1 - test_for_pred[j]) *
      (1 - Se[test_id_for_pred[j]]) * pi[ind_p_test[j]] /
      ((1 - Se[test_id_for_pred[j]]) * pi[ind_p_test[j]] + Sp[test_id_for_pred[j]] * (1 - pi[ind_p_test[j]]))

    predicted_status[ind_last_is_test[j]] ~ dbern(predicted_proba[ind_last_is_test[j]])

  }


  ## Loop for statuses to predict when test result is not available
  for(k in 1:n_pred_no_test){

    predicted_proba[ind_last_is_not_test[k]] <- tau1 * (1 - Status[ind_p_no_test[k] - 1]) +
                                                tau2 * Status[ind_f[ind_last_is_not_test[k]]]

    predicted_status[ind_last_is_not_test[k]] ~ dbern(predicted_proba[ind_last_is_not_test[k]])

  }


##############################################################################
###  Priors
##############################################################################

  ## Priors for sensitivities and specificities
  for(i_test in 1:n_tests){

  Se[i_test] ~ dbeta(Se_beta_a[i_test], Se_beta_b[i_test])
  Sp[i_test] ~ dbeta(Sp_beta_a[i_test], Sp_beta_b[i_test])

  }

  ## Probability of not eliminating the infection
  tau1 ~ dbeta(tau1_beta_a, tau1_beta_b)
  tau2 ~ dbeta(tau2_beta_a, tau2_beta_b)

}',
      file =  "JAGS_model.txt")
}




#' @export
write_JAGS_model.herd_ntests_rf <- function(data){

  cat('model{

##############################################################################
###  Inference from historical data
##############################################################################

  ## Loop for infection dynamics up to test before last
  for(h in 1:n_herds){

    # First month
    pi[ind_i[h]] ~ dbeta(pi1_beta_a, pi1_beta_b)
    Status[ind_i[h]] ~ dbern(pi[ind_i[h]])

    # Months 2 to 1 minus final
    for(t in ind_j[h]:ind_f[h]){

      logit(tau1[t]) <- inprod(risk_factors[t,], theta)

      pi[t] <- tau1[t] * (1 - Status[t-1]) +
               tau2 * Status[t-1]

      Status[t] ~ dbern(pi[t])

    } #t

    ## tau1 for test to predict
    logit(tau1[ind_p[h]]) <- inprod(risk_factors[ind_p[h],], theta)

  }


  ## Loop for test results - historical data
  for(i in 1:n_tests_perf){

    pTestPos[i] <-  Se[test_id[i]] * Status[ind_test[i]] +
                    (1 - Sp[test_id[i]]) * (1 - Status[ind_test[i]])

    test_res[i] ~ dbern(pTestPos[i])

  }


##############################################################################
###  Prediction of probability of infection
##############################################################################

  ## Loop for statuses to predict when test result is available
  for(j in 1:n_pred_test){

    ## Indices in loops:
    # ind_p_test: relates to status index
    # ind_last_is_test relates to herd_id -> 1 to number of herds between this loop and the next
    # test_for_pred relates to the cases when test is available for prediction -> 1 to number of such cases

    # After having observed the risk factors and the test result
    pi[ind_p_test[j]] <- tau1[ind_p_test[j]] * (1 - Status[ind_p_test[j] - 1]) +
                               tau2 * Status[ind_p_test[j] - 1]

    predicted_proba[ind_last_is_test[j]] <-
    test_for_pred[j] *
      (Se[test_id_for_pred[j]] * pi[ind_p_test[j]]) /
      (Se[test_id_for_pred[j]] * pi[ind_p_test[j]] + (1 - Sp[test_id_for_pred[j]]) * (1 - pi[ind_p_test[j]])) +
      (1 - test_for_pred[j]) *
      (1 - Se[test_id_for_pred[j]]) * pi[ind_p_test[j]] /
      ((1 - Se[test_id_for_pred[j]]) * pi[ind_p_test[j]] + Sp[test_id_for_pred[j]] * (1 - pi[ind_p_test[j]]))

    predicted_status[ind_last_is_test[j]] ~ dbern(predicted_proba[ind_last_is_test[j]])

  }


  ## Loop for statuses to predict when test result is not available
  for(k in 1:n_pred_no_test){

    predicted_proba[ind_last_is_not_test[k]] <- tau1[ind_p_no_test[k]] * (1 - Status[ind_p_no_test[k] - 1]) +
                                                tau2 * Status[ind_f[ind_last_is_not_test[k]]]

    predicted_status[ind_last_is_not_test[k]] ~ dbern(predicted_proba[ind_last_is_not_test[k]])

  }


##############################################################################
###  Priors
##############################################################################

  ## Priors for sensitivities and specificities
  for(i_test in 1:n_tests){

  Se[i_test] ~ dbeta(Se_beta_a[i_test], Se_beta_b[i_test])
  Sp[i_test] ~ dbeta(Sp_beta_a[i_test], Sp_beta_b[i_test])

  }

  ## Probability of not eliminating the infection
  tau2 ~ dbeta(tau2_beta_a, tau2_beta_b)

  ## Logistic regression coefficients
  for(i_rf in 1:n_risk_factors){

    theta[i_rf] ~ dnorm(theta_norm_mean[i_rf], theta_norm_prec[i_rf])

  }

}',
      file =  "JAGS_model.txt")
}


#' @export
write_JAGS_model.animal_1test <- function(data){

cat('model{

##############################################################################
###  Inference from historical data
##############################################################################

  ## Loop for infection dynamics up to test before last
  for(h in 1:n_herds){

    # First month
    pi[ind_i[h]] ~ dbeta(pi1_beta_a, pi1_beta_b)
    Status[ind_i[h]] ~ dbern(pi[ind_i[h]])

    # Months 2 to 1 minus final
    for(t in ind_j[h]:ind_f[h]){

      pi[t] <- tau1 * (1 - Status[t-1]) +
               tau2 * Status[t-1]

      Status[t] ~ dbern(pi[t])

    } #t


  }

  ## Loop for test results - historical data
  for(i in 1:n_tests_perf){

   p_infected[i] <- Status[ind_test[i]] * pi_within

   p_test_pos[i] <- Se * p_infected[i] +
                  (1 - Sp) * (1 - p_infected[i])

   n_pos[i] ~ dbin(p_test_pos[i], n_tested[i])

   }


##############################################################################
###  Prediction of probability of infection
##############################################################################

  ## Loop for statuses to predict when test result is available
  for(j in 1:n_pred_test){

    ## Indices in loops:
    # ind_p_test: relates to status index
    # ind_last_is_test relates to herd_id -> 1 to number of herds between this loop and the next
    # test_for_pred relates to the cases when test is available for prediction -> 1 to number of such cases

    ## Distribution of number of positive tests
    n_pos_for_pred[j] ~ dbinom(p_pos[j], n_tested_for_pred[j])
    p_pos[j] ~ dbeta(1, 1)

    # After having observed the risk factors and the test result
    pi[ind_p_test[j]] <- tau1 * (1 - Status[ind_p_test[j] - 1]) +
                         tau2 * Status[ind_p_test[j] - 1]

    predicted_proba[ind_last_is_test[j]] <-
    p_pos[j] *
      (Se * pi[ind_p_test[j]] * pi_within) /
      (Se * pi[ind_p_test[j]] * pi_within + (1 - Sp) * (1 - pi[ind_p_test[j]])) +
    (1 - p_pos[j]) *
      (1 - Se) * pi[ind_p_test[j]]  * pi_within /
      ((1 - Se) * pi[ind_p_test[j]] * pi_within + Sp * (1 - pi[ind_p_test[j]]))

    predicted_status[ind_last_is_test[j]] ~ dbern(predicted_proba[ind_last_is_test[j]])

  }


  ## Loop for statuses to predict when test result is not available
  for(k in 1:n_pred_no_test){

    predicted_proba[ind_last_is_not_test[k]] <- tau1 * (1 - Status[ind_p_no_test[k] - 1]) +
                                                tau2 * Status[ind_f[ind_last_is_not_test[k]]]

    predicted_status[ind_last_is_not_test[k]] ~ dbern(predicted_proba[ind_last_is_not_test[k]])

  }


##############################################################################
###  Priors
##############################################################################

  ## Priors for sensitivities and specificities
  Se ~ dbeta(Se_beta_a, Se_beta_b)
  Sp ~ dbeta(Sp_beta_a, Sp_beta_b)

  ## Prior for within herd prevalence in infected herds
  pi_within ~ dbeta(pi_within_a, pi_within_b)

  ## Probability of not eliminating the infection
  tau1 ~ dbeta(tau1_beta_a, tau1_beta_b)
  tau2 ~ dbeta(tau2_beta_a, tau2_beta_b)


}',
      file =  "JAGS_model.txt")
}



write_JAGS_model.animal_1test_rf <- function(data){

  cat('model{

##############################################################################
###  Inference from historical data
##############################################################################

  ## Loop for infection dynamics up to test before last
  for(h in 1:n_herds){

    # First month
    pi[ind_i[h]] ~ dbeta(pi1_beta_a, pi1_beta_b)
    Status[ind_i[h]] ~ dbern(pi[ind_i[h]])

    # Months 2 to 1 minus final
    for(t in ind_j[h]:ind_f[h]){

    logit(tau1[t]) <- inprod(risk_factors[t,], theta)

    pi[t] <- tau1[t] * (1 - Status[t-1]) +
               tau2 * Status[t-1]

    Status[t] ~ dbern(pi[t])

    } #t

  ## tau1 for test to predict
  logit(tau1[ind_p[h]]) <- inprod(risk_factors[ind_p[h],], theta)

  }

  ## Loop for test results - historical data
  for(i in 1:n_tests_perf){

   p_infected[i] <- Status[ind_test[i]] * pi_within

   p_test_pos[i] <- Se * p_infected[i] +
                  (1 - Sp) * (1 - p_infected[i])

   n_pos[i] ~ dbin(p_test_pos[i], n_tested[i])

   }


##############################################################################
###  Prediction of probability of infection
##############################################################################

  ## Loop for statuses to predict when test result is available
  for(j in 1:n_pred_test){

    ## Indices in loops:
    # ind_p_test: relates to status index
    # ind_last_is_test relates to herd_id -> 1 to number of herds between this loop and the next
    # test_for_pred relates to the cases when test is available for prediction -> 1 to number of such cases

    ## Distribution of number of positive tests
    n_pos_for_pred[j] ~ dbinom(p_pos[j], n_tested_for_pred[j])
    p_pos[j] ~ dbeta(1, 1)

    # After having observed the risk factors and the test result
    pi[ind_p_test[j]] <- tau1[ind_p_test[j]] * (1 - Status[ind_p_test[j] - 1]) +
                         tau2 * Status[ind_p_test[j] - 1]

    predicted_proba[ind_last_is_test[j]] <-
    p_pos[j] *
      (Se * pi[ind_p_test[j]] * pi_within) /
      (Se * pi[ind_p_test[j]] * pi_within + (1 - Sp) * (1 - pi[ind_p_test[j]])) +
    (1 - p_pos[j]) *
      (1 - Se) * pi[ind_p_test[j]]  * pi_within /
      ((1 - Se) * pi[ind_p_test[j]] * pi_within + Sp * (1 - pi[ind_p_test[j]]))

    predicted_status[ind_last_is_test[j]] ~ dbern(predicted_proba[ind_last_is_test[j]])

  }


  ## Loop for statuses to predict when test result is not available
  for(k in 1:n_pred_no_test){

    predicted_proba[ind_last_is_not_test[k]] <- tau1[ind_last_is_not_test[k]] * (1 - Status[ind_p_no_test[k] - 1]) +
                                                tau2 * Status[ind_f[ind_last_is_not_test[k]]]

    predicted_status[ind_last_is_not_test[k]] ~ dbern(predicted_proba[ind_last_is_not_test[k]])

  }


##############################################################################
###  Priors
##############################################################################

  ## Priors for sensitivities and specificities
  Se ~ dbeta(Se_beta_a, Se_beta_b)
  Sp ~ dbeta(Sp_beta_a, Sp_beta_b)

  ## Prior for within herd prevalence in infected herds
  pi_within ~ dbeta(pi_within_a, pi_within_b)

  ## Probability of not eliminating the infection
  tau2 ~ dbeta(tau2_beta_a, tau2_beta_b)

  ## Logistic regression coefficients
  for(i_rf in 1:n_risk_factors){

  theta[i_rf] ~ dnorm(theta_norm_mean[i_rf], theta_norm_prec[i_rf])

 }

}',
      file =  "JAGS_model.txt")
}
