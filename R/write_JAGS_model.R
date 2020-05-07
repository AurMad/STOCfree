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

write_JAGS_model.herd <- function(data){

  cat('model{
##############################################################################
###  Inference from historical data
##############################################################################

  ## First test in a herd - status type = 1
  for(i1 in 1:n_status_typ1){

  ## prior probability of infection on first test
  pi1_init[i1] ~ dbeta(pi1_beta_a, pi1_beta_b)

  ## prior probability updated with test result
  pi[status_typ1[i1]] <- test_res_typ1[i1] * (
     Se[test_id_typ1[i1]] * pi1_init[i1] /
     (Se[test_id_typ1[i1]] * pi1_init[i1] +
      (1 - Sp[test_id_typ1[i1]]) * (1 - pi1_init[i1]))
    ) +
     (1 - test_res_typ1[i1]) * (
     (1 - Se[test_id_typ1[i1]]) * pi1_init[i1] /
     ((1 - Se[test_id_typ1[i1]]) * pi1_init[i1] +
       Sp[test_id_typ1[i1]] * (1 - pi1_init[i1])
     )
     )

  # status sampled from probability
  Status[status_typ1[i1]] ~ dbern(pi[status_typ1[i1]])

  }

  ## 2: first test on a month which is not first test in herd
  ## status type = 2
  for(i2 in 1:n_status_typ2){

  # probability of infection given previous status and dynamics
  pi2_init[i2] <- tau1 * (1 - Status[status_typ2[i2] - 1]) +
                  tau2 * Status[status_typ2[i2] - 1]

  # probability of infection updated with test result
  pi[status_typ2[i2]] <- test_res_typ2[i2] * (
     Se[test_id_typ2[i2]] * pi2_init[i2] /
     (Se[test_id_typ2[i2]] * pi2_init[i2] +
      (1 - Sp[test_id_typ2[i2]]) * (1 - pi2_init[i2]))
    ) +
     (1 - test_res_typ2[i2]) * (
     (1 - Se[test_id_typ2[i2]]) * pi2_init[i2] /
     ((1 - Se[test_id_typ2[i2]]) * pi2_init[i2] +
       Sp[test_id_typ2[i2]] * (1 - pi2_init[i2])
     )
     )

  # status sampled from probability
  Status[status_typ2[i2]] ~ dbern(pi[status_typ2[i2]])

  }

  ## 3: test > 1 on a month
  ## status type = 3
  for(i3 in 1:n_status_typ3){

  # probability of infection is given by previous test result on the same month
  pi3_init[i3] <- pi[status_typ3[i3] - 1]

  # probability of infection updated with test result
  pi[status_typ3[i3]] <- test_res_typ3[i3] * (
     Se[test_id_typ3[i3]] * pi3_init[i3] /
     (Se[test_id_typ3[i3]] * pi3_init[i3] +
      (1 - Sp[test_id_typ3[i3]]) * (1 - pi3_init[i3]))
    ) +
     (1 - test_res_typ3[i3]) * (
     (1 - Se[test_id_typ3[i3]]) * pi3_init[i3] /
     ((1 - Se[test_id_typ3[i3]]) * pi3_init[i3] +
       Sp[test_id_typ3[i3]] * (1 - pi3_init[i3])
     )
     )

  # status sampled from probability
  Status[status_typ3[i3]] ~ dbern(pi[status_typ3[i3]])

  }

  ## no test
  ## no_test
  for(j in 1:n_no_test){

  # probability of infection given previous status and dynamics
  pi[status_no_test[j]] <- tau1 * (1 - Status[status_no_test[j] - 1]) +
                           tau2 * Status[status_no_test[j] - 1]

  # status sampled from probability
  Status[status_no_test[j]] ~ dbern(pi[status_no_test[j]])

  }

##############################################################################
###  Prediction of probability of infection
##############################################################################
## 4: status to predict without test result
for(i4 in 1:n_status_typ4){

  # probability of infection given previous status and dynamics
  predicted_proba[herd_id_pr4[i4]] <- tau1 * (1 - Status[status_typ4[i4] - 1]) +
                 tau2 * Status[status_typ4[i4] - 1]

  # status sampled from probability
  predicted_status[herd_id_pr4[i4]] ~ dbern(predicted_proba[herd_id_pr4[i4]])

  }

## 5: status to predict with a single test performed
for(i5 in 1:n_status_typ5){

  # probability of infection given previous status and dynamics
  pi5_init[i5] <- tau1 * (1 - Status[status_typ5[i5] - 1]) +
                  tau2 * Status[status_typ5[i5] - 1]

  # probability of infection updated with test result
  predicted_proba[herd_id_pr5[i5]] <- test_res_typ5[i5] * (
     Se[test_id_typ5[i5]] * pi5_init[i5] /
     (Se[test_id_typ5[i5]] * pi5_init[i5] +
      (1 - Sp[test_id_typ5[i5]]) * (1 - pi5_init[i5]))
    ) +
     (1 - test_res_typ5[i5]) * (
     (1 - Se[test_id_typ5[i5]]) * pi5_init[i5] /
     ((1 - Se[test_id_typ5[i5]]) * pi5_init[i5] +
       Sp[test_id_typ5[i5]] * (1 - pi5_init[i5])
     )
     )

  # status sampled from probability
  predicted_status[herd_id_pr5[i5]] ~ dbern(predicted_proba[herd_id_pr5[i5]])

  }

## 6: status to predict with several tests on this month
## same as above except that pi_init is probability of infection after previous test
## no dynamics
for(i6 in 1:n_status_typ6){

  # probability of infection given previous status and dynamics
  pi6_init[i6] <- Status[status_typ6[i6] - 1]

  # probability of infection updated with test result
  predicted_proba[herd_id_pr6[i6]] <- test_res_typ6[i6] * (
     Se[test_id_typ6[i6]] * pi6_init[i6] /
     (Se[test_id_typ6[i6]] * pi6_init[i6] +
      (1 - Sp[test_id_typ6[i6]]) * (1 - pi6_init[i6]))
    ) +
     (1 - test_res_typ6[i6]) * (
     (1 - Se[test_id_typ6[i6]]) * pi6_init[i6] /
     ((1 - Se[test_id_typ6[i6]]) * pi6_init[i6] +
       Sp[test_id_typ6[i6]] * (1 - pi6_init[i6])
     )
     )

  # status sampled from probability
  predicted_status[herd_id_pr6[i6]] ~ dbern(predicted_proba[herd_id_pr6[i6]])

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
write_JAGS_model.herd_rf <- function(data){

  test_data <- data$test_data
  ## several tests on the same month
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])


  cat('model{
##############################################################################
###  Inference from historical data
##############################################################################

  ## First test in a herd - status type = 1
  for(i1 in 1:n_status_typ1){

  ## prior probability of infection on first test
  pi1_init[i1] ~ dbeta(pi1_beta_a, pi1_beta_b)

  ## prior probability updated with test result
  pi[status_typ1[i1]] <- test_res_typ1[i1] * (
     Se[test_id_typ1[i1]] * pi1_init[i1] /
     (Se[test_id_typ1[i1]] * pi1_init[i1] +
      (1 - Sp[test_id_typ1[i1]]) * (1 - pi1_init[i1]))
    ) +
     (1 - test_res_typ1[i1]) * (
     (1 - Se[test_id_typ1[i1]]) * pi1_init[i1] /
     ((1 - Se[test_id_typ1[i1]]) * pi1_init[i1] +
       Sp[test_id_typ1[i1]] * (1 - pi1_init[i1])
     )
     )

  # status sampled from probability
  Status[status_typ1[i1]] ~ dbern(pi[status_typ1[i1]])

  }

  ## 2: first test on a month which is not first test in herd
  ## status type = 2
  for(i2 in 1:n_status_typ2){

  # probability of new infection
  logit(tau1[status_typ2[i2]]) <- inprod(risk_factors[status_typ2[i2],], theta)

  # probability of infection given previous status and dynamics
  pi2_init[i2] <- tau1[status_typ2[i2]] * (1 - Status[status_typ2[i2] - 1]) +
                  tau2 * Status[status_typ2[i2] - 1]

  # probability of infection updated with test result
  pi[status_typ2[i2]] <- test_res_typ2[i2] * (
     Se[test_id_typ2[i2]] * pi2_init[i2] /
     (Se[test_id_typ2[i2]] * pi2_init[i2] +
      (1 - Sp[test_id_typ2[i2]]) * (1 - pi2_init[i2]))
    ) +
     (1 - test_res_typ2[i2]) * (
     (1 - Se[test_id_typ2[i2]]) * pi2_init[i2] /
     ((1 - Se[test_id_typ2[i2]]) * pi2_init[i2] +
       Sp[test_id_typ2[i2]] * (1 - pi2_init[i2])
     )
     )

  # status sampled from probability
  Status[status_typ2[i2]] ~ dbern(pi[status_typ2[i2]])

  }',

  if(n_status_typ3 == 0){''} else {'
  ## 3: test > 1 on a month
  ## status type = 3
  for(i3 in 1:n_status_typ3){

  # probability of infection is given by previous test result on the same month
  pi3_init[i3] <- pi[status_typ3[i3] - 1]

  # probability of infection updated with test result
  pi[status_typ3[i3]] <- test_res_typ3[i3] * (
     Se[test_id_typ3[i3]] * pi3_init[i3] /
     (Se[test_id_typ3[i3]] * pi3_init[i3] +
      (1 - Sp[test_id_typ3[i3]]) * (1 - pi3_init[i3]))
    ) +
     (1 - test_res_typ3[i3]) * (
     (1 - Se[test_id_typ3[i3]]) * pi3_init[i3] /
     ((1 - Se[test_id_typ3[i3]]) * pi3_init[i3] +
       Sp[test_id_typ3[i3]] * (1 - pi3_init[i3])
     )
     )

  # status sampled from probability
  Status[status_typ3[i3]] ~ dbern(pi[status_typ3[i3]])

  }'},
  '## no test
  ## no_test
  for(j in 1:n_no_test){

  # probability of new infection
  logit(tau1[status_no_test[j]]) <- inprod(risk_factors[status_no_test[j],], theta)

  # probability of infection given previous status and dynamics
  pi[status_no_test[j]] <- tau1[status_no_test[j]] * (1 - Status[status_no_test[j] - 1]) +
                           tau2 * Status[status_no_test[j] - 1]

  # status sampled from probability
  Status[status_no_test[j]] ~ dbern(pi[status_no_test[j]])

  }

##############################################################################
###  Prediction of probability of infection
##############################################################################
## 4: status to predict without test result
for(i4 in 1:n_status_typ4){

  # probability of new infection
  logit(tau1[status_typ4[i4]]) <- inprod(risk_factors[status_typ4[i4],], theta)

  # probability of infection given previous status and dynamics
  predicted_proba[herd_id_pr4[i4]] <- tau1[status_typ4[i4]] * (1 - Status[status_typ4[i4] - 1]) +
                 tau2 * Status[status_typ4[i4] - 1]

  # status sampled from probability
  predicted_status[herd_id_pr4[i4]] ~ dbern(predicted_proba[herd_id_pr4[i4]])

  }

## 5: status to predict with a single test performed
for(i5 in 1:n_status_typ5){

  # probability of new infection
  logit(tau1[status_typ5[i5]]) <- inprod(risk_factors[status_typ5[i5],], theta)

  # probability of infection given previous status and dynamics
  pi5_init[i5] <- tau1[status_typ5[i5]] * (1 - Status[status_typ5[i5] - 1]) + tau2 * Status[status_typ5[i5] - 1]

  # probability of infection updated with test result
  predicted_proba[herd_id_pr5[i5]] <- test_res_typ5[i5] * (
     Se[test_id_typ5[i5]] * pi5_init[i5] /
     (Se[test_id_typ5[i5]] * pi5_init[i5] +
      (1 - Sp[test_id_typ5[i5]]) * (1 - pi5_init[i5]))
    ) +
     (1 - test_res_typ5[i5]) * (
     (1 - Se[test_id_typ5[i5]]) * pi5_init[i5] /
     ((1 - Se[test_id_typ5[i5]]) * pi5_init[i5] +
       Sp[test_id_typ5[i5]] * (1 - pi5_init[i5])
     )
     )

  # status sampled from probability
  predicted_status[herd_id_pr5[i5]] ~ dbern(predicted_proba[herd_id_pr5[i5]])

  }',
  if(n_status_typ6 == 0){''} else {'
## 6: status to predict with several tests on this month
## same as above except that pi_init is probability of infection after previous test
## no dynamics
for(i6 in 1:n_status_typ6){

  # probability of infection given previous status and dynamics
  pi6_init[i6] <- Status[status_typ6[i6] - 1]

  # probability of infection updated with test result
  predicted_proba[herd_id_pr6[i6]] <- test_res_typ6[i6] * (
     Se[test_id_typ6[i6]] * pi6_init[i6] /
     (Se[test_id_typ6[i6]] * pi6_init[i6] +
      (1 - Sp[test_id_typ6[i6]]) * (1 - pi6_init[i6]))
    ) +
     (1 - test_res_typ6[i6]) * (
     (1 - Se[test_id_typ6[i6]]) * pi6_init[i6] /
     ((1 - Se[test_id_typ6[i6]]) * pi6_init[i6] +
       Sp[test_id_typ6[i6]] * (1 - pi6_init[i6])
     )
     )

  # status sampled from probability
  predicted_status[herd_id_pr6[i6]] ~ dbern(predicted_proba[herd_id_pr6[i6]])

  }'},
  '
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
write_JAGS_model.animal <- function(data){

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

write_JAGS_model.animal_rf <- function(data){

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

   p_test_pos[i] <- Se[test_id[i]] * p_infected[i] +
                  (1 - Sp[test_id[i]]) * (1 - p_infected[i])

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
      (Se[test_id_for_pred[j]] * pi[ind_p_test[j]] * pi_within) /
      (Se[test_id_for_pred[j]] * pi[ind_p_test[j]] * pi_within + (1 - Sp[test_id_for_pred[j]]) * (1 - pi[ind_p_test[j]])) +
    (1 - p_pos[j]) *
      (1 - Se[test_id_for_pred[j]]) * pi[ind_p_test[j]]  * pi_within /
      ((1 - Se[test_id_for_pred[j]]) * pi[ind_p_test[j]] * pi_within + Sp[test_id_for_pred[j]] * (1 - pi[ind_p_test[j]]))

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
  for(i_test in 1:n_tests){

  Se[i_test] ~ dbeta(Se_beta_a[i_test], Se_beta_b[i_test])
  Sp[i_test] ~ dbeta(Sp_beta_a[i_test], Sp_beta_b[i_test])

  }

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
