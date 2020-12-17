#' Writes a text file with the JAGS model. For internal use.
#'
#' This function is for internal use.
#'
#' @return a text file with the JAGS model
#'
write_JAGS_model <- function(data){

  UseMethod("write_JAGS_model")

  }


#' @export
write_JAGS_model.default <- function(data){

  print("No method defined for this type of data")

}

## JAGS model: herd level, no risk factor
write_JAGS_model.herd <- function(data){

  test_data <- data$test_data
  ## several tests on the same month
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ4 <- nrow(test_data[test_data$status_type == 4,])
  n_status_typ5 <- nrow(test_data[test_data$status_type == 5,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])

model <- paste('model{
##############################################################################
###  Inference from historical data
##############################################################################

  ## First test in a herd - status type = 1
  for(i1 in 1:n_status_typ1){

  ## prior probability of infection on first test
  pi[status_typ1[i1]] ~ dbeta(pi1_beta_a, pi1_beta_b)

  ## status sampled from probability - not sampled in this version
  Status[status_typ1[i1]] ~ dbern(pi[status_typ1[i1]])

  ## test result
  p_test_pos_typ1[i1] <- Se[test_id_typ1[i1]] * Status[status_typ1[i1]] +
                          (1 - Sp[test_id_typ1[i1]]) * (1 - Status[status_typ1[i1]])

  test_res_typ1[i1] ~ dbern(p_test_pos_typ1[i1])

  }

  ## 2: first test on a month which is not first test in herd
  ## status type = 2
  for(i2 in 1:n_status_typ2){

  # probability of infection given previous status and dynamics
  pi[status_typ2[i2]] <- tau1 * (1 - Status[status_typ2[i2] - 1]) +
                         tau2 * Status[status_typ2[i2] - 1]

  ## status sampled from probability
  Status[status_typ2[i2]] ~ dbern(pi[status_typ2[i2]])

  ## test result
  p_test_pos_typ2[i2] <- Se[test_id_typ2[i2]] * Status[status_typ2[i2]] +
                          (1 - Sp[test_id_typ2[i2]]) * (1 - Status[status_typ2[i2]])

  test_res_typ2[i2] ~ dbern(p_test_pos_typ2[i2])

  }',

      if(n_status_typ3 == 0){''} else {'
  ## 3: test > 1 on a month
  ## status type = 3
  for(i3 in 1:n_status_typ3){

  # probability of infection is given by previous test result on the same month
  pi[status_typ3[i3]] <- pi[status_typ3[i3] - 1]

  ## status sampled from probability
  Status[status_typ3[i3]] ~ dbern(pi[status_typ3[i3]])

  ## test result
  p_test_pos_typ3[i3] <- Se[test_id_typ3[i3]] * Status[status_typ3[i3]] +
                          (1 - Sp[test_id_typ3[i3]]) * (1 - Status[status_typ3[i3]])

  test_res_typ3[i3] ~ dbern(p_test_pos_typ3[i3])

  }
'},
  '## no test
  ## no_test
  for(j in 1:n_no_test){

  # probability of infection given previous status and dynamics
  pi[status_no_test[j]] <- tau1 * (1 - Status[status_no_test[j] - 1]) +
                           tau2 * Status[status_no_test[j] - 1]

  # # status sampled from probability
  Status[status_no_test[j]] ~ dbern(pi[status_no_test[j]])

  }

##############################################################################
###  Prediction of probability of infection
##############################################################################

',

if(n_status_typ4 == 0){''} else {'
## 4: status to predict without test result
for(i4 in 1:n_status_typ4){

  # probability of infection given previous status and dynamics
  predicted_proba[herd_id_pr4[i4]] <- tau1 * (1 - Status[status_typ4[i4] - 1]) +
                 tau2 * Status[status_typ4[i4] - 1]

  # status sampled from probability
  predicted_status[herd_id_pr4[i4]] ~ dbern(predicted_proba[herd_id_pr4[i4]])

  }'},

if(n_status_typ5 == 0){''} else {'

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

  }'},

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
###  Loop for monthly prevalences
##############################################################################
for(i_month in 1:month_max){

 month_prev[i_month] <- sum(Status[month_mat[1:month_N[i_month], i_month]]) /
                            month_N[i_month]

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

}
')

 }


## JAGS model: herd level, no risk factor
## Dynamics parameters modelled on the logit scale
write_JAGS_model.herd_dynLogit <- function(data){

  test_data <- data$test_data
  ## several tests on the same month
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ4 <- nrow(test_data[test_data$status_type == 4,])
  n_status_typ5 <- nrow(test_data[test_data$status_type == 5,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])

  model <- paste('model{
##############################################################################
###  Inference from historical data
##############################################################################

  ## First test in a herd - status type = 1
  for(i1 in 1:n_status_typ1){

  ## prior probability of infection on first test
  logit_pi[status_typ1[i1]] ~ dnorm(logit_pi1_mean, logit_pi1_prec)

  ## status sampled from probability - not sampled in this version
  Status[status_typ1[i1]] ~ dbern(ilogit(logit_pi[status_typ1[i1]]))

  ## test result
  p_test_pos_typ1[i1] <- Se[test_id_typ1[i1]] * Status[status_typ1[i1]] +
                          (1 - Sp[test_id_typ1[i1]]) * (1 - Status[status_typ1[i1]])

  test_res_typ1[i1] ~ dbern(p_test_pos_typ1[i1])

  }

  ## 2: first test on a month which is not first test in herd
  ## status type = 2
  for(i2 in 1:n_status_typ2){

  # probability of infection given previous status and dynamics
  pi[status_typ2[i2]] <- tau1 * (1 - Status[status_typ2[i2] - 1]) +
                         tau2 * Status[status_typ2[i2] - 1]

  ## status sampled from probability
  Status[status_typ2[i2]] ~ dbern(pi[status_typ2[i2]])

  ## test result
  p_test_pos_typ2[i2] <- Se[test_id_typ2[i2]] * Status[status_typ2[i2]] +
                          (1 - Sp[test_id_typ2[i2]]) * (1 - Status[status_typ2[i2]])

  test_res_typ2[i2] ~ dbern(p_test_pos_typ2[i2])

  }',

                 if(n_status_typ3 == 0){''} else {'
  ## 3: test > 1 on a month
  ## status type = 3
  for(i3 in 1:n_status_typ3){

  # probability of infection is given by previous test result on the same month
  pi[status_typ3[i3]] <- pi[status_typ3[i3] - 1]

  ## status sampled from probability
  Status[status_typ3[i3]] ~ dbern(pi[status_typ3[i3]])

  ## test result
  p_test_pos_typ3[i3] <- Se[test_id_typ3[i3]] * Status[status_typ3[i3]] +
                          (1 - Sp[test_id_typ3[i3]]) * (1 - Status[status_typ3[i3]])

  test_res_typ3[i3] ~ dbern(p_test_pos_typ3[i3])

  }
'},
                 '## no test
  ## no_test
  for(j in 1:n_no_test){

  # probability of infection given previous status and dynamics
  pi[status_no_test[j]] <- tau1 * (1 - Status[status_no_test[j] - 1]) +
                           tau2 * Status[status_no_test[j] - 1]

  # # status sampled from probability
  Status[status_no_test[j]] ~ dbern(pi[status_no_test[j]])

  }

##############################################################################
###  Prediction of probability of infection
##############################################################################

',

                 if(n_status_typ4 == 0){''} else {'
## 4: status to predict without test result
for(i4 in 1:n_status_typ4){

  # probability of infection given previous status and dynamics
  predicted_proba[herd_id_pr4[i4]] <- tau1 * (1 - Status[status_typ4[i4] - 1]) +
                 tau2 * Status[status_typ4[i4] - 1]

  # status sampled from probability
  predicted_status[herd_id_pr4[i4]] ~ dbern(predicted_proba[herd_id_pr4[i4]])

  }'},

                 if(n_status_typ5 == 0){''} else {'

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

  }'},

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
###  Loop for monthly prevalences
##############################################################################
for(i_month in 1:month_max){

 month_prev[i_month] <- sum(Status[month_mat[1:month_N[i_month], i_month]]) /
                            month_N[i_month]

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
  logit_tau1 ~ dnorm(logit_tau1_mean, logit_tau1_prec)
  logit_tau2 ~ dnorm(logit_tau2_mean, logit_tau2_prec)

  tau1 <- ilogit(logit_tau1)
  tau2 <- ilogit(logit_tau2)

  }
')

}




## JAGS model: herd level, with risk factors
write_JAGS_model.herd_rf <- function(data){

  test_data <- data$test_data
  ## several tests on the same month
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ4 <- nrow(test_data[test_data$status_type == 4,])
  n_status_typ5 <- nrow(test_data[test_data$status_type == 5,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])

model <- paste('model{
##############################################################################
###  Inference from historical data
##############################################################################

  ## First test in a herd - status type = 1
  for(i1 in 1:n_status_typ1){

  ## prior probability of infection on first test
  pi[status_typ1[i1]] ~ dbeta(pi1_beta_a, pi1_beta_b)

  ## status sampled from probability - not sampled in this version
  Status[status_typ1[i1]] ~ dbern(pi[status_typ1[i1]])

  ## test result
  p_test_pos_typ1[i1] <- Se[test_id_typ1[i1]] * Status[status_typ1[i1]] +
                          (1 - Sp[test_id_typ1[i1]]) * (1 - Status[status_typ1[i1]])

  test_res_typ1[i1] ~ dbern(p_test_pos_typ1[i1])

  }

  ## 2: first test on a month which is not first test in herd
  ## status type = 2
  for(i2 in 1:n_status_typ2){

  # probability of new infection
  logit(tau1[status_typ2[i2]]) <- inprod(risk_factors[status_typ2[i2],], theta)

  # probability of infection given previous status and dynamics
  pi[status_typ2[i2]] <- tau1[status_typ2[i2]] * (1 - Status[status_typ2[i2] - 1]) +
                         tau2 * Status[status_typ2[i2] - 1]

  ## status sampled from probability
  Status[status_typ2[i2]] ~ dbern(pi[status_typ2[i2]])

  ## test result
  p_test_pos_typ2[i2] <- Se[test_id_typ2[i2]] * Status[status_typ2[i2]] +
                          (1 - Sp[test_id_typ2[i2]]) * (1 - Status[status_typ2[i2]])

  test_res_typ2[i2] ~ dbern(p_test_pos_typ2[i2])

  }',

  if(n_status_typ3 == 0){''} else {'
  ## 3: test > 1 on a month
  ## status type = 3
  for(i3 in 1:n_status_typ3){

  # probability of infection is given by previous test result on the same month
  pi[status_typ3[i3]] <- pi[status_typ3[i3] - 1]

  ## status sampled from probability
  Status[status_typ3[i3]] ~ dbern(pi[status_typ3[i3]])

  ## test result
  p_test_pos_typ3[i3] <- Se[test_id_typ3[i3]] * Status[status_typ3[i3]] +
                          (1 - Sp[test_id_typ3[i3]]) * (1 - Status[status_typ3[i3]])

  test_res_typ3[i3] ~ dbern(p_test_pos_typ3[i3])

  }'},
  '## no test
  ## no_test
  for(j in 1:n_no_test){

  # probability of new infection
  logit(tau1[status_no_test[j]]) <- inprod(risk_factors[status_no_test[j],], theta)

  # probability of infection given previous status and dynamics
  pi[status_no_test[j]] <- tau1[status_no_test[j]] * (1 - Status[status_no_test[j] - 1]) +
                           tau2 * Status[status_no_test[j] - 1]

  ## status sampled from probability
  Status[status_no_test[j]] ~ dbern(pi[status_no_test[j]])

  }

##############################################################################
###  Prediction of probability of infection
##############################################################################
',
if(n_status_typ4 == 0){''} else {'
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

'},
if(n_status_typ5 == 0){''} else {'
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

  }'},
if(n_status_typ6 == 0){''} else {'
## 6: status to predict with several tests on this month
## same as above except that pi_init is probability of infection after previous test
## no dynamics
for(i6 in 1:n_status_typ6){

  # probability of infection given previous status and dynamics
  pi6_init[i6] <- pi[status_typ6[i6] - 1]

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
###  Loop for monthly prevalences
##############################################################################
for(i_month in 1:month_max){

 month_prev[i_month] <- sum(Status[month_mat[1:month_N[i_month], i_month]]) /
                            month_N[i_month]

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

}')

}


## JAGS model: herd level, with risk factors
## Dynamics parameters modelled on the logit scale
write_JAGS_model.herd_dynLogit_rf <- function(data){

  test_data <- data$test_data
  ## several tests on the same month
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ4 <- nrow(test_data[test_data$status_type == 4,])
  n_status_typ5 <- nrow(test_data[test_data$status_type == 5,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])

  model <- paste('model{
##############################################################################
###  Inference from historical data
##############################################################################

  ## First test in a herd - status type = 1
  for(i1 in 1:n_status_typ1){

  ## prior probability of infection on first test
  logit_pi[status_typ1[i1]] ~ dnorm(logit_pi1_mean, logit_pi1_prec)

  ## status sampled from probability - not sampled in this version
  Status[status_typ1[i1]] ~ dbern(ilogit(logit_pi[status_typ1[i1]]))

  ## test result
  p_test_pos_typ1[i1] <- Se[test_id_typ1[i1]] * Status[status_typ1[i1]] +
                          (1 - Sp[test_id_typ1[i1]]) * (1 - Status[status_typ1[i1]])

  test_res_typ1[i1] ~ dbern(p_test_pos_typ1[i1])

  }

  ## 2: first test on a month which is not first test in herd
  ## status type = 2
  for(i2 in 1:n_status_typ2){

  # probability of new infection
  logit(tau1[status_typ2[i2]]) <- inprod(risk_factors[status_typ2[i2],], theta)

  # probability of infection given previous status and dynamics
  pi[status_typ2[i2]] <- tau1[status_typ2[i2]] * (1 - Status[status_typ2[i2] - 1]) +
                         tau2 * Status[status_typ2[i2] - 1]

  ## status sampled from probability
  Status[status_typ2[i2]] ~ dbern(pi[status_typ2[i2]])

  ## test result
  p_test_pos_typ2[i2] <- Se[test_id_typ2[i2]] * Status[status_typ2[i2]] +
                          (1 - Sp[test_id_typ2[i2]]) * (1 - Status[status_typ2[i2]])

  test_res_typ2[i2] ~ dbern(p_test_pos_typ2[i2])

  }',

                 if(n_status_typ3 == 0){''} else {'
  ## 3: test > 1 on a month
  ## status type = 3
  for(i3 in 1:n_status_typ3){

  # probability of infection is given by previous test result on the same month
  pi[status_typ3[i3]] <- pi[status_typ3[i3] - 1]

  ## status sampled from probability
  Status[status_typ3[i3]] ~ dbern(pi[status_typ3[i3]])

  ## test result
  p_test_pos_typ3[i3] <- Se[test_id_typ3[i3]] * Status[status_typ3[i3]] +
                          (1 - Sp[test_id_typ3[i3]]) * (1 - Status[status_typ3[i3]])

  test_res_typ3[i3] ~ dbern(p_test_pos_typ3[i3])

  }'},
                 '## no test
  ## no_test
  for(j in 1:n_no_test){

  # probability of new infection
  logit(tau1[status_no_test[j]]) <- inprod(risk_factors[status_no_test[j],], theta)

  # probability of infection given previous status and dynamics
  pi[status_no_test[j]] <- tau1[status_no_test[j]] * (1 - Status[status_no_test[j] - 1]) +
                           tau2 * Status[status_no_test[j] - 1]

  ## status sampled from probability
  Status[status_no_test[j]] ~ dbern(pi[status_no_test[j]])

  }

##############################################################################
###  Prediction of probability of infection
##############################################################################
',
                 if(n_status_typ4 == 0){''} else {'
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

'},
                 if(n_status_typ5 == 0){''} else {'
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

  }'},
                 if(n_status_typ6 == 0){''} else {'
## 6: status to predict with several tests on this month
## same as above except that pi_init is probability of infection after previous test
## no dynamics
for(i6 in 1:n_status_typ6){

  # probability of infection given previous status and dynamics
  pi6_init[i6] <- pi[status_typ6[i6] - 1]

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
###  Loop for monthly prevalences
##############################################################################
for(i_month in 1:month_max){

 month_prev[i_month] <- sum(Status[month_mat[1:month_N[i_month], i_month]]) /
                            month_N[i_month]

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
  logit_tau2 ~ dnorm(logit_tau2_mean, logit_tau2_prec)
  tau2 <- ilogit(logit_tau2)

  ## Logistic regression coefficients
  for(i_rf in 1:n_risk_factors){

    theta[i_rf] ~ dnorm(theta_norm_mean[i_rf], theta_norm_prec[i_rf])

  }

}')

}


## JAGS model: animal level, no risk factor
write_JAGS_model.animal <- function(data){

  test_data <- data$test_data
  ## several tests on the same month
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ4 <- nrow(test_data[test_data$status_type == 4,])
  n_status_typ5 <- nrow(test_data[test_data$status_type == 5,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])

  # Formula for the prediction of posterior probability of infection:
  #   D+: herd is infected
  #   d+: animal is infected
  #   pi_h: herd level prevalence
  #   pi_w: animal level prevalence in infected herds
  #   T+: test positive
  #
  #   p(D+) = p(D+|d+)*[p(d+|T+)p(T+) + p(d+|T-)p(T-)] +
  #           p(D+|d-)*[p(d-|T+)p(T+) + p(d-|T-)p(T-)]
  #
  #   p(D+|d+) = 1 -> as soon as 1 animal infected, the herd is infected
  #   p(D+|d-) = (1 - pi_w)*pi_h / [(1 - pi_w)*pi_h + (1 - pi_h)

model <- paste('model{

##############################################################################
###  Inference from historical data
##############################################################################

  ## First test in a herd - status type = 1
  for(i1 in 1:n_status_typ1){

  ## prior probability of herd infection on first test
  pi[status_typ1[i1]] ~ dbeta(pi1_beta_a, pi1_beta_b)

  ## status sampled from probability - not sampled in this version
  Status[status_typ1[i1]] ~ dbern(pi[status_typ1[i1]])

  ## test result
  p_test_pos_typ1[i1] <- Se[test_id_typ1[i1]] * Status[status_typ1[i1]] * pi_within +
                          (1 - Sp[test_id_typ1[i1]]) * (1 - Status[status_typ1[i1]] * pi_within)

  n_pos_typ1[i1] ~ dbin(p_test_pos_typ1[i1], n_tested_typ1[i1])

  }


  ## 2: first test on a month which is not first test in herd
  ## status type = 2
  for(i2 in 1:n_status_typ2){

  # probability of infection given previous status and dynamics
  pi[status_typ2[i2]] <- tau1 * (1 - Status[status_typ2[i2] - 1]) +
                         tau2 * Status[status_typ2[i2] - 1]

  ## status sampled from probability - not sampled in this version
  Status[status_typ2[i2]] ~ dbern(pi[status_typ2[i2]])

  ## test result
  p_test_pos_typ2[i2] <- Se[test_id_typ2[i2]] * Status[status_typ2[i2]] * pi_within +
                          (1 - Sp[test_id_typ2[i2]]) * (1 - Status[status_typ2[i2]] * pi_within)

  n_pos_typ2[i2] ~ dbin(p_test_pos_typ2[i2], n_tested_typ2[i2])

  }
',

  if(n_status_typ3 == 0){''} else {'
  ## 3: test > 1 on a month
  ## status type = 3
  for(i3 in 1:n_status_typ3){

  # probability of infection is given by previous test result on the same month
  pi[status_typ3[i3]] <- pi[status_typ3[i3] - 1]

  ## status sampled from probability - not sampled in this version
  Status[status_typ3[i3]] ~ dbern(pi[status_typ3[i3]])

  ## test result
  p_test_pos_typ3[i3] <- Se[test_id_typ3[i3]] * pi[status_typ3[i3]] * pi_within +
                          (1 - Sp[test_id_typ3[i3]]) * (1 - pi[status_typ2[i3]] * pi_within)

  n_pos_typ3[i3] ~ dbin(p_test_pos_typ3[i3], n_tested_typ3[i3])

  }'}, '

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
',

  if(n_status_typ4 == 0){''} else {'
## 4: status to predict without test result
for(i4 in 1:n_status_typ4){

  # probability of infection given previous status and dynamics
  predicted_proba[herd_id_pr4[i4]] <- tau1 * (1 - Status[status_typ4[i4] - 1]) +
                 tau2 * Status[status_typ4[i4] - 1]

  # status sampled from probability
  predicted_status[herd_id_pr4[i4]] ~ dbern(predicted_proba[herd_id_pr4[i4]])

  }'},

if(n_status_typ5 == 0){''} else {'

## 5: status to predict with a single test performed
for(i5 in 1:n_status_typ5){

  # probability of infection given previous status and dynamics
  pi5_init[i5] <- tau1 * (1 - Status[status_typ5[i5] - 1]) +
                  tau2 * Status[status_typ5[i5] - 1]

  ## estimation of proportion of positives from Binomial
  p_T1_i5[i5] ~ dbeta(1, 1)
  n_pos_typ5[i5] ~ dbin(p_T1_i5[i5], n_tested_typ5[i5])

  ## posterior probability of herd infection
  ### Overall prior probability of infection
  prev_i5[i5] <- pi5_init[i5] * pi_within
  ### p(D+|d-) probability that herd is infected even if no animal detected
  D1_d0_i5[i5] <- (1 - pi_within) * pi5_init[i5] /
                   (1 - prev_i5[i5])
  ### p(d+|T+)
  d1_T1_i5[i5] <- Se[test_id_typ5[i5]] * prev_i5[i5] /
                  (Se[test_id_typ5[i5]] * prev_i5[i5] +
                  (1 - Sp[test_id_typ5[i5]]) * (1 - prev_i5[i5]))
  ### p(d+|T-)
  d1_T0_i5[i5] <- (1 - Se[test_id_typ5[i5]]) * prev_i5[i5] /
                  ((1 - Se[test_id_typ5[i5]]) * prev_i5[i5] +
                  Sp[test_id_typ5[i5]] * (1- prev_i5[i5]))
  ### p(d-|T+)
  d0_T1_i5[i5] <- (1 - Sp[test_id_typ5[i5]]) * (1 - prev_i5[i5]) /
                  ((1 - Sp[test_id_typ5[i5]]) * (1 - prev_i5[i5]) +
                  Se[test_id_typ5[i5]] * prev_i5[i5])
  ### p(d-|T-)
  d0_T0_i5[i5] <- Sp[test_id_typ5[i5]] * (1 - prev_i5[i5]) /
                  (Sp[test_id_typ5[i5]] * (1 - prev_i5[i5]) +
                   (1 - Se[test_id_typ5[i5]]) * prev_i5[i5])

  ### herd level posterior probability of infection
  predicted_proba[herd_id_pr5[i5]] <- d1_T1_i5[i5] * p_T1_i5[i5] +
                         d1_T0_i5[i5] * (1 - p_T1_i5[i5]) +
                         D1_d0_i5[i5] * (
                         d0_T1_i5[i5] * p_T1_i5[i5] +
                         d0_T0_i5[i5] * (1 - p_T1_i5[i5])  )

  # status sampled from probability
  predicted_status[herd_id_pr5[i5]] ~ dbern(predicted_proba[herd_id_pr5[i5]])

  }'},

if(n_status_typ6 == 0){''} else {'

## 6: status to predict with several tests on this month
## same as above except that pi_init is probability of infection after previous test
## no dynamics
for(i6 in 1:n_status_typ6){

  # probability of infection given previous status and dynamics
  pi6_init[i6] <- pi[status_typ6[i6] - 1]

  ## estimation of proportion of positives from Binomial
  p_T1_i6[i6] ~ dbeta(1, 1)
  n_pos_typ6[i6] ~ dbin(p_T1_i6[i6], n_tested_typ6[i6])

  ## posterior probability of herd infection
  ### Overall prior probability of infection
  prev_i6[i6] <- pi6_init[i6] * pi_within
  ### p(D+|d-) probability that herd is infected even if no animal detected
  D1_d0_i6[i6] <- (1 - pi_within) * pi6_init[i6] /
                   (1 - prev_i6[i6])
  ### p(d+|T+)
  d1_T1_i6[i6] <- Se[test_id_typ6[i6]] * prev_i6[i6] /
                  (Se[test_id_typ6[i6]] * prev_i6[i6] +
                  (1 - Sp[test_id_typ6[i6]]) * (1 - prev_i6[i6]))
  ### p(d+|T-)
  d1_T0_i6[i6] <- (1 - Se[test_id_typ6[i6]]) * prev_i6[i6] /
                  ((1 - Se[test_id_typ6[i6]]) * prev_i6[i6] +
                  Sp[test_id_typ6[i6]] * (1- prev_i6[i6]))
  ### p(d-|T+)
  d0_T1_i6[i6] <- (1 - Sp[test_id_typ6[i6]]) * (1 - prev_i6[i6]) /
                  ((1 - Sp[test_id_typ6[i6]]) * (1 - prev_i6[i6]) +
                  Se[test_id_typ6[i6]] * prev_i6[i6])
  ### p(d-|T-)
  d0_T0_i6[i6] <- Sp[test_id_typ6[i6]] * (1 - prev_i6[i6]) /
                  (Sp[test_id_typ6[i6]] * (1 - prev_i6[i6]) +
                   (1 - Se[test_id_typ6[i6]]) * prev_i6[i6])

  ### herd level posterior probability of infection
  predicted_proba[herd_id_pr6[i6]] <- d1_T1_i6[i6] * p_T1_i6[i6] +
                         d1_T0_i6[i6] * (1 - p_T1_i6[i6]) +
                         D1_d0_i6[i6] * (
                         d0_T1_i6[i6] * p_T1_i6[i6] +
                         d0_T0_i6[i6] * (1 - p_T1_i6[i6])  )

  # status sampled from probability
  predicted_status[herd_id_pr6[i6]] ~ dbern(predicted_proba[herd_id_pr6[i6]])

  }'},
'
##############################################################################
###  Loop for monthly prevalences
##############################################################################
for(i_month in 1:month_max){

 month_prev[i_month] <- sum(Status[month_mat[1:month_N[i_month], i_month]]) /
                            month_N[i_month]

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
  tau1 ~ dbeta(tau1_beta_a, tau1_beta_b)
  tau2 ~ dbeta(tau2_beta_a, tau2_beta_b)


}')

}



## JAGS model: animal level, no risk factor
write_JAGS_model.animal_dynLogit <- function(data){

  test_data <- data$test_data
  ## several tests on the same month
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ4 <- nrow(test_data[test_data$status_type == 4,])
  n_status_typ5 <- nrow(test_data[test_data$status_type == 5,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])

  # Formula for the prediction of posterior probability of infection:
  #   D+: herd is infected
  #   d+: animal is infected
  #   pi_h: herd level prevalence
  #   pi_w: animal level prevalence in infected herds
  #   T+: test positive
  #
  #   p(D+) = p(D+|d+)*[p(d+|T+)p(T+) + p(d+|T-)p(T-)] +
  #           p(D+|d-)*[p(d-|T+)p(T+) + p(d-|T-)p(T-)]
  #
  #   p(D+|d+) = 1 -> as soon as 1 animal infected, the herd is infected
  #   p(D+|d-) = (1 - pi_w)*pi_h / [(1 - pi_w)*pi_h + (1 - pi_h)

  model <- paste('model{

##############################################################################
###  Inference from historical data
##############################################################################

  ## First test in a herd - status type = 1
  for(i1 in 1:n_status_typ1){

  ## prior probability of infection on first test
  logit_pi[status_typ1[i1]] ~ dnorm(logit_pi1_mean, logit_pi1_prec)

  ## status sampled from probability - not sampled in this version
  Status[status_typ1[i1]] ~ dbern(ilogit(logit_pi[status_typ1[i1]]))

  ## test result
  p_test_pos_typ1[i1] <- Se[test_id_typ1[i1]] * Status[status_typ1[i1]] * pi_within +
                          (1 - Sp[test_id_typ1[i1]]) * (1 - Status[status_typ1[i1]] * pi_within)

  n_pos_typ1[i1] ~ dbin(p_test_pos_typ1[i1], n_tested_typ1[i1])

  }


  ## 2: first test on a month which is not first test in herd
  ## status type = 2
  for(i2 in 1:n_status_typ2){

  # probability of infection given previous status and dynamics
  pi[status_typ2[i2]] <- tau1 * (1 - Status[status_typ2[i2] - 1]) +
                         tau2 * Status[status_typ2[i2] - 1]

  ## status sampled from probability - not sampled in this version
  Status[status_typ2[i2]] ~ dbern(pi[status_typ2[i2]])

  ## test result
  p_test_pos_typ2[i2] <- Se[test_id_typ2[i2]] * Status[status_typ2[i2]] * pi_within +
                          (1 - Sp[test_id_typ2[i2]]) * (1 - Status[status_typ2[i2]] * pi_within)

  n_pos_typ2[i2] ~ dbin(p_test_pos_typ2[i2], n_tested_typ2[i2])

  }
',

                 if(n_status_typ3 == 0){''} else {'
  ## 3: test > 1 on a month
  ## status type = 3
  for(i3 in 1:n_status_typ3){

  # probability of infection is given by previous test result on the same month
  pi[status_typ3[i3]] <- pi[status_typ3[i3] - 1]

  ## status sampled from probability - not sampled in this version
  Status[status_typ3[i3]] ~ dbern(pi[status_typ3[i3]])

  ## test result
  p_test_pos_typ3[i3] <- Se[test_id_typ3[i3]] * pi[status_typ3[i3]] * pi_within +
                          (1 - Sp[test_id_typ3[i3]]) * (1 - pi[status_typ2[i3]] * pi_within)

  n_pos_typ3[i3] ~ dbin(p_test_pos_typ3[i3], n_tested_typ3[i3])

  }'}, '

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
',

                 if(n_status_typ4 == 0){''} else {'
## 4: status to predict without test result
for(i4 in 1:n_status_typ4){

  # probability of infection given previous status and dynamics
  predicted_proba[herd_id_pr4[i4]] <- tau1 * (1 - Status[status_typ4[i4] - 1]) +
                 tau2 * Status[status_typ4[i4] - 1]

  # status sampled from probability
  predicted_status[herd_id_pr4[i4]] ~ dbern(predicted_proba[herd_id_pr4[i4]])

  }'},

                 if(n_status_typ5 == 0){''} else {'

## 5: status to predict with a single test performed
for(i5 in 1:n_status_typ5){

  # probability of infection given previous status and dynamics
  pi5_init[i5] <- tau1 * (1 - Status[status_typ5[i5] - 1]) +
                  tau2 * Status[status_typ5[i5] - 1]

  ## estimation of proportion of positives from Binomial
  p_T1_i5[i5] ~ dbeta(1, 1)
  n_pos_typ5[i5] ~ dbin(p_T1_i5[i5], n_tested_typ5[i5])

  ## posterior probability of herd infection
  ### Overall prior probability of infection
  prev_i5[i5] <- pi5_init[i5] * pi_within
  ### p(D+|d-) probability that herd is infected even if no animal detected
  D1_d0_i5[i5] <- (1 - pi_within) * pi5_init[i5] /
                   (1 - prev_i5[i5])
  ### p(d+|T+)
  d1_T1_i5[i5] <- Se[test_id_typ5[i5]] * prev_i5[i5] /
                  (Se[test_id_typ5[i5]] * prev_i5[i5] +
                  (1 - Sp[test_id_typ5[i5]]) * (1 - prev_i5[i5]))
  ### p(d+|T-)
  d1_T0_i5[i5] <- (1 - Se[test_id_typ5[i5]]) * prev_i5[i5] /
                  ((1 - Se[test_id_typ5[i5]]) * prev_i5[i5] +
                  Sp[test_id_typ5[i5]] * (1- prev_i5[i5]))
  ### p(d-|T+)
  d0_T1_i5[i5] <- (1 - Sp[test_id_typ5[i5]]) * (1 - prev_i5[i5]) /
                  ((1 - Sp[test_id_typ5[i5]]) * (1 - prev_i5[i5]) +
                  Se[test_id_typ5[i5]] * prev_i5[i5])
  ### p(d-|T-)
  d0_T0_i5[i5] <- Sp[test_id_typ5[i5]] * (1 - prev_i5[i5]) /
                  (Sp[test_id_typ5[i5]] * (1 - prev_i5[i5]) +
                   (1 - Se[test_id_typ5[i5]]) * prev_i5[i5])

  ### herd level posterior probability of infection
  predicted_proba[herd_id_pr5[i5]] <- d1_T1_i5[i5] * p_T1_i5[i5] +
                         d1_T0_i5[i5] * (1 - p_T1_i5[i5]) +
                         D1_d0_i5[i5] * (
                         d0_T1_i5[i5] * p_T1_i5[i5] +
                         d0_T0_i5[i5] * (1 - p_T1_i5[i5])  )

  # status sampled from probability
  predicted_status[herd_id_pr5[i5]] ~ dbern(predicted_proba[herd_id_pr5[i5]])

  }'},

                 if(n_status_typ6 == 0){''} else {'

## 6: status to predict with several tests on this month
## same as above except that pi_init is probability of infection after previous test
## no dynamics
for(i6 in 1:n_status_typ6){

  # probability of infection given previous status and dynamics
  pi6_init[i6] <- pi[status_typ6[i6] - 1]

  ## estimation of proportion of positives from Binomial
  p_T1_i6[i6] ~ dbeta(1, 1)
  n_pos_typ6[i6] ~ dbin(p_T1_i6[i6], n_tested_typ6[i6])

  ## posterior probability of herd infection
  ### Overall prior probability of infection
  prev_i6[i6] <- pi6_init[i6] * pi_within
  ### p(D+|d-) probability that herd is infected even if no animal detected
  D1_d0_i6[i6] <- (1 - pi_within) * pi6_init[i6] /
                   (1 - prev_i6[i6])
  ### p(d+|T+)
  d1_T1_i6[i6] <- Se[test_id_typ6[i6]] * prev_i6[i6] /
                  (Se[test_id_typ6[i6]] * prev_i6[i6] +
                  (1 - Sp[test_id_typ6[i6]]) * (1 - prev_i6[i6]))
  ### p(d+|T-)
  d1_T0_i6[i6] <- (1 - Se[test_id_typ6[i6]]) * prev_i6[i6] /
                  ((1 - Se[test_id_typ6[i6]]) * prev_i6[i6] +
                  Sp[test_id_typ6[i6]] * (1- prev_i6[i6]))
  ### p(d-|T+)
  d0_T1_i6[i6] <- (1 - Sp[test_id_typ6[i6]]) * (1 - prev_i6[i6]) /
                  ((1 - Sp[test_id_typ6[i6]]) * (1 - prev_i6[i6]) +
                  Se[test_id_typ6[i6]] * prev_i6[i6])
  ### p(d-|T-)
  d0_T0_i6[i6] <- Sp[test_id_typ6[i6]] * (1 - prev_i6[i6]) /
                  (Sp[test_id_typ6[i6]] * (1 - prev_i6[i6]) +
                   (1 - Se[test_id_typ6[i6]]) * prev_i6[i6])

  ### herd level posterior probability of infection
  predicted_proba[herd_id_pr6[i6]] <- d1_T1_i6[i6] * p_T1_i6[i6] +
                         d1_T0_i6[i6] * (1 - p_T1_i6[i6]) +
                         D1_d0_i6[i6] * (
                         d0_T1_i6[i6] * p_T1_i6[i6] +
                         d0_T0_i6[i6] * (1 - p_T1_i6[i6])  )

  # status sampled from probability
  predicted_status[herd_id_pr6[i6]] ~ dbern(predicted_proba[herd_id_pr6[i6]])

  }'},
                 '
##############################################################################
###  Loop for monthly prevalences
##############################################################################
for(i_month in 1:month_max){

 month_prev[i_month] <- sum(Status[month_mat[1:month_N[i_month], i_month]]) /
                            month_N[i_month]

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
  logit_pi_within ~ dnorm(logit_pi_within_mean, logit_pi_within_prec)
  pi_within <- ilogit(logit_pi_within)

  ## Probability of not eliminating the infection
  logit_tau1 ~ dnorm(logit_tau1_mean, logit_tau1_prec)
  logit_tau2 ~ dnorm(logit_tau2_mean, logit_tau2_prec)

  tau1 <- ilogit(logit_tau1)
  tau2 <- ilogit(logit_tau2)


}')

}


## JAGS model: animal level, with risk factors
write_JAGS_model.animal_rf <- function(data){

  test_data <- data$test_data
  ## several tests on the same month
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ4 <- nrow(test_data[test_data$status_type == 4,])
  n_status_typ5 <- nrow(test_data[test_data$status_type == 5,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])

model <- paste('model{
  # Formula for the prediction of posterior probability of infection:
  #   D+: herd is infected
  #   d+: animal is infected
  #   pi_h: herd level prevalence
  #   pi_w: animal level prevalence in infected herds
  #   T+: test positive
  #
  #   p(D+) = p(D+|d+)*[p(d+|T+)p(T+) + p(d+|T-)p(T-)] +
  #           p(D+|d-)*[p(d-|T+)p(T+) + p(d-|T-)p(T-)]
  #
  #   p(D+|d+) = 1 -> as soon as 1 animal infected, the herd is infected
  #   p(D+|d-) = (1 - pi_w)*pi_h / [(1 - pi_w)*pi_h + (1 - pi_h)

  ##############################################################################
  ###  Inference from historical data
  ##############################################################################

  ## First test in a herd - status type = 1
  for(i1 in 1:n_status_typ1){

  ## prior probability of herd infection on first test
  pi[status_typ1[i1]] ~ dbeta(pi1_beta_a, pi1_beta_b)

  ## status sampled from probability - not sampled in this version
  Status[status_typ1[i1]] ~ dbern(pi[status_typ1[i1]])

  ## test result
  p_test_pos_typ1[i1] <- Se[test_id_typ1[i1]] * Status[status_typ1[i1]] * pi_within +
                          (1 - Sp[test_id_typ1[i1]]) * (1 - Status[status_typ1[i1]] * pi_within)

  n_pos_typ1[i1] ~ dbin(p_test_pos_typ1[i1], n_tested_typ1[i1])
 }


  ## 2: first test on a month which is not first test in herd
  ## status type = 2
  for(i2 in 1:n_status_typ2){

  # probability of new infection
  logit(tau1[status_typ2[i2]]) <- inprod(risk_factors[status_typ2[i2],], theta)

  # probability of infection given previous status and dynamics
  pi[status_typ2[i2]] <- tau1[status_typ2[i2]] * (1 - Status[status_typ2[i2] - 1]) +
                         tau2 * Status[status_typ2[i2] - 1]

    ## status sampled from probability - not sampled in this version
    Status[status_typ2[i2]] ~ dbern(pi[status_typ2[i2]])

    ## test result
    p_test_pos_typ2[i2] <- Se[test_id_typ2[i2]] * Status[status_typ2[i2]] * pi_within +
      (1 - Sp[test_id_typ2[i2]]) * (1 - Status[status_typ2[i2]] * pi_within)

    n_pos_typ2[i2] ~ dbin(p_test_pos_typ2[i2], n_tested_typ2[i2])

  }

',

if(n_status_typ3 == 0){''} else {'
  ## 3: test > 1 on a month
  ## status type = 3
  for(i3 in 1:n_status_typ3){

    # probability of infection is given by previous test result on the same month
    pi[status_typ3[i3]] <- pi[status_typ3[i3] - 1]

    ## status sampled from probability - not sampled in this version
    Status[status_typ3[i3]] ~ dbern(pi[status_typ3[i3]])

    ## test result
    p_test_pos_typ3[i3] <- Se[test_id_typ3[i3]] * pi[status_typ3[i3]] * pi_within +
      (1 - Sp[test_id_typ3[i3]]) * (1 - pi[status_typ2[i3]] * pi_within)

    n_pos_typ3[i3] ~ dbin(p_test_pos_typ3[i3], n_tested_typ3[i3])

  }'},
'
  ## no test
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
  ',

if(n_status_typ4 == 0){''} else {'
## 4: status to predict without test result
  for(i4 in 1:n_status_typ4){

  # probability of new infection
  logit(tau1[status_no_test[j]]) <- inprod(risk_factors[status_no_test[j],], theta)

  # probability of infection given previous status and dynamics
  pi[status_no_test[j]] <- tau1[status_no_test[j]] * (1 - Status[status_no_test[j] - 1]) +
    tau2 * Status[status_no_test[j] - 1]

    # status sampled from probability
    predicted_status[herd_id_pr4[i4]] ~ dbern(predicted_proba[herd_id_pr4[i4]])

  }'},

  if(n_status_typ5 == 0){''} else {'

  ## 5: status to predict with a single test performed
  for(i5 in 1:n_status_typ5){

  # probability of new infection
  logit(tau1[status_typ5[i5]]) <- inprod(risk_factors[status_typ5[i5],], theta)

  # probability of infection given previous status and dynamics
  pi5_init[i5] <- tau1[status_typ5[i5]] * (1 - Status[status_typ5[i5] - 1]) + tau2 * Status[status_typ5[i5] - 1]

    ## estimation of proportion of positives from Binomial
    p_T1_i5[i5] ~ dbeta(1, 1)
    n_pos_typ5[i5] ~ dbin(p_T1_i5[i5], n_tested_typ5[i5])

    ## posterior probability of herd infection
    ### Overall prior probability of infection
    prev_i5[i5] <- pi5_init[i5] * pi_within
    ### p(D+|d-) probability that herd is infected even if no animal detected
    D1_d0_i5[i5] <- (1 - pi_within) * pi5_init[i5] /
      (1 - prev_i5[i5])
    ### p(d+|T+)
    d1_T1_i5[i5] <- Se[test_id_typ5[i5]] * prev_i5[i5] /
      (Se[test_id_typ5[i5]] * prev_i5[i5] +
         (1 - Sp[test_id_typ5[i5]]) * (1 - prev_i5[i5]))
    ### p(d+|T-)
    d1_T0_i5[i5] <- (1 - Se[test_id_typ5[i5]]) * prev_i5[i5] /
      ((1 - Se[test_id_typ5[i5]]) * prev_i5[i5] +
         Sp[test_id_typ5[i5]] * (1- prev_i5[i5]))
    ### p(d-|T+)
    d0_T1_i5[i5] <- (1 - Sp[test_id_typ5[i5]]) * (1 - prev_i5[i5]) /
      ((1 - Sp[test_id_typ5[i5]]) * (1 - prev_i5[i5]) +
         Se[test_id_typ5[i5]] * prev_i5[i5])
    ### p(d-|T-)
    d0_T0_i5[i5] <- Sp[test_id_typ5[i5]] * (1 - prev_i5[i5]) /
      (Sp[test_id_typ5[i5]] * (1 - prev_i5[i5]) +
         (1 - Se[test_id_typ5[i5]]) * prev_i5[i5])

  ### herd level posterior probability of infection
  predicted_proba[herd_id_pr5[i5]] <- d1_T1_i5[i5] * p_T1_i5[i5] +
                         d1_T0_i5[i5] * (1 - p_T1_i5[i5]) +
                         D1_d0_i5[i5] * (
                         d0_T1_i5[i5] * p_T1_i5[i5] +
                         d0_T0_i5[i5] * (1 - p_T1_i5[i5])  )

  # status sampled from probability
  predicted_status[herd_id_pr5[i5]] ~ dbern(predicted_proba[herd_id_pr5[i5]])

  }'},

 if(n_status_typ6 == 0){''} else {'

  ## 6: status to predict with several tests on this month
  ## same as above except that pi_init is probability of infection after previous test
  ## no dynamics
  for(i6 in 1:n_status_typ6){

    # probability of infection given previous status and dynamics
    pi6_init[i6] <- pi[status_typ6[i6] - 1]

    ## estimation of proportion of positives from Binomial
    p_T1_i6[i6] ~ dbeta(1, 1)
    n_pos_typ6[i6] ~ dbin(p_T1_i6[i6], n_tested_typ6[i6])

    ## posterior probability of herd infection
    ### Overall prior probability of infection
    prev_i6[i6] <- pi6_init[i6] * pi_within
    ### p(D+|d-) probability that herd is infected even if no animal detected
    D1_d0_i6[i6] <- (1 - pi_within) * pi6_init[i6] /
      (1 - prev_i6[i6])
    ### p(d+|T+)
    d1_T1_i6[i6] <- Se[test_id_typ6[i6]] * prev_i6[i6] /
      (Se[test_id_typ6[i6]] * prev_i6[i6] +
         (1 - Sp[test_id_typ6[i6]]) * (1 - prev_i6[i6]))
    ### p(d+|T-)
    d1_T0_i6[i6] <- (1 - Se[test_id_typ6[i6]]) * prev_i6[i6] /
      ((1 - Se[test_id_typ6[i6]]) * prev_i6[i6] +
         Sp[test_id_typ6[i6]] * (1- prev_i6[i6]))
    ### p(d-|T+)
    d0_T1_i6[i6] <- (1 - Sp[test_id_typ6[i6]]) * (1 - prev_i6[i6]) /
      ((1 - Sp[test_id_typ6[i6]]) * (1 - prev_i6[i6]) +
         Se[test_id_typ6[i6]] * prev_i6[i6])
    ### p(d-|T-)
    d0_T0_i6[i6] <- Sp[test_id_typ6[i6]] * (1 - prev_i6[i6]) /
      (Sp[test_id_typ6[i6]] * (1 - prev_i6[i6]) +
         (1 - Se[test_id_typ6[i6]]) * prev_i6[i6])

  ### herd level posterior probability of infection
  predicted_proba[herd_id_pr6[i6]] <- d1_T1_i6[i6] * p_T1_i6[i6] +
                         d1_T0_i6[i6] * (1 - p_T1_i6[i6]) +
                         D1_d0_i6[i6] * (
                         d0_T1_i6[i6] * p_T1_i6[i6] +
                         d0_T0_i6[i6] * (1 - p_T1_i6[i6])  )

  # status sampled from probability
  predicted_status[herd_id_pr6[i6]] ~ dbern(predicted_proba[herd_id_pr6[i6]])

  }'},
'
##############################################################################
###  Loop for monthly prevalences
##############################################################################
for(i_month in 1:month_max){

 month_prev[i_month] <- sum(Status[month_mat[1:month_N[i_month], i_month]]) /
                            month_N[i_month]

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

}')

}




## JAGS model: animal level, with risk factors
write_JAGS_model.animal_dynLogit_rf <- function(data){

  test_data <- data$test_data
  ## several tests on the same month
  n_status_typ3 <- nrow(test_data[test_data$status_type == 3,])
  n_status_typ4 <- nrow(test_data[test_data$status_type == 4,])
  n_status_typ5 <- nrow(test_data[test_data$status_type == 5,])
  n_status_typ6 <- nrow(test_data[test_data$status_type == 6,])

  model <- paste('model{
  # Formula for the prediction of posterior probability of infection:
  #   D+: herd is infected
  #   d+: animal is infected
  #   pi_h: herd level prevalence
  #   pi_w: animal level prevalence in infected herds
  #   T+: test positive
  #
  #   p(D+) = p(D+|d+)*[p(d+|T+)p(T+) + p(d+|T-)p(T-)] +
  #           p(D+|d-)*[p(d-|T+)p(T+) + p(d-|T-)p(T-)]
  #
  #   p(D+|d+) = 1 -> as soon as 1 animal infected, the herd is infected
  #   p(D+|d-) = (1 - pi_w)*pi_h / [(1 - pi_w)*pi_h + (1 - pi_h)

  ##############################################################################
  ###  Inference from historical data
  ##############################################################################

  ## First test in a herd - status type = 1
  for(i1 in 1:n_status_typ1){

  ## prior probability of infection on first test
  logit_pi[status_typ1[i1]] ~ dnorm(logit_pi1_mean, logit_pi1_prec)

  ## status sampled from probability - not sampled in this version
  Status[status_typ1[i1]] ~ dbern(ilogit(logit_pi[status_typ1[i1]]))

  ## test result
  p_test_pos_typ1[i1] <- Se[test_id_typ1[i1]] * Status[status_typ1[i1]] * pi_within +
                          (1 - Sp[test_id_typ1[i1]]) * (1 - Status[status_typ1[i1]] * pi_within)

  n_pos_typ1[i1] ~ dbin(p_test_pos_typ1[i1], n_tested_typ1[i1])
 }


  ## 2: first test on a month which is not first test in herd
  ## status type = 2
  for(i2 in 1:n_status_typ2){

  # probability of new infection
  logit(tau1[status_typ2[i2]]) <- inprod(risk_factors[status_typ2[i2],], theta)

  # probability of infection given previous status and dynamics
  pi[status_typ2[i2]] <- tau1[status_typ2[i2]] * (1 - Status[status_typ2[i2] - 1]) +
                         tau2 * Status[status_typ2[i2] - 1]

    ## status sampled from probability - not sampled in this version
    Status[status_typ2[i2]] ~ dbern(pi[status_typ2[i2]])

    ## test result
    p_test_pos_typ2[i2] <- Se[test_id_typ2[i2]] * Status[status_typ2[i2]] * pi_within +
      (1 - Sp[test_id_typ2[i2]]) * (1 - Status[status_typ2[i2]] * pi_within)

    n_pos_typ2[i2] ~ dbin(p_test_pos_typ2[i2], n_tested_typ2[i2])

  }

',

                 if(n_status_typ3 == 0){''} else {'
  ## 3: test > 1 on a month
  ## status type = 3
  for(i3 in 1:n_status_typ3){

    # probability of infection is given by previous test result on the same month
    pi[status_typ3[i3]] <- pi[status_typ3[i3] - 1]

    ## status sampled from probability - not sampled in this version
    Status[status_typ3[i3]] ~ dbern(pi[status_typ3[i3]])

    ## test result
    p_test_pos_typ3[i3] <- Se[test_id_typ3[i3]] * pi[status_typ3[i3]] * pi_within +
      (1 - Sp[test_id_typ3[i3]]) * (1 - pi[status_typ2[i3]] * pi_within)

    n_pos_typ3[i3] ~ dbin(p_test_pos_typ3[i3], n_tested_typ3[i3])

  }'},
                 '
  ## no test
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
  ',

                 if(n_status_typ4 == 0){''} else {'
## 4: status to predict without test result
  for(i4 in 1:n_status_typ4){

  # probability of new infection
  logit(tau1[status_no_test[j]]) <- inprod(risk_factors[status_no_test[j],], theta)

  # probability of infection given previous status and dynamics
  pi[status_no_test[j]] <- tau1[status_no_test[j]] * (1 - Status[status_no_test[j] - 1]) +
    tau2 * Status[status_no_test[j] - 1]

    # status sampled from probability
    predicted_status[herd_id_pr4[i4]] ~ dbern(predicted_proba[herd_id_pr4[i4]])

  }'},

                 if(n_status_typ5 == 0){''} else {'

  ## 5: status to predict with a single test performed
  for(i5 in 1:n_status_typ5){

  # probability of new infection
  logit(tau1[status_typ5[i5]]) <- inprod(risk_factors[status_typ5[i5],], theta)

  # probability of infection given previous status and dynamics
  pi5_init[i5] <- tau1[status_typ5[i5]] * (1 - Status[status_typ5[i5] - 1]) + tau2 * Status[status_typ5[i5] - 1]

    ## estimation of proportion of positives from Binomial
    p_T1_i5[i5] ~ dbeta(1, 1)
    n_pos_typ5[i5] ~ dbin(p_T1_i5[i5], n_tested_typ5[i5])

    ## posterior probability of herd infection
    ### Overall prior probability of infection
    prev_i5[i5] <- pi5_init[i5] * pi_within
    ### p(D+|d-) probability that herd is infected even if no animal detected
    D1_d0_i5[i5] <- (1 - pi_within) * pi5_init[i5] /
      (1 - prev_i5[i5])
    ### p(d+|T+)
    d1_T1_i5[i5] <- Se[test_id_typ5[i5]] * prev_i5[i5] /
      (Se[test_id_typ5[i5]] * prev_i5[i5] +
         (1 - Sp[test_id_typ5[i5]]) * (1 - prev_i5[i5]))
    ### p(d+|T-)
    d1_T0_i5[i5] <- (1 - Se[test_id_typ5[i5]]) * prev_i5[i5] /
      ((1 - Se[test_id_typ5[i5]]) * prev_i5[i5] +
         Sp[test_id_typ5[i5]] * (1- prev_i5[i5]))
    ### p(d-|T+)
    d0_T1_i5[i5] <- (1 - Sp[test_id_typ5[i5]]) * (1 - prev_i5[i5]) /
      ((1 - Sp[test_id_typ5[i5]]) * (1 - prev_i5[i5]) +
         Se[test_id_typ5[i5]] * prev_i5[i5])
    ### p(d-|T-)
    d0_T0_i5[i5] <- Sp[test_id_typ5[i5]] * (1 - prev_i5[i5]) /
      (Sp[test_id_typ5[i5]] * (1 - prev_i5[i5]) +
         (1 - Se[test_id_typ5[i5]]) * prev_i5[i5])

  ### herd level posterior probability of infection
  predicted_proba[herd_id_pr5[i5]] <- d1_T1_i5[i5] * p_T1_i5[i5] +
                         d1_T0_i5[i5] * (1 - p_T1_i5[i5]) +
                         D1_d0_i5[i5] * (
                         d0_T1_i5[i5] * p_T1_i5[i5] +
                         d0_T0_i5[i5] * (1 - p_T1_i5[i5])  )

  # status sampled from probability
  predicted_status[herd_id_pr5[i5]] ~ dbern(predicted_proba[herd_id_pr5[i5]])

  }'},

                 if(n_status_typ6 == 0){''} else {'

  ## 6: status to predict with several tests on this month
  ## same as above except that pi_init is probability of infection after previous test
  ## no dynamics
  for(i6 in 1:n_status_typ6){

    # probability of infection given previous status and dynamics
    pi6_init[i6] <- pi[status_typ6[i6] - 1]

    ## estimation of proportion of positives from Binomial
    p_T1_i6[i6] ~ dbeta(1, 1)
    n_pos_typ6[i6] ~ dbin(p_T1_i6[i6], n_tested_typ6[i6])

    ## posterior probability of herd infection
    ### Overall prior probability of infection
    prev_i6[i6] <- pi6_init[i6] * pi_within
    ### p(D+|d-) probability that herd is infected even if no animal detected
    D1_d0_i6[i6] <- (1 - pi_within) * pi6_init[i6] /
      (1 - prev_i6[i6])
    ### p(d+|T+)
    d1_T1_i6[i6] <- Se[test_id_typ6[i6]] * prev_i6[i6] /
      (Se[test_id_typ6[i6]] * prev_i6[i6] +
         (1 - Sp[test_id_typ6[i6]]) * (1 - prev_i6[i6]))
    ### p(d+|T-)
    d1_T0_i6[i6] <- (1 - Se[test_id_typ6[i6]]) * prev_i6[i6] /
      ((1 - Se[test_id_typ6[i6]]) * prev_i6[i6] +
         Sp[test_id_typ6[i6]] * (1- prev_i6[i6]))
    ### p(d-|T+)
    d0_T1_i6[i6] <- (1 - Sp[test_id_typ6[i6]]) * (1 - prev_i6[i6]) /
      ((1 - Sp[test_id_typ6[i6]]) * (1 - prev_i6[i6]) +
         Se[test_id_typ6[i6]] * prev_i6[i6])
    ### p(d-|T-)
    d0_T0_i6[i6] <- Sp[test_id_typ6[i6]] * (1 - prev_i6[i6]) /
      (Sp[test_id_typ6[i6]] * (1 - prev_i6[i6]) +
         (1 - Se[test_id_typ6[i6]]) * prev_i6[i6])

  ### herd level posterior probability of infection
  predicted_proba[herd_id_pr6[i6]] <- d1_T1_i6[i6] * p_T1_i6[i6] +
                         d1_T0_i6[i6] * (1 - p_T1_i6[i6]) +
                         D1_d0_i6[i6] * (
                         d0_T1_i6[i6] * p_T1_i6[i6] +
                         d0_T0_i6[i6] * (1 - p_T1_i6[i6])  )

  # status sampled from probability
  predicted_status[herd_id_pr6[i6]] ~ dbern(predicted_proba[herd_id_pr6[i6]])

  }'},
                 '
##############################################################################
###  Loop for monthly prevalences
##############################################################################
for(i_month in 1:month_max){

 month_prev[i_month] <- sum(Status[month_mat[1:month_N[i_month], i_month]]) /
                            month_N[i_month]

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
  logit_pi_within ~ dnorm(logit_pi_within_mean, logit_pi_within_prec)
  pi_within <- ilogit(logit_pi_within)

  ## Probability of not eliminating the infection
  logit_tau2 ~ dnorm(logit_tau2_mean, logit_tau2_prec)
  tau2 <- ilogit(logit_tau2)

  ## Logistic regression coefficients
  for(i_rf in 1:n_risk_factors){

  theta[i_rf] ~ dnorm(theta_norm_mean[i_rf], theta_norm_prec[i_rf])

 }

}')

}
