#################################################################################
### The code used for the Stan model is adapted from Damiano et al. (2017)
### See: https://github.com/luisdamiano/stancon18
### Here, because there are only 4 possible transitions between 2 states,
### the possible transitions are modelled explicitly
### also, for most consecutive time step, there will be no observation (emission)
### This is dealt with using an if statement whereby the latent state is only
### driven by the probability of transition
#################################################################################

#' Writes the code for the Stan model. For internal use.
#'
#' This function is for internal use.
#'
#' @return creates the Stan model
#'
write_Stan_model <- function(data){

  UseMethod("write_Stan_model")
}


write_Stan_model.default <- function(data){

  print("No method defined for this type of data")

}


## Stan model: herd level, no risk factor
write_Stan_model.herd <- function(data){
"
data{

  int<lower=1> n_herds;
  int<lower=1> herds_t1[n_herds];
  int<lower=1> herds_t2[n_herds];
  int<lower=1> herds_T[n_herds];
  int<lower=1> N;
  int<lower=0, upper=3> test_res[N];
  int<lower = 1> n_tests;
  int<lower = 0> test_id[N];
  real<lower = 0> Se_beta_a[n_tests];
  real<lower = 0> Se_beta_b[n_tests];
  real<lower = 0> Sp_beta_a[n_tests];
  real<lower = 0> Sp_beta_b[n_tests];
  real<lower = 0> pi1_beta_a;
  real<lower = 0> pi1_beta_b;
  real<lower = 0> tau1_beta_a;
  real<lower = 0> tau1_beta_b;
  real<lower = 0> tau2_beta_a;
  real<lower = 0> tau2_beta_b;
  int<lower = 1> n_month;
  int<lower = 1> month_N[n_month];
  int month_mat [n_herds, n_month];

}
parameters{

  real<lower = 0, upper = 1> Se[n_tests];
  real<lower = 0, upper = 1> Sp[n_tests];
  real<lower = 0, upper = 1> pi1;
  real<lower = 0, upper = 1> tau1;
  real<lower = 0, upper = 1> tau2;

}
transformed parameters{

  // logalpha needs to be accessible to toher blocks
  matrix[N, 2] logalpha;

  {

    // accumulator used at each time step
    real accumulator[2];

    // looping over all herds

    for(h in 1:n_herds){

      // first test in sequence
      // negative status
      logalpha[herds_t1[h], 1] = log(1 - pi1) + bernoulli_lpmf(test_res[herds_t1[h]] | 1 - Sp[test_id[herds_t1[h]]]);
      // positive status
      logalpha[herds_t1[h], 2] = log(pi1) + bernoulli_lpmf(test_res[herds_t1[h]] | Se[test_id[herds_t1[h]]]);

      // tests 2 in T in sequence
      for(t in herds_t2[h]:herds_T[h]){

        if(test_res[t] == 3){

          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(1 - tau1);
          // transition from status positive to status negative (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(1 - tau2);

          logalpha[t, 1] = log_sum_exp(accumulator);


          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(tau1);
          // transition from status positive to status positive (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(tau2);

          logalpha[t, 2] = log_sum_exp(accumulator);

        } else {

          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(1 - tau1) + bernoulli_lpmf(test_res[t] | 1 - Sp[test_id[t]]);
          // transition from status positive to status negative (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(1 - tau2) + bernoulli_lpmf(test_res[t] | 1 - Sp[test_id[t]]);

          logalpha[t, 1] = log_sum_exp(accumulator);

          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(tau1) + bernoulli_lpmf(test_res[t] | Se[test_id[t]]);
          // transition from status positive to status positive (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(tau2) + bernoulli_lpmf(test_res[t] | Se[test_id[t]]);

          logalpha[t, 2] = log_sum_exp(accumulator);

        } // if

      } // time sequence loop

    } // herd loop

  } //local

} // end of block
model{

  for(i_test in 1:n_tests){

     Se[i_test] ~ beta(Se_beta_a[i_test], Se_beta_b[i_test]);
     Sp[i_test] ~ beta(Sp_beta_a[i_test], Sp_beta_b[i_test]);

  }

  pi1 ~ beta(pi1_beta_a, pi1_beta_b);
  tau1 ~ beta(tau1_beta_a, tau1_beta_b);
  tau2 ~ beta(tau2_beta_a, tau2_beta_b);

  // update based only on last logalpha of each herd
  for(i in 1:n_herds)
    target += log_sum_exp(logalpha[herds_T[i]]);

}
generated quantities{

  // variable in which predictions are stored
  real pred[n_herds];

  // monthly prevalences
  real<lower = 0> month_prev[n_month];


  // loop in which the probabilities of infection are predicted
  {
    matrix[n_herds, 2] alpha;
  // predicted statuses on each month
    matrix[N, 2] proba_pred_0;
    real proba_pred[N];
    int status_pred[N];

  // loop in which the probabilities of infection are predicted
    for(i in 1:n_herds){
      alpha[i] = softmax(logalpha[herds_T[i],]')';
                         pred[i] = alpha[i, 2];
    }

   // predicted porbabilities of status positive for each ehrd on each month
    for(j in 1:N){

     proba_pred_0[j] = softmax(logalpha[j,]')';
     proba_pred[j] = proba_pred_0[j, 2];
     status_pred[j] = bernoulli_rng(proba_pred[j]);

    }

  // estimation of monthly prevalences
    for(i_month in 1:n_month){

    month_prev[i_month] = 0;

     for(ii_month in 1:month_N[i_month]){

     month_prev[i_month] = month_prev[i_month] + status_pred[month_mat[ii_month, i_month]];

     }

     month_prev[i_month] =  month_prev[i_month] / month_N[i_month];

  }

  }
}"
}




## Stan model: herd level, no risk factor
## Dynamics parameters modelled on the logit scale
write_Stan_model.herd_dynLogit <- function(data){
"
data{

  int<lower=1> n_herds;
  int<lower=1> herds_t1[n_herds];
  int<lower=1> herds_t2[n_herds];
  int<lower=1> herds_T[n_herds];
  int<lower=1> N;
  int<lower=0, upper=3> test_res[N];
  int<lower = 1> n_tests;
  int<lower = 0> test_id[N];
  real<lower = 0> Se_beta_a[n_tests];
  real<lower = 0> Se_beta_b[n_tests];
  real<lower = 0> Sp_beta_a[n_tests];
  real<lower = 0> Sp_beta_b[n_tests];
  real logit_pi1_mean;
  real logit_pi1_sd;
  real logit_tau1_mean;
  real logit_tau1_sd;
  real logit_tau2_mean;
  real logit_tau2_sd;
  int<lower = 1> n_month;
  int<lower = 1> month_N[n_month];
  int month_mat [n_herds, n_month];

}
parameters{

  real<lower = 0, upper = 1> Se[n_tests];
  real<lower = 0, upper = 1> Sp[n_tests];
  real<lower = 0, upper = 1> pi1;
  real<lower = 0, upper = 1> tau1;
  real<lower = 0, upper = 1> tau2;

}
transformed parameters{

  // logalpha needs to be accessible to toher blocks
  matrix[N, 2] logalpha;

  {

    // accumulator used at each time step
    real accumulator[2];

    // looping over all herds

    for(h in 1:n_herds){

      // first test in sequence
      // negative status
      logalpha[herds_t1[h], 1] = log(1 - pi1) + bernoulli_lpmf(test_res[herds_t1[h]] | 1 - Sp[test_id[herds_t1[h]]]);
      // positive status
      logalpha[herds_t1[h], 2] = log(pi1) + bernoulli_lpmf(test_res[herds_t1[h]] | Se[test_id[herds_t1[h]]]);

      // tests 2 in T in sequence
      for(t in herds_t2[h]:herds_T[h]){

        if(test_res[t] == 3){

          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(1 - tau1);
          // transition from status positive to status negative (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(1 - tau2);

          logalpha[t, 1] = log_sum_exp(accumulator);


          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(tau1);
          // transition from status positive to status positive (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(tau2);

          logalpha[t, 2] = log_sum_exp(accumulator);

        } else {

          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(1 - tau1) + bernoulli_lpmf(test_res[t] | 1 - Sp[test_id[t]]);
          // transition from status positive to status negative (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(1 - tau2) + bernoulli_lpmf(test_res[t] | 1 - Sp[test_id[t]]);

          logalpha[t, 1] = log_sum_exp(accumulator);

          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(tau1) + bernoulli_lpmf(test_res[t] | Se[test_id[t]]);
          // transition from status positive to status positive (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(tau2) + bernoulli_lpmf(test_res[t] | Se[test_id[t]]);

          logalpha[t, 2] = log_sum_exp(accumulator);

        } // if

      } // time sequence loop

    } // herd loop

  } //local

} // end of block
model{

  for(i_test in 1:n_tests){

     Se[i_test] ~ beta(Se_beta_a[i_test], Se_beta_b[i_test]);
     Sp[i_test] ~ beta(Sp_beta_a[i_test], Sp_beta_b[i_test]);

  }

  logit(pi1) ~ normal(logit_pi1_mean, logit_pi1_sd);
  logit(tau1) ~ normal(logit_tau1_mean, logit_tau1_sd);
  logit(tau2) ~ normal(logit_tau2_mean, logit_tau2_sd);

  // update based only on last logalpha of each herd
  for(i in 1:n_herds)
    target += log_sum_exp(logalpha[herds_T[i]]);

}
generated quantities{

  // variable in which predictions are stored
  real pred[n_herds];

  // monthly prevalences
  real<lower = 0> month_prev[n_month];


  // loop in which the probabilities of infection are predicted
  {
    matrix[n_herds, 2] alpha;
  // predicted statuses on each month
    matrix[N, 2] proba_pred_0;
    real proba_pred[N];
    int status_pred[N];

  // loop in which the probabilities of infection are predicted
    for(i in 1:n_herds){
      alpha[i] = softmax(logalpha[herds_T[i],]')';
                         pred[i] = alpha[i, 2];
    }

   // predicted porbabilities of status positive for each ehrd on each month
    for(j in 1:N){

     proba_pred_0[j] = softmax(logalpha[j,]')';
     proba_pred[j] = proba_pred_0[j, 2];
     status_pred[j] = bernoulli_rng(proba_pred[j]);

    }

  // estimation of monthly prevalences
    for(i_month in 1:n_month){

    month_prev[i_month] = 0;

     for(ii_month in 1:month_N[i_month]){

     month_prev[i_month] = month_prev[i_month] + status_pred[month_mat[ii_month, i_month]];

     }

     month_prev[i_month] =  month_prev[i_month] / month_N[i_month];

  }

  }
}"
}


## Stan model: herd level, with risk factors
write_Stan_model.herd_rf <- function(data){
  "
data{

  int<lower=1> n_herds;
  int<lower=1> herds_t1[n_herds];
  int<lower=1> herds_t2[n_herds];
  int<lower=1> herds_T[n_herds];
  int<lower=1> N;
  int<lower=0, upper=3> test_res[N];
  int<lower = 1> n_tests;
  int<lower = 0> test_id[N];
  real<lower = 0> Se_beta_a[n_tests];
  real<lower = 0> Se_beta_b[n_tests];
  real<lower = 0> Sp_beta_a[n_tests];
  real<lower = 0> Sp_beta_b[n_tests];
  real<lower = 0> pi1_beta_a;
  real<lower = 0> pi1_beta_b;
  real<lower = 0> tau2_beta_a;
  real<lower = 0> tau2_beta_b;
  int<lower = 0> n_risk_factors;
  real theta_norm_mean[n_risk_factors];
  real theta_norm_sd[n_risk_factors];
  matrix[N, n_risk_factors] risk_factors;
  int<lower = 1> n_month;
  int<lower = 1> month_N[n_month];
  int month_mat [n_herds, n_month];

}
parameters{

  real<lower = 0, upper = 1> Se[n_tests];
  real<lower = 0, upper = 1> Sp[n_tests];
  real<lower = 0, upper = 1> pi1;
  real<lower = 0, upper = 1> tau2;
  vector[n_risk_factors] theta;

}
transformed parameters{

  // logalpha needs to be accessible to toher blocks
  matrix[N, 2] logalpha;

  {

    // accumulator used at each time step
    real tau1[N];
    real accumulator[2];

    // logistic regression for tau1
  for(n in 1:N){

  tau1[n] = inv_logit(risk_factors[n,] * theta);

  }


    // looping over all herds
    for(h in 1:n_herds){

      // first test in sequence
      // negative status
      logalpha[herds_t1[h], 1] = log(1 - pi1) + bernoulli_lpmf(test_res[herds_t1[h]] | 1 - Sp[test_id[herds_t1[h]]]);
      // positive status
      logalpha[herds_t1[h], 2] = log(pi1) + bernoulli_lpmf(test_res[herds_t1[h]] | Se[test_id[herds_t1[h]]]);

      // tests 2 in T in sequence
      for(t in herds_t2[h]:herds_T[h]){

        if(test_res[t] == 3){

          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(1 - tau1[t]);
          // transition from status positive to status negative (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(1 - tau2);

          logalpha[t, 1] = log_sum_exp(accumulator);


          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(tau1[t]);
          // transition from status positive to status positive (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(tau2);

          logalpha[t, 2] = log_sum_exp(accumulator);

        } else {

          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(1 - tau1[t]) + bernoulli_lpmf(test_res[t] | 1 - Sp[test_id[t]]);
          // transition from status positive to status negative (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(1 - tau2) + bernoulli_lpmf(test_res[t] | 1 - Sp[test_id[t]]);

          logalpha[t, 1] = log_sum_exp(accumulator);

          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(tau1[t]) + bernoulli_lpmf(test_res[t] | Se[test_id[t]]);
          // transition from status positive to status positive (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(tau2) + bernoulli_lpmf(test_res[t] | Se[test_id[t]]);

          logalpha[t, 2] = log_sum_exp(accumulator);

        } // if

      } // time sequence loop

    } // herd loop

  } //local

} // end of block
model{


// priors for test characteristics
  for(i_test in 1:n_tests){

     Se[i_test] ~ beta(Se_beta_a[i_test], Se_beta_b[i_test]);
     Sp[i_test] ~ beta(Sp_beta_a[i_test], Sp_beta_b[i_test]);

  }

// priors for status dynamics
  pi1 ~ beta(pi1_beta_a, pi1_beta_b);
  tau2 ~ beta(tau2_beta_a, tau2_beta_b);

// priors for the logistic regression coefficients
  for(k in 1:n_risk_factors){

  theta[k] ~ normal(theta_norm_mean[k], theta_norm_sd[k]);

  }

// update based only on last logalpha of each herd
  for(i in 1:n_herds)
    target += log_sum_exp(logalpha[herds_T[i]]);

}
generated quantities{

  // variable in which predictions are stored
  real pred[n_herds];

  // monthly prevalences
  real<lower = 0> month_prev[n_month];


  {
    matrix[n_herds, 2] alpha;
  // predicted statuses on each month
    matrix[N, 2] proba_pred_0;
    real proba_pred[N];
    int status_pred[N];

  // loop in which the probabilities of infection are predicted
    for(i in 1:n_herds){
      alpha[i] = softmax(logalpha[herds_T[i],]')';
                         pred[i] = alpha[i, 2];
    }

   // predicted porbabilities of status positive for each ehrd on each month
    for(j in 1:N){

     proba_pred_0[j] = softmax(logalpha[j,]')';
     proba_pred[j] = proba_pred_0[j, 2];
     status_pred[j] = bernoulli_rng(proba_pred[j]);

    }

  // estimation of monthly prevalences
    for(i_month in 1:n_month){

    month_prev[i_month] = 0;

     for(ii_month in 1:month_N[i_month]){

     month_prev[i_month] = month_prev[i_month] + status_pred[month_mat[ii_month, i_month]];

     }

     month_prev[i_month] =  month_prev[i_month] / month_N[i_month];

  }

  }
}"
}



## Stan model: herd level, with risk factors
write_Stan_model.herd_dynLogit_rf <- function(data){
  "
data{

  int<lower=1> n_herds;
  int<lower=1> herds_t1[n_herds];
  int<lower=1> herds_t2[n_herds];
  int<lower=1> herds_T[n_herds];
  int<lower=1> N;
  int<lower=0, upper=3> test_res[N];
  int<lower = 1> n_tests;
  int<lower = 0> test_id[N];
  real<lower = 0> Se_beta_a[n_tests];
  real<lower = 0> Se_beta_b[n_tests];
  real<lower = 0> Sp_beta_a[n_tests];
  real<lower = 0> Sp_beta_b[n_tests];
  real logit_pi1_mean;
  real logit_pi1_sd;
  real logit_tau2_mean;
  real logit_tau2_sd;
  int<lower = 0> n_risk_factors;
  real theta_norm_mean[n_risk_factors];
  real theta_norm_sd[n_risk_factors];
  matrix[N, n_risk_factors] risk_factors;
  int<lower = 1> n_month;
  int<lower = 1> month_N[n_month];
  int month_mat [n_herds, n_month];

}
parameters{

  real<lower = 0, upper = 1> Se[n_tests];
  real<lower = 0, upper = 1> Sp[n_tests];
  real<lower = 0, upper = 1> pi1;
  real<lower = 0, upper = 1> tau2;
  vector[n_risk_factors] theta;

}
transformed parameters{

  // logalpha needs to be accessible to toher blocks
  matrix[N, 2] logalpha;

  {

    // accumulator used at each time step
    real tau1[N];
    real accumulator[2];

    // logistic regression for tau1
  for(n in 1:N){

  tau1[n] = inv_logit(risk_factors[n,] * theta);

  }


    // looping over all herds
    for(h in 1:n_herds){

      // first test in sequence
      // negative status
      logalpha[herds_t1[h], 1] = log(1 - pi1) + bernoulli_lpmf(test_res[herds_t1[h]] | 1 - Sp[test_id[herds_t1[h]]]);
      // positive status
      logalpha[herds_t1[h], 2] = log(pi1) + bernoulli_lpmf(test_res[herds_t1[h]] | Se[test_id[herds_t1[h]]]);

      // tests 2 in T in sequence
      for(t in herds_t2[h]:herds_T[h]){

        if(test_res[t] == 3){

          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(1 - tau1[t]);
          // transition from status positive to status negative (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(1 - tau2);

          logalpha[t, 1] = log_sum_exp(accumulator);


          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(tau1[t]);
          // transition from status positive to status positive (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(tau2);

          logalpha[t, 2] = log_sum_exp(accumulator);

        } else {

          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(1 - tau1[t]) + bernoulli_lpmf(test_res[t] | 1 - Sp[test_id[t]]);
          // transition from status positive to status negative (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(1 - tau2) + bernoulli_lpmf(test_res[t] | 1 - Sp[test_id[t]]);

          logalpha[t, 1] = log_sum_exp(accumulator);

          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(tau1[t]) + bernoulli_lpmf(test_res[t] | Se[test_id[t]]);
          // transition from status positive to status positive (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(tau2) + bernoulli_lpmf(test_res[t] | Se[test_id[t]]);

          logalpha[t, 2] = log_sum_exp(accumulator);

        } // if

      } // time sequence loop

    } // herd loop

  } //local

} // end of block
model{


// priors for test characteristics
  for(i_test in 1:n_tests){

     Se[i_test] ~ beta(Se_beta_a[i_test], Se_beta_b[i_test]);
     Sp[i_test] ~ beta(Sp_beta_a[i_test], Sp_beta_b[i_test]);

  }

// priors for status dynamics
  logit(pi1) ~ normal(logit_pi1_mean, logit_pi1_sd);
  logit(tau2) ~ normal(logit_tau2_mean, logit_tau2_sd);

// priors for the logistic regression coefficients
  for(k in 1:n_risk_factors){

  theta[k] ~ normal(theta_norm_mean[k], theta_norm_sd[k]);

  }

// update based only on last logalpha of each herd
  for(i in 1:n_herds)
    target += log_sum_exp(logalpha[herds_T[i]]);

}
generated quantities{

  // variable in which predictions are stored
  real pred[n_herds];

  // monthly prevalences
  real<lower = 0> month_prev[n_month];


  {
    matrix[n_herds, 2] alpha;
  // predicted statuses on each month
    matrix[N, 2] proba_pred_0;
    real proba_pred[N];
    int status_pred[N];

  // loop in which the probabilities of infection are predicted
    for(i in 1:n_herds){
      alpha[i] = softmax(logalpha[herds_T[i],]')';
                         pred[i] = alpha[i, 2];
    }

   // predicted porbabilities of status positive for each ehrd on each month
    for(j in 1:N){

     proba_pred_0[j] = softmax(logalpha[j,]')';
     proba_pred[j] = proba_pred_0[j, 2];
     status_pred[j] = bernoulli_rng(proba_pred[j]);

    }

  // estimation of monthly prevalences
    for(i_month in 1:n_month){

    month_prev[i_month] = 0;

     for(ii_month in 1:month_N[i_month]){

     month_prev[i_month] = month_prev[i_month] + status_pred[month_mat[ii_month, i_month]];

     }

     month_prev[i_month] =  month_prev[i_month] / month_N[i_month];

  }

  }
}"
}


