Stan_code <- function(){
"
// code adapted from Damiano 2017
data{

  int<lower=0> n_herds;
  int<lower=0> herds_t1[n_herds];
  int<lower=0> herds_t2[n_herds];
  int<lower=0> herds_T[n_herds];
  int<lower=0> N;
  int<lower=0, upper=3> test_res[N];
  real<lower = 0> Se_beta_a;
  real<lower = 0> Se_beta_b;
  real<lower = 0> Sp_beta_a;
  real<lower = 0> Sp_beta_b;
  real<lower = 0> pi1_beta_a;
  real<lower = 0> pi1_beta_b;
  real<lower = 0> tau1_beta_a;
  real<lower = 0> tau1_beta_b;
  real<lower = 0> tau2_beta_a;
  real<lower = 0> tau2_beta_b;


}
parameters{

  real<lower = 0, upper = 1> Se;
  real<lower = 0, upper = 1> Sp;
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
      logalpha[herds_t1[h], 1] = log(1 - pi1) + bernoulli_lpmf(test_res[herds_t1[h]] | 1 - Sp);
      // positive status
      logalpha[herds_t1[h], 2] = log(pi1) + bernoulli_lpmf(test_res[herds_t1[h]] | Se);

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
          accumulator[1] = logalpha[t-1, 1] + log(1 - tau1) + bernoulli_lpmf(test_res[t] | 1 - Sp);
          // transition from status positive to status negative (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(1 - tau2) + bernoulli_lpmf(test_res[t] | 1 - Sp);

          logalpha[t, 1] = log_sum_exp(accumulator);

          // transition from status negative to status negative (j = 1; i = 1)
          accumulator[1] = logalpha[t-1, 1] + log(tau1) + bernoulli_lpmf(test_res[t] | Se);
          // transition from status positive to status positive (j = 1; i = 1)
          accumulator[2] = logalpha[t-1, 2] + log(tau2) + bernoulli_lpmf(test_res[t] | Se);

          logalpha[t, 2] = log_sum_exp(accumulator);

        } // if

      } // time sequence loop

    } // herd loop

  } //local

} // end of block
model{

  Se ~ beta(Se_beta_a, Se_beta_b);
  Sp ~ beta(Sp_beta_a, Sp_beta_b);
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

  // loop in which the probabilities of infection are predicted
  // could be done in R
  {
    matrix[n_herds, 2] alpha;

    for(i in 1:n_herds){
      alpha[i] = softmax(logalpha[herds_T[i],]')';
                         pred[i] = alpha[i, 2];
    }
  }
}"
}
