#' Shows the parameters asssociated with priors related to infection dynamics
#'
#' @param x a STOC_free data object
#'
#' @return a table with the alpha and beta parameters of the Beta distributions associated with the prior probabilities of infection on the first test and the probability of remaining infected between consecutive tests
#' @export
#'
#' @examples
show_inf_dyn <- function(x = STOCfree_data()){

  x$inf_dyn_priors

}


#' Sets the parameters asssociated with priors related to infection dynamics
#'
#' @param data a STOCfree_data object
#' @param pi1_a alpha parameter for the prior distribution of the probability of being infected on the first test
#' @param pi1_b beta parameter for the prior distribution of the probability of being infected on the first test
#' @param tau2_a alpha parameter for the prior distribution of the probability of remaining infected between consecutive months
#' @param tau2_b beta parameter for the prior distribution of the probability of remaining infected between consecutive months
#'
#' @return
#' @export
#'
#' @examples
set_priors_inf_dyn <- function(data = STOCfree_data(),
                               pi1_a = 1, pi1_b = 1,
                               tau1_a = 1, tau1_b = 1,
                               tau2_a = 1, tau2_b = 1){

  data$inf_dyn_priors["pi1_a"] = pi1_a
  data$inf_dyn_priors["pi1_b"] = pi1_b
  data$inf_dyn_priors["tau1_a"] = tau1_a
  data$inf_dyn_priors["tau1_b"] = tau1_b
  data$inf_dyn_priors["tau2_a"] = tau2_a
  data$inf_dyn_priors["tau2_b"] = tau2_b

  data

}

#' Plots the prior distributions for the parameters related to infection dynamics
#'
#' @param data a STOCfree_data object
#'
#' @return
#' @export
#'
#' @examples
plot_priors_inf_dyn <- function(data = STOCfree_data()){

  priors <- data$inf_dyn_priors
  n_risk_factors <- attributes(data, "number of risk factors")
  plot_n_rows <- ifelse(n_risk_factors == 0, 2, 1)

  par(mfrow = c(plot_n_rows, 2))
  curve(dbeta(x, priors["pi1_a"], priors["pi1_b"]),
        from = 0, to = 1,
        main = "Probability of infection\non the first test",
        xlab = "Probability",
        ylab = "Density")
  if(n_risk_factors == 0){

    curve(dbeta(x, priors["tau1_a"], priors["tau1_b"]),
          from = 0, to = 1,
          main = "Probability of new\ninfection between to 2 months",
          xlab = "tau1",
          ylab = "Density")

  }
  curve(dbeta(x, priors["tau2_a"], priors["tau2_b"]),
        from = 0, to = 1,
        main = "Probability of not eliminating the\ninfection between to 2 months",
        xlab = "tau2",
        ylab = "Density")

 }
