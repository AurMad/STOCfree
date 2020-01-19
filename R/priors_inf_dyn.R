#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
show_inf_dyn <- function(x = STOCfree_data()){

  x$inf_dyn_priors

}


#' Title
#'
#' @param data
#' @param pi1_a
#' @param pi1_b
#' @param tau2_a
#' @param tau2_b
#'
#' @return
#' @export
#'
#' @examples
set_priors_inf_dyn <- function(data = STOCfree_data(),
                               pi1_a = 1, pi1_b = 1,
                               tau2_a = 1, tau2_b = 1){

  data$inf_dyn_priors["pi1_a"] = pi1_a
  data$inf_dyn_priors["pi1_b"] = pi1_b
  data$inf_dyn_priors["tau2_a"] = tau2_a
  data$inf_dyn_priors["tau2_b"] = tau2_b

  data

}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
plot_priors_inf_dyn <- function(data = STOCfree_data()){

  priors <- data$inf_dyn_priors

  par(mfrow = c(1, 2))
  curve(dbeta(x, priors["pi1_a"], priors["pi1_b"]),
        from = 0, to = 1,
        main = "Probability of infection\non the first test",
        xlab = "Probability",
        ylab = "Density")
  curve(dbeta(x, priors["tau2_a"], priors["tau2_b"]),
        from = 0, to = 1,
        main = "Probability of not eliminating the\ninfection between to 2 months",
        xlab = "tau2",
        ylab = "Density")

 }
