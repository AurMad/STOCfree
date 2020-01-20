#' Shows parameter of the Beta distributions used to represent priors on test charactersitics
#'
#' @param x a STOCfree_data object
#'
#' @return
#' @export
#'
#' @examples
show_tests <- function(x = STOCfree_data()){

  x$test_perf_prior

}


#' Set parameters for the Beta distributions used to represent priors on test charactersitics
#'
#' @param x a STOCfree_data object
#' @param Se_a alpha parameter for the prior on sensitivity
#' @param Se_b beta parameter for the prior on sensitivity
#' @param Sp_a alpha parameter for the prior on specificity
#' @param Sp_b beta parameter for the prior on specificity
#'
#' @return
#' @export
#'
#' @examples
set_priors_tests <- function(x = STOCfree_data(),
                 Se_a = 1, Se_b = 1,
                 Sp_a = 1, Sp_b = 1){

  x$test_perf_prior$Se_a = Se_a
  x$test_perf_prior$Se_b = Se_b
  x$test_perf_prior$Sp_a = Sp_a
  x$test_perf_prior$Sp_b = Sp_b

  x

}

#' Plot the distribution for the prior distributions on test characteristics
#'
#' @param data an object of class STOCfree_data
#'
#' @return
#' @export
#'
#' @examples
plot_priors_tests <- function(data = STOCfree_data()){

  priors <- data$test_perf_prior

  par(mfrow = c(1, 2))
  curve(dbeta(x, priors$Se_a, priors$Se_b),
        from = 0, to = 1,
        main = "Sensitivity",
        xlab = "Sensitivity",
        ylab = "Density")
  curve(dbeta(x, priors$Sp_a, priors$Sp_b),
        from = 0, to = 1,
        main = "Specificity",
        xlab = "Specificity",
        ylab = "Density")

 }
