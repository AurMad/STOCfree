#' Plot prior densities for the probability of new infection in the herds unexposed and exposed to the risk factors included in the STOCfree model
#'
#' @param x a list with the first element containing the mean(s) and the second element containing standard deviation(s) of normal distributions for the logistic regression
#' @param n number of draws from the parameter distributions used for plotting
#'
#' @return a plot of the prior distributions for the probabilities of new infection in the unexposed (intercept) and exposed to the risk factors evaluated
#' @export
#'
#' @examples
plot_priors_rf <- function(x, n = 100000){

  ## Number of samples to draw
  n <- 100000
  ## Number of variables (including intercept)
  lg_x <- length(x[[1]])
  ## Number of rows of the plotting window (assuming 2 columns
  plot_nrows <- ceiling(lg_x / 2)


  dist_samples <- matrix(rep(NA, n * lg_x),
                         ncol = lg_x)
  ## sampling values from prior normal dsitributions for all variables
  for(i in 1:lg_x){

    dist_samples[, i] <- rnorm(n = n,
                               mean = x[[1]][i],
                               sd = x[[2]][i])

  }

  ## Intercept sample value added to risk factor value row wise
  for(j in 2:lg_x){

    dist_samples[, j] <- invlogit(dist_samples[, 1] + dist_samples[, j])

  }

  ## Inverse logit of the sum to get a probability
  dist_samples[, 1] <- invlogit(dist_samples[, 1])

  par(mfrow = c(plot_nrows, 2))
  plot(density(dist_samples[, 1]),
       main = "Reference category",
       xlab = "Probability")

  for(k in 2:lg_x){

    plot(density(dist_samples[, k]),
         main = paste("Category exposed to\nrisk factor", k - 1),
         xlab = "Probability")


  }

}
