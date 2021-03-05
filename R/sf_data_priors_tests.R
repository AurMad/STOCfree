#' Shows parameter of the Beta distributions used to represent priors on test charactersitics
#'
#' @param x a STOCfree_data object
#'
#' @return
#' @export
#'
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
set_priors_tests <- function(x = STOCfree_data(),
                 test = character(),
                 Se_a = 1, Se_b = 1,
                 Sp_a = 1, Sp_b = 1){

    if(length(test) == 0 & nrow(x$test_perf_prior) == 1){

    x$test_perf_prior$Se_a <- Se_a
    x$test_perf_prior$Se_b <- Se_b
    x$test_perf_prior$Sp_a <- Sp_a
    x$test_perf_prior$Sp_b <- Sp_b

    } else {

    row_id <- match(test, x$test_perf_prior$test)
    x$test_perf_prior$Se_a[row_id] <- Se_a
    x$test_perf_prior$Se_b[row_id] <- Se_b
    x$test_perf_prior$Sp_a[row_id] <- Sp_a
    x$test_perf_prior$Sp_b[row_id] <- Sp_b


  }

  x

}

#' Plot the distribution for the prior distributions on test characteristics
#'
#' @param data an object of class STOCfree_data
#' @param test name of at least one of the tests used. Use show_tests() to see the different tests. If omitted, the function will lot the prior distributions of all tests.
#'
#' @return
#' @export
#'
plot_priors_tests <- function(data = STOCfree_data(),
                              test = character()){

  priors <- data$test_perf_prior
  if(length(test) > 0){

    priors <- priors[match(test, priors$test),]

  }

  n_tests <- nrow(priors)

  par(mfrow = c(n_tests, 2))
  for(n in 1:n_tests){
  curve(dbeta(x, priors$Se_a[n], priors$Se_b[n]),
        from = 0, to = 1,
        main = priors$test[n],
        xlab = "Sensitivity",
        ylab = "Density")
  curve(dbeta(x, priors$Sp_a[n], priors$Sp_b[n]),
        from = 0, to = 1,
        main = priors$test[n],
        xlab = "Specificity",
        ylab = "Density")
  }

 }
