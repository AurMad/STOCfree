#' List of risk factors of new infection in a STOCfree_data
#'
#' @param x a STOC_free data object
#'
#' @return Displays a table with the list of risk factors considered as well as the associated priors
#' @export
#'
show_rf <- function(x = STOCfree_data()){

  if(attr(x, "number of risk factors") == 0) stop("No risk factor defined")
  x$risk_factors

}

#' Set priors associated with risk factors of new infection
#'
#' @param x an object of class STOCfree_data
#' @param risk_factor name of the risk factor for which priors are defined
#' @param modality name of the modality for which priors are defined. Only applies to categorical risk factors.
#' @param mean mean of the normal distribution used as a prior
#' @param sd  standard deviation of the normal distribution used as a prior
#'
#' @return
#' @export
#'
set_priors_rf <- function(x = STOCfree_data(),
                          risk_factor = character(),
                          modality = NULL,
                          mean = NULL,
                          sd = NULL){

  if(attr(x, "number of risk factors") == 0) stop("No risk factor defined")

  if(length(risk_factor) > 1 || (exists("modality") && length(modality) > 1) || length(mean) > 1 || length(sd) > 1){

    stop("Define priors for 1 variable / modality at a time")

  }


  ## Risk factor extracted from the STOCfree_data object
  risk_factors <- x$risk_factors

  ## creation of columns for parameters of prior distributions
  if(!"mean_prior" %in% colnames(risk_factors)){

    risk_factors$mean_prior <- rep(NA)
    risk_factors$sd_prior   <- rep(NA)

  }


  rf_rows <- which(risk_factors$risk_factor == risk_factor)
  if(risk_factor == "Intercept"){

    risk_factors$mean_prior[rf_rows] <- mean
    risk_factors$sd_prior[rf_rows] <- sd

  }


  if(unique(risk_factors$type[rf_rows]) == "continuous"){

    risk_factors$mean_prior[rf_rows] <- mean
    risk_factors$sd_prior[rf_rows] <- sd

  }

  if(unique(risk_factors$type[rf_rows]) == "categorical"){


    if(is.null(modality)){

      risk_factors$mean_prior[rf_rows] <- mean
      risk_factors$sd_prior[rf_rows] <- sd

    } else {

    if(! modality %in% risk_factors$modality[rf_rows]) stop("Modality not found")

    md_row <- which(risk_factors$modality == modality)

    row <- intersect(rf_rows, md_row)

    risk_factors$mean_prior[row] <- mean
    risk_factors$sd_prior[row] <- sd

       }

  }

  x$risk_factors <- risk_factors

  x

}



#' Plotting prior distributions for the risk factors of new infection
#'
#' @param x a dataset in the STOCfree_data format
#' @param n number of samples to draw from the prior distribution
#'
#' @return
#' @export
#'
plot_priors_rf <- function(x = STOCfree_data(), n = 100000){

  if(attr(x, "number of risk factors") == 0) stop("No risk factor defined")

  ## Risk factor data
  risk_factors <- x$risk_factors

  ## Number of variables (including intercept)
  lg_x <- nrow(risk_factors)
  ## Number of rows of the plotting window (assuming 2 columns)
  plot_nrows <- ceiling(lg_x / 2)

  dist_samples <- matrix(rep(NA, n * lg_x),
                         ncol = lg_x)
  ## sampling values from prior normal distributions for all variables
  for(i in 1:lg_x){

    dist_samples[, i] <- rnorm(n = n,
                               mean = risk_factors$mean_prior[i],
                               sd = risk_factors$sd_prior[i])

  }

  ## plots
  par(mfrow = c(plot_nrows, 2))

  ## Plot for intercept
  plot(density(invlogit(dist_samples[, 1])),
       main = "Intercept",
       xlab = "Probability")

  for(j in 2:lg_x){

    rf_j <- risk_factors$risk_factor[j]

    if(risk_factors$type[j] == "continuous"){

      range_x <- range(x$risk_factor_data[, rf_j])
      seq_x   <- matrix(
        seq(range_x[1], range_x[2], length.out = 100),
        nrow = 1)

      y_val <- dist_samples[, j] %*% seq_x

      for(k in 1:ncol(y_val)){

        y_val[, k] <- invlogit(dist_samples[, 1] + y_val[, k])

      }

      to_plot <- data.frame(
        x = as.double(seq_x),
        lb  = unlist(apply(y_val, 2, quantile, .025, na.rm = TRUE)),
        mid = unlist(apply(y_val, 2, quantile, .5, na.rm = TRUE)),
        ub  = unlist(apply(y_val, 2, quantile, .975, na.rm = TRUE))
      )

      rm(y_val)

      plot(mid ~ x, data = to_plot,
           type = "l", lwd = 2,
           ylim = c(0, 1), ylab = "Probability",
           xlab = rf_j, main = rf_j)
      lines(lb ~ x, data = to_plot, lty = 2)
      lines(ub ~ x, data = to_plot, lty = 2)

    } else {

      dist_samples[, j] <- invlogit(dist_samples[, 1] + dist_samples[, j])

      plot(density(dist_samples[, j]),
           main = paste("Risk factor", rf_j, "- modality ", risk_factors$modality[j]),
           xlab = "Probability")

    }
  }
}

