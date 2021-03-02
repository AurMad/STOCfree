## class and methods for parameters estimated by the STOCfree_model

validate_STOCfree_param <- function(x){

  ## input should be a data.frame
  if(!is.data.frame(x)) stop("x should be a data.frame")

  ## columns associated with mcmc draws should be present
  columns_x <- colnames(x)
  mcmc_cols <- match(c(".chain", ".iteration", ".draw"), columns_x)

  if(length(mcmc_cols) != 3) stop("wrong type of input")

  ## names of parameters
  param_names <- grep("Se|Sp|tau|theta", colnames(x), value = TRUE)

  if(length(param_names) < 4) stop("wrong type of input")

  }


new_STOCfree_param <- function(x){

  validate_STOCfree_param(x)

  class(x) <- c("STOCfree_param", class(x))
  attr(x,"parameters") <- grep("Se|Sp|tau|theta|pi_within", colnames(x), value = TRUE)

  return(x)

  }


#' Extracting parameter MCMC samples from the results of a STOC free model
#'
#' @param x results returned by the STOC free model
#'
#' @return
#' @export
extract_STOCfree_param <- function(x){

  if("CmdStanMCMC" %in% class(x)){ # if Stan model

    ### extracting samples for parameters
    param_list <- c("Se", "Sp", "pi1", "tau1", "tau2")

    parameters <- x$draws(param_list)
    parameters <- posterior::as_draws_df(parameters)

    ## reordering columns for consistency with JAGS output
    id_cols1 <- match(c(".chain", ".iteration", ".draw"),
                      colnames(parameters))
    id_cols2 <- (1:length(parameters))[-id_cols1]

    parameters <- parameters[, c(id_cols1, id_cols2)]

  } else { # if not Stan model

  ## extracting MCMC draws from the STOCfree_model
  samples <- x$mcmc

  ## tidying the results for model parameters
  n_tests  <- length(grep("Se", colnames(samples[[1]])))
  herd_lev <- ifelse("pi_within" %in% colnames(samples[[1]]), 0, 1)
  rf       <- ifelse(length(grep("theta", colnames(samples[[1]]))) > 0, 1, 0)

  ## no risk factors
  if(herd_lev == 1 & n_tests == 1 & rf == 0){

    parameters <- tidybayes::spread_draws(samples, Se, Sp, tau1, tau2)

  }

  if(herd_lev == 1 & n_tests > 1 & rf == 0){

    parameters <- tidybayes::spread_draws(samples, Se[..], Sp[..], tau1, tau2)

  }

  if(herd_lev == 0 & n_tests == 1 & rf == 0){

    parameters <- tidybayes::spread_draws(samples, Se, Sp, tau1, tau2, pi_within)

  }

  if(herd_lev == 0 & n_tests > 1 & rf == 0){

    parameters <- tidybayes::spread_draws(samples, Se[..], Sp[..], tau1, tau2, pi_within)

  }

  ## with risk factors
  if(herd_lev == 1 & n_tests == 1 & rf == 1){

    parameters <- tidybayes::spread_draws(samples, Se, Sp, theta[..], tau2)

  }

  if(herd_lev == 1 & n_tests > 1 & rf == 1){

    parameters <- tidybayes::spread_draws(samples, Se[..], Sp[..], theta[..], tau2)

  }

  if(herd_lev == 0 & n_tests == 1 & rf == 1){

    parameters <- tidybayes::spread_draws(samples, Se, Sp, theta[..], tau2, pi_within)

  }

  if(herd_lev == 0 & n_tests > 1 & rf == 1){

    parameters <- tidybayes::spread_draws(samples, Se[..], Sp[..], theta[..], tau2, pi_within)

  }

  } # end if JAGS model

  parameters <- new_STOCfree_param(parameters)

  return(parameters)

  }


#' Importing MCMC samples for model parameters from the results of a STOC free model
#'
#' @param out_path folder in which the model, data and result files will be stored
#'
#' @return
#' @export
read_STOCfree_param <- function(out_path = "STOCfree_files"){

  nm    <- paste0(out_path, "/parameters.csv")
  param <- new_STOCfree_param(read.csv(nm))

  return(param)

  }


#' print method for STOCfree parameters
#'
#' @param x an object of class STOCfree_param created with extract_STOCfree_param() or read_STOCfree_param()
#'
#' @return
#' @export
print.STOCfree_param <- function(x){

  cat("MCMC samples from STOC free model parameters\n\n")
  cat(paste("Parameters:", paste(attr(x,"parameters"), collapse = ", "), "\n"))
  cat(paste("Number of chains:", max(x$.chain), "\n"))
  cat(paste("Number of iterations per chain:", max(x$.iteration), "\n"))
  cat(paste("Number of draws:", max(x$.iteration) * max(x$.chain), "\n"))

  }


#' plot method for STOCfree parameters
#'
#' @param x an object of class STOCfree_param created with extract_STOCfree_param() or read_STOCfree_param()
#' @param parameter parameter name
#' @param type type of plot. Either 'traceplot' or 'density'
#'
#' @return
#' @export
plot.STOCfree_param <- function(x,
                                parameter = NULL,
                                type = "density"){

  ## list of parameters
  ls_param <- colnames(x)[4:length(x)]
  ## checking if parameter in list of parameters
  if(is.null(parameter) | ! parameter %in% ls_param){

    stop(paste("Please provide the name of a parameter to plot from: ",
               paste0(ls_param, collapse = ", ")))

  }

  ## checking plot type
  if(!type %in% c("density", "traceplot") | length(type) > 1){

    stop("'type' argument is missing. 'type' should be either 'density' or 'traceplot'")

  }

  ## density plot
  if(type == "density"){

    sfp <- ggplot2::ggplot(x, ggplot2::aes_string(x = parameter)) +
      ggplot2::geom_density()

  }

  ## traceplot
  if(type == "traceplot"){

    x$.chain <- as.factor(x$.chain)

    sfp <- ggplot2::ggplot(x, ggplot2::aes_string(x = ".iteration", y = parameter, col = ".chain")) +
      ggplot2::geom_line() +
      ggplot2::xlab("Iteration") +
      ggplot2::scale_color_discrete(name = "Chain")

  }

  sfp

}


#' Effective sample sizes of STOC free model parameters
#'
#' @param x estimated model parameters
#'
#' @return
#' @export
ess_STOCfree_param <- function(x){

  ## checking class
  if(! "STOCfree_param" %in% class(x)) stop("Input should be of class 'STOCfree_param'")

  ## number of chains
  n_chains <- max(x$.chain)

  ## creating empty list
  mcmc <- list()

  ## populating object
  for(i in 1:n_chains){

    mcmc[[i]] <- as.matrix(x[x$.chain == i, -(1:3)])

  }

  ## adding class to mcmc object in order to be able to pass it to the effectiveSize function
  class(mcmc) <- "mcmc.list"

  ## effective sample sizes calculated using a function from the coda package
  ess <- coda::effectiveSize(mcmc)

  return(ess)

}


## function that calculates summary statistics
f_summary <- function(x, quantiles){

  z <- c(mean = mean(x),
         sd = sd(x),
         median = median(x),
         quantile(x, quantiles))

  return(z)

}


#' summary method for the STOCfree_param class
#'
#' @param x an object of class STOCfree_param created with extract_STOCfree_param() or read_STOCfree_param()
#' @param quantiles quantiles of the posterior distributions to display. Default values are 0.025 and 0.975
#'
#' @return
#' @export
summary.STOCfree_param <- function(x, quantiles = c(.025, .975)){

  z <- t(apply(x[,-(1:3)], 2, f_summary, quantiles = quantiles))

  ess <- ess_STOCfree_param(x)

  z <- as.data.frame(
    cbind(z, ess = ess[match(names(ess), rownames(z))]))

  return(z)

  }


