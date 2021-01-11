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

  ## extracting MCMC draws from the STOCfree_model
  samples <- x$mcmc

  ## tidying the results for model parameters
  n_tests <- length(grep("Se", colnames(samples[[1]])))
  herd_lev <- ifelse("pi_within" %in% colnames(samples[[1]]), 0, 1)

  if(herd_lev == 1 & n_tests == 1){

    param <- tidybayes::spread_draws(samples, Se, Sp, tau1, tau2)

  }

  if(herd_lev == 1 & n_tests > 1){

    param <- tidybayes::spread_draws(samples, Se[..], Sp[..], tau1, tau2)

  }

  if(herd_lev == 0 & n_tests == 1){

    param <- tidybayes::spread_draws(samples, Se, Sp, tau1, tau2, pi_within)

  }

  if(herd_lev == 0 & n_tests > 1){

    param <- tidybayes::spread_draws(samples, Se[..], Sp[..], tau1, tau2, pi_within)

  }

  param <- new_STOCfree_param(param)

  return(param)

  }


#' Importing MCMC samples from the results of a STOC free model
#'
#' @param out_path
#'
#' @return
#' @export
read_STOCfree_param <- function(out_path = "STOCfree_files"){

  nm    <- paste0(out_path, "/parameters.csv")
  param <- new_STOCfree_param(read.csv(nm))

  return(param)

  }


#' Title
#'
#' @param x
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


#' Title
#'
#' @param x
#' @param parameter
#' @param type
#'
#' @return
#' @export
plot.STOCfree_param <- function(x,
                                parameter = NULL,
                                type = "density"){

  ## list of parameters
  ls_param <- colnames(x)[4:length(x)]
  ## chacking if paramter in list of parameters
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

