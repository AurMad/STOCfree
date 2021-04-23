#' logit function
#'
#' @param p a proportion
#'
#' @return the logit of a proportion
#' @export
#'
#' @examples logit(.5)
#'
logit <- function(p){

  if(length(p[p < 0]) > 0 | length(p[p > 1]) > 0) stop("p should be between 0 and 1")

  log(p / (1 - p))

}

#' Inverse logit function
#'
#' @param x a number
#'
#' @return a proportion
#' @export
#'
#' @examples
#' invlogit(0)
#' invlogit(logit(.5)) == .5
invlogit <- function(x){

  exp(x) / (1 + exp(x))

}

#' The normal distribution on the logit scale
#'
#' @param x vector of quantiles on the probability scale
#' @param mean_logit vector of means on the logit scale
#' @param sd_logit vector of standard deviations on the logit scale
#'
#' @return
#' @export
#'
#' @examples
#' curve(dnorm_logit(x, 0, 1))
dnorm_logit <- function(x, mean_logit = .5, sd_logit = 1){

  z <- 1 / (sd_logit *sqrt(2 * pi)) * 1 / (x * (1 - x)) * exp(-(logit(x) - mean_logit)^2 / (2 * sd_logit^2))

  return(z)

}
