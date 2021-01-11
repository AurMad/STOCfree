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

  if(p < 0 | p > 1) stop("p should be between 0 and 1")

  log(p / (1 - p))

}

#' Inverse logit function
#'
#' @param x a number
#'
#' @return a proportion
#' @export
#'
#' @examples invlogit(logit(.5))
invlogit <- function(x){

  exp(x) / (1 + exp(x))

  }
