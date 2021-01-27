## class and methods for parameters estimated by the STOCfree_model

validate_STOCfree_pred <- function(x){

  ## input should be a data.frame
  if(!is.data.frame(x)) stop("x should be a data.frame")

  ## checking column names
  if(!all.equal(match(c(".chain", ".iteration", ".draw", "predicted_proba"), colnames(x)), 2:5)) stop("wrong type of input")

  }


new_STOCfree_pred <- function(x){

  validate_STOCfree_pred(x)

  class(x) <- c("STOCfree_pred", class(x))
  attr(x,"number of herds") <- length(unique(x$herd))

  return(x)

  }

#' Title
#'
#' @param x
#'
#' @return
#' @export
extract_STOCfree_pred <- function(x){

  samples <- x$mcmc
  herd_id_corresp <- x$herd_id_corresp

  ## tidying the results for predicted probabilities
  predictions <- tidybayes::spread_draws(samples,
                                         predicted_proba[herd_id])

  predictions <- predictions[, c("herd_id", ".chain", ".iteration", ".draw", "predicted_proba")]

  ## adding original herd ids
  predictions <- merge(predictions, herd_id_corresp, all.x = TRUE)
  predictions <- predictions[, c(6, 2:5)]

  colnames(predictions)[1] <- "herd"

  predictions <- new_STOCfree_pred(predictions)

  return(predictions)

  }


#' Title
#'
#' @param out_path
#'
#' @return
#' @export
read_STOCfree_pred <- function(out_path = "STOCfree_files"){

  nm   <- paste0(out_path, "/predictions.csv")
  pred <- new_STOCfree_pred(read.csv(nm))

  return(pred)

  }


#' Title
#'
#' @param x
#'
#' @return
#' @export
print.STOCfree_pred <- function(x){

  cat("MCMC samples from STOC free model herd level predicted probabilities of infection\n\n")
  cat("Number of herds:", attr(x,"number of herds"), "\n\n")
  cat("Herds:", unique(x[,1]))

  }

#' Title
#'
#' @param x
#' @param herd
#' @param type
#'
#' @return
#' @export
plot.STOCfree_pred <- function(x, herd = "all", type = "aggregated", legend = FALSE){

  ## all herds plotted using a single line
  if(length(herd) == 1 && herd == "all" && type == "aggregated"){

    p <- ggplot2::ggplot(x, ggplot2::aes(x = predicted_proba))

  }

  ## all herds plotted with one line for each herd
  if(length(herd) == 1 && herd == "all" && type == "individual"){

    p <- ggplot2::ggplot(x, ggplot2::aes(x = predicted_proba, color = factor(herd)))

  }

  ## one or a subset of herds plotted
  if(length(herd) > 1 || herd != "all"){

    list_herds <- unique(x$herd)

    ## Are there requested herds that are missing?
    missing_herds <- which(!herd %in% list_herds)
    if(length(missing_herds) > 0) stop(cat("No predicted probabilities of infection for the following herd(s): ", missing_herds))

    if(type == "aggregated"){

      p <- ggplot2::ggplot(x[x$herd %in% herd,], ggplot2::aes(x = predicted_proba))

    }

    if(type == "individual"){

      p <- ggplot2::ggplot(x[x$herd %in% herd,], ggplot2::aes(x = predicted_proba, color = factor(herd)))

    }

  }

  ## if legend included, legend title is changed
  if(legend == TRUE){

    p <- p + ggplot2::scale_color_discrete(name = "Herd")

  }

  ## legend = FALSE
  if(legend == FALSE){

    p <- p + ggplot2::theme(legend.position = "none")

  }

  p +
    ggplot2::geom_density() +
    ggplot2::xlab("Predicted probability of infection")

  }
