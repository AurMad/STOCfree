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

#' Extracting MCMC samples for the predicted probabilities of being latent status positive from the results of a STOC free model
#'
#' @param x results returned by the STOC free model
#'
#' @return
#' @export
extract_STOCfree_pred <- function(x, STOCfree_data = NULL){

  if("CmdStanMCMC" %in% class(x)){

    ## extracting predicted probabilities from Stan output
    predictions <- x$draws("pred")
    ## formatting the output as a data.frame
    predictions <- posterior::as_draws_df(predictions)
    ## columns are reordered to match JAGS output
    id_cols1 <- match(c(".chain", ".iteration", ".draw"),
                      colnames(predictions))

    id_cols2 <- (1:length(predictions))[-id_cols1]

    predictions <- predictions[, c(id_cols1, id_cols2)]

    ## data put in long format
    predictions <- tidybayes::gather_variables(predictions)

    ## columns renamed for consistency with JAGS output
    colnames(predictions)[match(c(".variable", ".value"), colnames(predictions))] <- c("herd", "predicted_proba")

    ## herd ids extracted
    predictions$herd <- gsub("pred\\[", "", predictions$herd)
    predictions$herd <- gsub("\\]", "", predictions$herd)

    if(is.null(STOCfree_data)){

      warning("No STOCfree_data object supplied. Internal herd ids will be used.")

    } else {

      herd_id_corresp <- STOCfree_data$herd_id_corresp

      predictions$herd <- as.integer(predictions$herd)

      predictions$herd <- herd_id_corresp[predictions$herd, 1]

    }

    ## columns reordered
    predictions <- as.data.frame(predictions[, c("herd", ".chain", ".iteration", ".draw", "predicted_proba")])

  } else {

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

  }

  predictions <- new_STOCfree_pred(predictions)
  return(predictions)

  }


#' Importing MCMC samples for the predicted probabilities of being latent status positive
#'
#' @param out_path path to the folder where the results are stored
#'
#' @return
#' @export
read_STOCfree_pred <- function(out_path = "STOCfree_files"){

  nm   <- paste0(out_path, "/predictions.csv")
  pred <- new_STOCfree_pred(read.csv(nm))

  return(pred)

  }


#' print method for predicted probabilities of status positive
#'
#' @param x an object of class STOCfree_pred
#'
#' @return
#' @export
print.STOCfree_pred <- function(x){

  cat("MCMC samples from STOC free model herd level predicted probabilities of infection\n\n")
  cat("Number of herds:", attr(x,"number of herds"), "\n\n")
  cat("Herds:", unique(x[,1]))

  }

#' plot method for predicted probabilities of status positive
#'
#' @param x an object of class STOCfree_pred
#' @param herd either one or several herd ids or 'all'. all' will display the probabilities for all herds
#' @param type if 'aggregated' a single density plot will be displayed for all herds. If 'individual', lines of different colors are displayed for the different herds.
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