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

#' plot method for herd-level predicted probabilities of being status positive
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



#' summary method for herd-level predicted probabilities of being status positive
#'
#' @param x a STOCfree_pred object
#' @param herd a vector containing the herds for which a summary statistics a required
#' @param quantiles quantiles to be passed to the quantile() function
#' @param digits number of digits used in the object returned
#'
#' @return
#' @export
summary.STOCfree_pred <- function(x, herd = "all", quantiles = c(.025, .975), digits = 3){

  ## checking that quantiles are between 0 and 1
  rng_qtls <- range(quantiles)
  if(min(quantiles) < 0 | max(quantiles) > 1) stop("Quantiles should be between 0 and 1")


  ## selection of herds when only a subset is used
  if(length(herd) > 1 ||  herd != "all"){

    ## checking that requested herds are in the STOCfree_pred data
    herd_pos <- match(herd, x$herd)
    n_miss_herds <- length(herd_pos[is.na(herd_pos)])
    if(n_miss_herds > 0) stop(paste0("The following herds are missing from ", deparse(substitute(x)), ": ", herd[which(is.na(herd_pos))]))


    x <- x[x$herd %in% herd, ]

  } else { # if all herds, vector of unique herds

    herd <- unique(x$herd)

  }

  ## apply summary function to each herd
  means   <- tapply(x$pred, x$herd, mean)
  means   <- means[match(herd, names(means))]

  sds     <- tapply(x$pred, x$herd, sd)
  sds     <- sds[match(herd, names(sds))]

  medians <- tapply(x$pred, x$herd, median)
  medians <- medians[match(herd, names(medians))]

  qtls <- tapply(x$pred, x$herd, function(y) as.data.frame(quantile(y, quantiles)))
  qtls <- qtls[match(herd, names(qtls))]
  qtls <- do.call("rbind", qtls)
  colnames(qtls) <- paste0(quantiles * 100, "%")

  ## different summaries gathered in a single dataset
  z <- data.frame(
    herd = herd,
    mean = as.double(means),
    sd = as.double(sds),
    median = as.double(medians)
  )

  ## adding quantiles
  z <- as.data.frame(cbind(z, qtls))

  z[, 2:length(z)] <- round(z[, 2:length(z)], digits)

  rownames(z) <- 1:nrow(z)

  return(z)

}

