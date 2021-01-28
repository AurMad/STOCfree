# Function that generates output in the tidy format
# from run.jags output
STOCfree_tidy_output <- function(STOCfree_model_output,
                                 STOCfree_data){

  samples <- STOCfree_model_output$mcmc

  ## tidying the results for predicted probabilities
  predictions <- tidybayes::spread_draws(samples,
                                         predicted_proba[herd_id])

  predictions <- predictions[, c("herd_id", ".chain", ".iteration", ".draw", "predicted_proba")]

  if(!is.null(STOCfree_data)){

    predictions <- merge(predictions, STOCfree_data$herd_id_corresp, all.x = TRUE)
    predictions <- predictions[, c(6, 2:5)]

  }

  ## tidying the results for model parameters
  n_tests  <- length(grep("Se", colnames(samples[[1]])))
  herd_lev <- ifelse("pi_within" %in% colnames(samples[[1]]), 0, 1)
  rf       <- ifelse(length(grep("theta", colnames(samples[[1]]))) > 0, 1, 0)

  print(colnames(samples[[1]]))
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

  ## tidying the results for monthly prevalences
  ## model output made tidy
  month_prev <- tidybayes::spread_draws(samples,
                                        month_prev[..])
  ## reconstructing the list of all months
  months_all <- data.frame(
    date = seq(as.Date(paste0(attr(STOCfree_data, "month first test"), "-01")),
               as.Date(paste0(attr(STOCfree_data, "month last test"), "-01")),
               by = "1 month"),
    month_abb = NA,
    tidybayes_names = NA
  )

  months_all <- months_all[-nrow(months_all), ]
  ## month names abbreviated to get meaningful column names
  months_all$month_abb <- format(months_all$date, "%b_%Y")

  ## month names in the month_prev dataset
  months_all$tidybayes_names <- paste0("month_prev.", 1:nrow(months_all))

  ## columns reordered
  cln <- match(months_all$tidybayes_names, colnames(month_prev))
  month_prev <- month_prev[, c(1:3, cln)]

  ## column names changed
  cln1 <- match(months_all$tidybayes_names, colnames(month_prev))
  colnames(month_prev)[cln1] <- months_all$month_abb

  list(
    parameters = parameters,
    predictions = predictions,
    month_prev = month_prev)

}



## Folder where models are saved
STOCfree_files <- function(out_path = out_path){

  out_path <- paste0("./", out_path)

  ## Does the folder already exists?
  ls_dirs <- list.dirs()

  if(!out_path %in% ls_dirs){

    dir.create(paste0(getwd(), out_path))

  }

  return(out_path)

}




