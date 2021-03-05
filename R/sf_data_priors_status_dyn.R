#' Shows the parameters associated with priors related to infection dynamics
#'
#' @param x a STOC_free data object
#'
#' @return a table with the alpha and beta parameters of the Beta distributions associated with the prior probabilities of infection on the first test and the probability of remaining infected between consecutive tests
#' @export
show_priors_status_dyn <- function(data){

  ## checking that data is a STOCfree_data object
  data_nm <- deparse(substitute(data))
  STOCfree_data_check(data = data, data_name = data_nm)

  print(data$inf_dyn_priors)

  ## checking for missing and impossible parameter values
  check_priors_status_dyn(data = data, data_name = data_nm)

}


#' Sets the parameters associated with priors related to infection dynamics
#'
#' @param data a STOCfree_data object
#' @param pi1_a alpha parameter for the prior distribution of the probability of being infected on the first test
#' @param pi1_b beta parameter for the prior distribution of the probability of being infected on the first test
#' @param tau2_a alpha parameter for the prior distribution of the probability of remaining infected between consecutive months
#' @param tau2_b beta parameter for the prior distribution of the probability of remaining infected between consecutive months
#'
#' @return
#' @export
set_priors_status_dyn <- function(data = NULL, pi1_a = NULL, pi1_b = NULL, logit_pi1_mean = NULL, logit_pi1_sd = NULL,
                               tau1_a = NULL, tau1_b = NULL, logit_tau1_mean = NULL, logit_tau1_sd = NULL,
                               tau2_a = NULL, tau2_b = NULL, logit_tau2_mean = NULL, logit_tau2_sd = NULL,
                               pi_within_a = NULL, pi_within_b = NULL, logit_pi_within_mean = NULL, logit_pi_within_sd = NULL){

  ## checking that data is a STOCfree_data object
  data_nm <- deparse(substitute(data))
  STOCfree_data_check(data = data, data_name = data_nm)

    ## parameters with required values
  inf_dyn_names <- names(data$inf_dyn_priors)

  ## supplied arguments
  arg_names <- names(as.list(match.call()))
  arg_names <- arg_names[-match("data", arg_names)]

  ## arguments effectively used
  args_used <- inf_dyn_names[which(inf_dyn_names %in% arg_names)]

  par_val <- as.list(match.call())[args_used]
  par_val <- do.call("c", par_val)

  data$inf_dyn_priors[args_used] <- par_val

  ## checking for missing and impossible parameter values
  check_priors_status_dyn(data = data, data_name = data_nm)

  data

}

#' Plots the prior distributions for the parameters related to infection dynamics
#'
#' @param data a STOCfree_data object
#'
#' @return
#' @export
plot_priors_status_dyn <- function(data){

  ## checking that data is a STOCfree_data object
  data_nm <- deparse(substitute(data))
  STOCfree_data_check(data = data, data_name = data_nm)

  ## scale on which dynamics is modelled
  ## either "proba" or "logit"
  dyn_scale <- attr(data, "status dynamics scale")

  priors <- data$inf_dyn_priors
  n_risk_factors <- attr(data, "number of risk factors")
  plot_n_rows <- ifelse(n_risk_factors == 0, 2, 1)

  par(mfrow = c(plot_n_rows, 2))

  if(dyn_scale == "proba"){

  curve(dbeta(x, priors["pi1_a"], priors["pi1_b"]),
        from = 0, to = 1,
        main = "Probability of status positive\non the first test",
        xlab = "Probability",
        ylab = "Density")
  if(n_risk_factors == 0){

    curve(dbeta(x, priors["tau1_a"], priors["tau1_b"]),
          from = 0, to = 1,
          main = "Probability of becoming\nstatus positive between to 2 months",
          xlab = "tau1",
          ylab = "Density")

  }
  curve(dbeta(x, priors["tau2_a"], priors["tau2_b"]),
        from = 0, to = 1,
        main = "Probability of remaining\nstatus positive between to 2 months",
        xlab = "tau2",
        ylab = "Density")
  }

  if(dyn_scale == "logit"){

    curve(dnorm_logit(x, priors["logit_pi1_mean"], priors["logit_pi1_sd"]),
          from = 0, to = 1,
          main = "Probability of status positive\non the first test",
          xlab = "Probability",
          ylab = "Density")
    if(n_risk_factors == 0){

      curve(dnorm_logit(x, priors["logit_tau1_mean"], priors["logit_tau1_sd"]),
            from = 0, to = 1,
            main = "Probability of becoming\nstatus positive between to 2 months",
            xlab = "tau1",
            ylab = "Density")

    }
    curve(dnorm_logit(x, priors["logit_tau2_mean"], priors["logit_tau2_sd"]),
          from = 0, to = 1,
          main = "Probability of remaining\nstatus positive between to 2 months",
          xlab = "tau2",
          ylab = "Density")
  }


 }


## This function checks the missing prior paramaters
## as well as standard deviations that are equal or smaller than 0
check_priors_status_dyn <- function(data, data_name){

  ## Missing parameter values
  missing_param <- names(data$inf_dyn_priors[is.na(data$inf_dyn_priors)])

  if(length(missing_param) > 0){

    message(paste0("Set prior distributions for status dynamics using:\n", data_name, " <- set_priors_status_dyn(", data_name, ", ", paste(missing_param, collapse = " = , "), " = )"))

  }

  ## Standard deviations < 0
  sd_par <- grep("_sd", names(data$inf_dyn_priors))
  sd_par_neg <- data$inf_dyn_priors[sd_par]
  sd_par_neg <- names(sd_par_neg[!is.na(sd_par_neg) & sd_par_neg <=0])

  if(length(sd_par_neg) > 0){

    warning(paste0("Standard deviations for prior distributions should be > 0. \n", data_name, " <- set_priors_status_dyn(", data_name, ", ", paste(sd_par_neg, collapse = " = , "), " = )"))

  }

}
