#' STOCfree_data class
#'
#' @param test_data a data.frame containing test results
#' @param test_herd_col name of the column with herd / farm identifiers
#' @param test_date_col name of the column with date of test
#' @param test_res_col name of the column with test results. Test results should be codes as 0 for negative results and 1 for positive results.
#' @param test_name_col when several tests are used, name of the column containing the test names
#' @param test_level level at which the tests are performed. Must be either herd of animal
#' @param test_N_anim name of the column with the number of animals to use as the denominator in the animal level model.
#' @param status_dynamics_scale scale on which priors for status dynamics are defined. The default is 'proba' and uses Beta priors. If 'logit' is used as an argument, the prior distributions will be normal distributions on the logit scale.
#' @param time_interval function used for risk factor aggregation. By default sum() is used, the values are added.
#' @param risk_factor_data a data.frame containing the risk factors
#' @param risk_herd_col name of the column with herd / farm identifiers
#' @param risk_date_col name of the column with date when the risk factors apply
#' @param risk_factor_col name of the column(s) with risk factor values
#' @param risk_factor_type risk factor type. Must be either continuous or categorical
#' @param lag1 start of the time interval for risk factor aggregation
#' @param lag2 end of the time interval for risk factor aggregation
#' @param FUN function used when aggregating the data. By default sum() is used, the values are added.
#'
#' @return A list of class STOCfree_data. A specific sub-class is defined based on the input.
#' @export
STOCfree_data <- function(test_data = data.frame(),
                          test_herd_col = NULL,
                          test_date_col = NULL,
                          test_res_col = NULL,
                          test_name_col = NULL,
                          test_level = c("herd", "animal"),
                          test_N_anim = NULL,
                          status_dynamics_scale = "proba",
                          risk_factor_data = NULL,
                          risk_herd_col = NULL,
                          risk_date_col = NULL,
                          risk_factor_col = NULL,
                          risk_factor_type = c("continuous", "categorical"),
                          lag1 = 0,
                          lag2 = 0,
                          time_interval = c("month", "year", "week"),
                          FUN = sum){


  time_interval <- match.arg(time_interval)
  time_in <- as.character(factor(time_interval,
                                 levels = c("month", "year", "week"),
                                 labels = c("%Y-%m", "%Y", "%Y-%W")))
  time_out <- as.character(factor(time_interval,
                                  levels = c("month", "year", "week"),
                                  labels = c("%Y-%m-%d", "%Y-%m-%d", "%Y-%W-%u")))
  ## Number of tests used
  if(is.null(test_res_col)) stop("Argument 'test_res_col' missing. Provide a column name for test results")

  if(is.null(test_name_col)){

    n_tests <- 1

  } else {

    n_tests <- length(na.omit(unique(test_data[[test_name_col]])))

  }

  ## Number of risk factors
  if(is.null(risk_factor_col)){

    n_risk_factors <- 0

  } else {

    n_risk_factors <- length(risk_factor_col)

  }

  ## Number of herds in the test dataset
  n_herds <- length(unique(test_data[[test_herd_col]]))
  ## Month of first test in the data
  month_first <- min(as.Date(test_data[[test_date_col]]))
  ## Month of last test in the data
  month_last <- max(as.Date(test_data[[test_date_col]]))
  ## List of all months in the dataset numbered
  all_months <- unique(format(seq(month_first, month_last, by = 1), time_in))
  ## Convert them back to a date and put them in a data.frame
  all_months_list <- data.frame(date__1 = as.Date(paste0(all_months, "-1-1"), format = time_out),
               month_id = seq_len(length(all_months)))
  month_first <- format(month_first, time_in)
  month_last <- format(month_last, time_in)

  ## List of all months in the dataset numbered
  date_min_herds <- all_months_list$date__1[which(all_months_list$month_id == (max(all_months_list$month_id) - 2))]
  ## First date in the dataset for each herd
  herd_first_date <- tapply(
    as.Date(test_data[[test_date_col]]),
    test_data[[test_herd_col]],
    function(x) as.character(min(x)))

  ## herds to discard because first date of recording is too close to the date to predict
  herds_discarded <- names(herd_first_date[herd_first_date >= date_min_herds])

  if(length(herds_discarded) > 0){

    test_data <- test_data[!test_data[[test_herd_col]] %in% herds_discarded, ]

    n_herds <- n_herds - length(herds_discarded)

    warning("The following herds were discarded because their first test result is too close to the month to predict:", paste(herds_discarded, collapse = ", "))

    }

  ## Error message if there are less than 2 herds remaining
  if(n_herds < 2) stop("Less than 2 herds. The model cannot be run.")

  ## old herd ids / new herd ids
  herd_id_corresp <- data.frame(
    old_herd_id = sort(unique(unlist(test_data[, test_herd_col]))),
    herd_id = 1:n_herds,
    stringsAsFactors = FALSE
  )

  colnames(herd_id_corresp)[1] <- test_herd_col

  ## renumbering herds
  test_data <- merge(herd_id_corresp, test_data,
                     by = test_herd_col)

  ## different tests used
  if(n_tests > 1){

    test_names <- sort(unique(na.omit(test_data[[test_name_col]])))
    test_data$test_id <- match(test_data[[test_name_col]], test_names)

  } else {

    test_names <- test_res_col
    test_data$test_id <- rep(1)

  }

  ## Priors for test performance
  test_perf_prior <- data.frame(
    test = test_names,
    test_id = 1:length(test_names),
    Se_a = NA,
    Se_b = NA,
    Sp_a = NA,
    Sp_b = NA,
    stringsAsFactors = FALSE
  )


  ## Adding month of test to the test_data
  test_data$date__1 <- as.Date(
    format(as.Date(test_data[[test_date_col]]), "%Y-%m-01"))

  ## ID of the last month in the dataset
  month_last_id <- max(all_months_list$month_id)

  ## Month ID added to the test_data dataset
  test_data <- merge(test_data, all_months_list)

  ## First test month of each herd
  herd_min_month <- tapply(test_data$month_id, test_data$herd_id,  min)

  ## Dataset with all possible months between first for each herd and last in dataset
  all_tests <- data.frame(
    herd_id = as.integer(names(herd_min_month)),
    month_min = herd_min_month
  )

  all_tests <- apply(all_tests, 1, function(x){
    data.frame(
      herd_id = rep(x["herd_id"],
                    month_last_id - as.integer(x["month_min"]) + 1),
      month_id = as.integer(x["month_min"]):month_last_id
    )
  })

  all_tests <- do.call("rbind", all_tests)

  ## Dataset with all test months
  test_data <- merge(test_data, all_tests,
                     by = c("herd_id", "month_id"),
                     all = TRUE)

  ## If level of testing is animal, test results is number pos / number tested
  if(test_level == "animal") test_data <- make_animal_test(test_data, test_res_col, test_N_anim)

  ## dataset to be merged with risk factor data
  ## numbers the status_id in the dataset
  rfd <- test_data[, c("herd_id", "month_id", "test_id")]
  rfd <- unique(rfd)
  rfd <- rfd[order(rfd$herd_id, rfd$month_id, rfd$test_id),]
  rfd$status_id <- 1:nrow(rfd)
  rfd <- rfd[, c("status_id", "herd_id", "month_id", "test_id")]

  ## Adding row status_id to test_data
  ## these are the status ids that will be used by JAGS
  test_data <- merge(test_data, rfd)

  ## dataset to be merged with risk factor data
  rfd$intercept <- rep(1)
  ## creation of risk factor list
  risk_factors <- data.frame(
    risk_factor = "Intercept",
    type = "intercept",
    modality = NA,
    stringsAsFactors = FALSE
  )

  ## going back to the data where test results are available
  if(test_level == "animal"){

    test_data <- test_data[!is.na(test_data[, "n_pos"]) | (is.na(test_data[, "n_pos"]) & test_data$month_id == month_last_id),]

  }

  if(test_level == "herd"){

    test_data <- test_data[!is.na(test_data[, test_res_col]) | (is.na(test_data[, test_res_col]) & test_data$month_id == month_last_id),]

  }

## ordering test_data by status_id and number of rows
  test_data <- test_data[order(test_data$status_id),]
  rownames(test_data) <- 1:nrow(test_data)
## status type - used in model
## default value
## other values added below
 test_data$status_type <- rep(2)

  ## Selecting the columns to keep from the test_data
  if(test_level == "animal"){

    cln <- c("status_id", "herd_id", "month_id", "status_type","test_id", "n_tested", "n_pos")

    test_data <- test_data[, cln]

    ## herd test data
    herd_test_data <- by(test_data, list(test_data$herd_id), function(x){
      data.frame(
        herd_id = as.integer(unique(x$herd_id)),
        ind_i = as.integer(min(x$status_id)),
        ind_p = as.integer(max(x$status_id))
      )
    })
    herd_test_data <- do.call("rbind", herd_test_data)


  }

  if(test_level == "herd"){

    test_data$test_res <- test_data[[test_res_col]]
    cln   <- c("status_id", "herd_id", "month_id", "status_type", "test_id", "test_res")

    test_data <- test_data[, cln]

    ## herd test data
    herd_test_data <- by(test_data, list(test_data$herd_id), function(x){
      data.frame(
        herd_id = as.integer(unique(x$herd_id)),
        ind_i = as.integer(min(x$status_id)),
        ind_p = as.integer(max(x$status_id))
      )
    })
    herd_test_data <- do.call("rbind", herd_test_data)

  }

  ## status type - used in model
  ## 1: first test in a herd
  ## 2: first test on a month which is not first test in herd
  ## 3: test > 1 on a month
  ## 4: status to predict without test result
  ## 5: status to predict with a single test performed
  ## 6: status to predict with several tests on this month
  dplctd_month <- which(duplicated(test_data[, c("herd_id", "month_id"),]))
  test_data$status_type[dplctd_month] <- rep(3)
  rm(dplctd_month)
  test_data$status_type[test_data$status_id %in% herd_test_data$ind_i] <- 1
  if(test_level == "herd")   test_data$status_type[test_data$status_id %in% herd_test_data$ind_p & is.na(test_data$test_res)]  <- 4
  if(test_level == "animal") test_data$status_type[test_data$status_id %in% herd_test_data$ind_p & is.na(test_data$n_pos)]  <- 4
  test_data$status_type[test_data$status_id %in% herd_test_data$ind_p & (is.na(test_data$test_res) | is.na(test_data$n_pos))]  <- 4
  test_data$status_type[test_data$status_id %in% herd_test_data$ind_p & test_data$status_type == 2] <- 5
  test_data$status_type[test_data$status_id %in% herd_test_data$ind_p & test_data$status_type == 3] <- 6


  ## Priors for infection dynamics
  ## all possibilities are listed and the unnecessary ones removed
  inf_dyn_priors <- make_dyn_prior_table(test_level = test_level,
                                   status_dynamics_scale = status_dynamics_scale)

  ## Variable names - for later use
  var_names <- c(test_herd_col = test_herd_col,
                 test_date_col = test_date_col,
                 test_res_col = test_res_col,
                 test_name_col = test_name_col,
                 risk_herd_col = risk_herd_col,
                 risk_date_col = risk_date_col)

  ## Creation of the STOCfree_data object
  sfd <-  list(
    var_names = var_names,
    herd_id_corresp = herd_id_corresp,
    test_data = test_data,
    herd_test_data = herd_test_data,
    test_perf_prior = test_perf_prior,
    risk_factors = risk_factors,
    risk_factor_data = rfd[, -match("test_id", colnames(rfd))],
    inf_dyn_priors = inf_dyn_priors
  )

  attr(sfd, "level")  <- test_level
  attr(sfd, "status dynamics scale")  <- status_dynamics_scale
  attr(sfd, "number of herds")  <- n_herds
  attr(sfd, "number of tests")  <- n_tests
  attr(sfd, "month first test") <- month_first
  attr(sfd, "month last test")  <- month_last
  attr(sfd, "number of risk factors")  <- n_risk_factors

  class(sfd) <- c(paste0(test_level,
                         ifelse(status_dynamics_scale == "logit", "_dynLogit", ""),
                         ifelse(n_risk_factors == 0, "", "_rf")),
                  "STOCfree_data")

  if(n_risk_factors > 0){

  sfd <- sf_add_risk_factor(sfd = sfd,
                     risk_factor_data = risk_factor_data,
                     risk_herd_col = risk_herd_col,
                     risk_date_col = risk_date_col,
                     risk_factor_col = risk_factor_col,
                     risk_factor_type = risk_factor_type,
                     lag1 = lag1,
                     lag2 = lag2,
                     FUN = FUN)

  }

  sfd

  }


## this function creates the table containing the parameters related to status dynamics
## it is called by STOCfree_data()
make_dyn_prior_table <- function(test_level, status_dynamics_scale){

  inf_dyn_priors <- NULL

  if(test_level == "herd" & status_dynamics_scale == "proba"){

    inf_dyn_priors <- c(
      pi1_a = NA,
      pi1_b = NA,
      tau1_a = NA,
      tau1_b = NA,
      tau2_a = NA,
      tau2_b = NA
    )
    }

  if(test_level == "herd" & status_dynamics_scale == "logit"){

    inf_dyn_priors <- c(
      logit_pi1_mean = NA,
      logit_pi1_sd = NA,
      logit_tau1_mean = NA,
      logit_tau1_sd = NA,
      logit_tau2_mean = NA,
      logit_tau2_sd = NA
    )
  }

  if(test_level == "animal" & status_dynamics_scale == "proba"){

    inf_dyn_priors <- c(
      pi1_a = NA,
      pi1_b = NA,
      pi_within_a = NA,
      pi_within_b = NA,
      tau1_a = NA,
      tau1_b = NA,
      tau2_a = NA,
      tau2_b = NA
    )
  }

  if(test_level == "animal" & status_dynamics_scale == "logit"){

    inf_dyn_priors <- c(
      logit_pi1_mean = NA,
      logit_pi1_sd = NA,
      logit_pi_within_mean = NA,
      logit_pi_within_sd = NA,
      logit_tau1_mean = NA,
      logit_tau1_sd = NA,
      logit_tau2_mean = NA,
      logit_tau2_sd = NA
    )
  }

  return(inf_dyn_priors)
}


#' Add a risk factor to STOCfree_data
#'
#' @param sfd STOC free data
#' @param risk_factor_data a data.frame with a risk factor
#' @param risk_herd_col name of the column with herd / farm identifiers
#' @param risk_date_col name of the column with date when the risk factors apply
#' @param risk_factor_col name of the column(s) with risk factor values
#' @param risk_factor_type risk factor type. Must be either continuous or categorical
#' @param lag1 start of the time interval for risk factor aggregation
#' @param lag2 end of the time interval for risk factor aggregation
#' @param FUN function used for risk factor aggregation. By default sum() is used, the values are added.
#'
#' @return
#' @export
#'
sf_add_risk_factor <- function(sfd,
                               risk_factor_data,
                               risk_herd_col = NULL,
                               risk_date_col = NULL,
                               risk_factor_col = NULL,
                               risk_factor_type = c("continuous", "categorical"),
                               risk_factor_ref = NULL,
                               lag1 = 0,
                               lag2 = 0,
                               FUN = sum){

  if(length(lag1) == 1 & length(risk_factor_col) > 1){ lag1 <- rep(lag1, length(risk_factor_col)) }
  if(length(lag2) == 1 & length(risk_factor_col) > 1){ lag2 <- rep(lag2, length(risk_factor_col)) }

  ## First month in the test data
  month_first <- attr(sfd, "month first test")

  ## Final risk factor data - to be populated
  rfd <- sfd$risk_factor_data

  ## Risk factor data to incorporate in the final data
  risk_factor_data <- risk_factor_data[, c(risk_herd_col, risk_date_col, risk_factor_col)]
  risk_factor_data <- merge(risk_factor_data,
                            sfd$herd_id_corresp)

  ## Adding month_id to risk factor data
  range_month <- range(as.Date(as.character(risk_factor_data[[risk_date_col]])))

  months_list <- make_month_id(
    start = range_month[1],
    end = range_month[2],
    orig = paste0(month_first, "-01")
  )

  risk_factor_data$date__1 <- as.Date(format(
    as.Date(risk_factor_data[[risk_date_col]]),
    "%Y-%m-01"))

  ## dates as months in the risk_factor_data made to correspond to the
  ## the dates in the STOCfree_data object
  risk_factor_data <- merge(risk_factor_data, months_list)

  ## the month_id column is renamed lagged_month
  ## this will be used to merge these data with the STOCfree data with a lag
  colnames(risk_factor_data)[match("month_id", colnames(risk_factor_data))] <- "lagged_month"

  ## Going through the list of all risk factors
  for(i in 1:length(risk_factor_col)){

    ## dataset for this risk factor
    risk_factor_data_i <- risk_factor_data[
      order(risk_factor_data$herd_id,
            risk_factor_data$lagged_month),
      c("herd_id", "lagged_month", risk_factor_col[i])]

    colnames(risk_factor_data_i)[3] <- "risk_factor"

    ## lagged months from lag1 to lag2
    rfd_i <- merge(rfd[, c("status_id", "herd_id", "month_id")],
                   data.frame(
                     lagged_month = lag1[i]:lag2[i]
                   ))

    ## each month_id in each herd is matched with the lagged months of interest
    rfd_i$lagged_month <- with(rfd_i, month_id - lagged_month)

    ## merging with risk factor data
    risk_factor_data_i <- merge(rfd_i,
                                risk_factor_data_i,
                                all.x = TRUE)

    ##########################
    ## Categorical risk factor
    ##########################
    if(risk_factor_type[i] == "categorical"){

      ## when several modalities observed for the same status (month)
      ## only the most frequent is kept
      risk_factor_data_i <- dplyr::group_by(risk_factor_data_i, status_id, herd_id, month_id)
      risk_factor_data_i <- dplyr::count(risk_factor_data_i, risk_factor)
      risk_factor_data_i <- dplyr::slice(risk_factor_data_i, which.max(n))

      ## list of modalities for this risk factor
      ls_modalities <- as.data.frame(
        table(risk_factor_data[[risk_factor_col[i]]]),
        stringsAsFactors = FALSE
      )

      ls_modalities <- ls_modalities[rev(order(ls_modalities$Freq)),]

      ## error message for categorical risk factors with 1 modality only
      if(nrow(ls_modalities) < 2) stop(paste("Cannot include", risk_factor_col[i], ": only 1 modality in the data"))

      ## removing the column with count of modality occurrences
      risk_factor_data_i <- risk_factor_data_i[, -match("n", colnames(risk_factor_data_i))]

      ## case where a reference modality is provided
      ## (when no reference provided, the most frequent modality is used as reference)
      if(!is.null(risk_factor_ref)){

        i_ref <- match(risk_factor_ref, ls_modalities[, 1])
        if(length(i_ref) == 0) stop(paste("The reference modality", risk_factor_ref, "does not exist"))
        mod_ref <- ls_modalities[i_ref, ]
        ls_modalities <- ls_modalities[-i_ref, ]
        ls_modalities <- ls_modalities[rev(order(ls_modalities$Freq)),]

        ls_modalities <- rbind(mod_ref, ls_modalities)

      }

      ls_modalities$modality <- 1:nrow(ls_modalities)

      colnames(ls_modalities) <- c("modality_name", "freq", "modality")
      ls_modalities <- ls_modalities[, c("modality", "modality_name", "freq")]

      ## updating risk factor list
      ### if risk factor already exists, it is removed
      if(paste(risk_factor_col[i], lag1[i], lag2[i], sep = "_") %in% sfd$risk_factors$risk_factor){

        i_row <- match(paste(risk_factor_col[i], lag1[i], lag2[i], sep = "_"), sfd$risk_factors$risk_factor)
        sfd$risk_factors <- sfd$risk_factors[-i_row,]

      }

      ## row with reference modality removed
      ls_modalities <- ls_modalities[-1, ]

      tab_rf_i <- data.frame(
        risk_factor = rep(paste(risk_factor_col[i], lag1[i], lag2[i], sep = "_")),
        type = rep(risk_factor_type[i]),
        modality = ls_modalities$modality_name,
        stringsAsFactors = FALSE)

      sfd$risk_factors <- rbind(sfd$risk_factors,
                                tab_rf_i)


      ## contrast matrix - each modality in a different column
      for(j in 1:nrow(ls_modalities)){

        risk_factor_data_i$tmp <- ifelse(risk_factor_data_i$risk_factor == ls_modalities$modality_name[j], 1, 0)
        ## when missing values, the reference category is used
        risk_factor_data_i$tmp[is.na(risk_factor_data_i$tmp)] <- 0

        ## column name for this modality
        ## if lag1 and lag2 = 0, lags not included in the name
        if(lag1[i] == 0 & lag2[i] == 0){

          nm_modality <- paste(risk_factor_col[i], ls_modalities$modality_name[j], sep = "_")

        } else {

          nm_modality <- paste(risk_factor_col[i], lag1[i], lag2[i], ls_modalities$modality_name[j], sep = "_")

        }

        colnames(risk_factor_data_i)[length(risk_factor_data_i)] <- nm_modality

        ## if the modality exists in the final risk factor dataset, it is removed
        if(nm_modality %in% colnames(rfd)){

          i_col <- match(nm_modality, colnames(rfd))
          rfd <- rfd[, -i_col]

        }

      }

      ## id of the column with risk factor
      j_col_rf <- match("risk_factor", colnames(risk_factor_data_i))
      ## risk factor column removed to keep contrasts only
      risk_factor_data_i <- risk_factor_data_i[, -j_col_rf]
      ## column names with the different modalities
      col_modalities <- colnames(risk_factor_data_i)[j_col_rf:length(risk_factor_data_i)]

    }
    ##########################
    ## continuous risk factors
    ##########################
    if(risk_factor_type[i] == "continuous"){

      ## Missing values values for continuous risk factors replaced by 0
      risk_factor_data_i$risk_factor[is.na(risk_factor_data_i$risk_factor)] <- 0

      ## aggregation by month id
      risk_factor_data_i <- dplyr::group_by(risk_factor_data_i, status_id, herd_id, month_id)
      risk_factor_data_i <- dplyr::summarise(risk_factor_data_i, risk_factor = FUN(risk_factor))



      if(lag1[i] == 0 & lag2[i] == 0){

        rf_name <- risk_factor_col[i]

      } else {

        rf_name <- paste(risk_factor_col[i], lag1[i], lag2[i], sep = "_")

      }

      ## renaming the column corresponding to the risk factor
      colnames(risk_factor_data_i)[match("risk_factor", colnames(risk_factor_data_i))] <- rf_name

      ## if the risk factor exists in the final risk factor dataset, it is removed
      if(rf_name %in% colnames(rfd)){

        i_col <- match(rf_name, colnames(rfd))
        rfd <- rfd[, -i_col]

      }

      ## updating risk factor list
      if(rf_name %in% sfd$risk_factors$risk_factor){

        i_row <- match(rf_name, sfd$risk_factors$risk_factor)
        sfd$risk_factors <- sfd$risk_factors[-i_row,]

      }
      sfd$risk_factors <- rbind(sfd$risk_factors,
                                data.frame(
                                  risk_factor = rf_name,
                                  type = risk_factor_type[i],
                                  modality = NA,
                                  stringsAsFactors = FALSE))


    }
  }

  ## Final risk factor data
  rfd <- merge(rfd,
               risk_factor_data_i,
               all.x = TRUE)


  rfd <- rfd[order(rfd$status),]
  rownames(rfd) <- 1:nrow(rfd)

  rfd

  sfd$risk_factor_data <- rfd

  ## If there are risk factors, there are no priors needed for tau1
  if(length(grep("tau1", names(sfd$inf_dyn_priors))) > 0){

    sfd$inf_dyn_priors <- sfd$inf_dyn_priors[-grep("tau1", names(sfd$inf_dyn_priors))]

  }
  #
  # ## Number of risk factors
  n_risk_factors <- nrow(sfd$risk_factors) - 1
  attr(sfd, "number of risk factors")  <- n_risk_factors

  test_level <- attr(sfd, "level")
  n_tests <- attr(sfd, "number of tests")

  class(sfd) <- c(paste0(test_level,
                         "_rf"),
                  "STOCfree_data")


  sfd

}


#' Remove risk factor from STOCfree_data
#'
#' @param sfd a STOCfree_data object
#' @param risk_factor risk factor name
#'
#' @return
#' @export
#'
sf_remove_risk_factor <- function(sfd,
                               risk_factor = character()){

  cln_to_rem <- match(risk_factor, colnames(sfd$risk_factor_data))
  sfd$risk_factor_data <- sfd$risk_factor_data[, -cln_to_rem]

  rf_to_rem <- match(risk_factor, sfd$risk_factors["risk_factor"])
  sfd$risk_factors <- sfd$risk_factors[-rf_to_rem, ]

 }


# Aggregate animal level data at the herd-month level
# Function called by STOCfree_data()
make_animal_test <- function(test_data, test_res_col, test_N_anim){

  ## if no denominator is provided for the binomial distribution, the the number of animals tested is used
  if(is.null(test_N_anim)){

    test_data <- by(test_data,
                    list(test_data$herd_id, test_data$month_id, test_data$test_id),
                    function(x){
                      data.frame(
                        herd_id = unique(x$herd_id),
                        month_id = unique(x$month_id),
                        test_id = unique(x$test_id),
                        n_tested = length(na.omit(x[[test_res_col]])),
                        n_pos = sum(x[[test_res_col]])
                      )
                    })

  } else {

    test_data <- by(test_data,
                    list(test_data$herd_id, test_data$month_id, test_data$test_id),
                    function(x){
                      data.frame(
                        herd_id = unique(x$herd_id),
                        month_id = unique(x$month_id),
                        test_id = unique(x$test_id),
                        n_tested = max(x[[test_N_anim]]),
                        n_pos = sum(x[[test_res_col]])
                      )
                    })

      }


  test_data <- do.call("rbind", test_data)
  test_data <- test_data[!is.na(test_data$n_tested),]
  # making sure the number of animals tested is at least as big as the number of positives
  test_data$n_tested <- ifelse(test_data$n_tested < test_data$n_pos, test_data$n_pos, test_data$n_tested)
  test_data <- test_data[order(test_data$herd_id, test_data$month_id),]
  rownames(test_data) <- 1:nrow(test_data)

  test_data

  }


## function that checks that an object is of the STOCfree_data class
STOCfree_data_check <- function(data, data_name){

  if(data_name == ""){

    stop(paste("Argument data is missing with no default."))

  }

  if(!"STOCfree_data" %in% class(data)){

     stop(paste0(data_name, " is not a STOCfree_data object."))
   }

 }
