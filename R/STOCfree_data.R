#' STOCfree_data class
#'
#' @param test_data a data.frame containing test results
#' @param test_herd_col name of the column with herd / farm identifiers
#' @param test_date_col name of the column with date of test
#' @param test_res_col name of the column with test results
#' @param test_level level at which the tests are performed. Must be either herd of animal
#' @param risk_factor_data a data.frame containing the risk factors
#' @param risk_herd_col name of the column with herd / farm identifiers
#' @param risk_date_col name of the column with date when the risk factors apply
#' @param risk_factor_col name of the column(s) with risk factor values
#' @param risk_factor_type risk factor type. Must be either continuous or categorical
#' @param test_name_col when several tests are used, name of the column containing the test names
#'
#' @return A list of class STOCfree_data.
#' @export
#'
#' @examples
STOCfree_data <- function(test_data = data.frame(),
                          test_herd_col = NULL,
                          test_date_col = NULL,
                          test_res_col = NULL,
                          test_name_col = NULL,
                          test_level = c("herd", "animal"),
                          risk_factor_data = NULL,
                          risk_herd_col = NULL,
                          risk_date_col = NULL,
                          risk_factor_col = NULL,
                          risk_factor_type = c("continuous", "categorical")){

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
  month_first <- format(min(as.Date(test_data[[test_date_col]])), "%Y-%m")
  ## Month of last test in the data
  month_last <- format(max(as.Date(test_data[[test_date_col]])), "%Y-%m")

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
  ## Adding month of test to the test_data
  test_data$date__1 <- format(as.Date(test_data[[test_date_col]]), "%Y-%m")

  ## List of all months in the dataset numbered
  all_months_list <- data.frame(
    date__1 = format(
      month_id = rep(NA),
      seq(from = as.Date(paste0(month_first, "-01")),
          to = as.Date(paste0(month_last, "-01")),
          by = "1 month"), "%Y-%m"),
    stringsAsFactors = FALSE  )

  all_months_list$month_id <- 1:nrow(all_months_list)

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


  ## Choice of the method to construct the STOCfree_data object
  if(test_level == "herd" & n_tests == 1 & n_risk_factors  > 0){

    sfd <- new_STOCfree_data.herd_1test_rf()

      }


  if(test_level == "herd" & n_tests > 1 & n_risk_factors  > 0){

    sfd <- new_STOCfree_data.herd_ntests_rf(test_data = test_data,
                                            test_herd_col = test_herd_col,
                                            test_date_col = test_date_col,
                                            test_res_col = test_res_col,
                                            test_name_col = test_name_col,
                                            test_level = test_level,
                                            risk_factor_data = risk_factor_data,
                                            risk_herd_col = risk_herd_col,
                                            risk_date_col = risk_date_col,
                                            risk_factor_col = risk_factor_col,
                                            risk_factor_type = risk_factor_type,
                                            all_tests = all_tests,
                                            month_last_id = month_last_id,
                                            herd_id_corresp = herd_id_corresp,
                                            all_months_list = all_months_list)

    }


  if(test_level == "animal" & n_risk_factors > 0){

    sfd <- new_STOCfree_data.animal_rf()

     }

  if(test_level == "herd" & n_tests == 1 & n_risk_factors  == 0){

    sfd <- new_STOCfree_data.herd_1test()

  }

  if(test_level == "animal" & n_risk_factors > 0){

    sfd <- new_STOCfree_data.animal_rf()

  }


  if(test_level == "animal") attr(sfd, "level")  <- "animal"

  attr(sfd, "number of herds")  <- n_herds
  attr(sfd, "number of tests")  <- n_tests
  attr(sfd, "month first test") <- month_first
  attr(sfd, "month last test")  <- month_last
  attr(sfd, "number of risk factors")  <- n_risk_factors

  sfd

  }

#############################################################################
##### Creation of STOCfree_data
##### herd level
##### 1 test
##### no risk factor
#############################################################################
new_STOCfree_data.herd_1test <- function(){

 stop("Function not defined when no risk factor provided")

}

#############################################################################
##### Creation of STOCfree_data
##### herd level
##### > 1 test
##### no risk factor
#############################################################################
new_STOCfree_data.herd_ntests <- function(){

 stop("Function not defined when no risk factor provided")


}

#############################################################################
##### Creation of STOCfree_data
##### animal level
##### 1 test
##### no risk factor
#############################################################################
new_STOCfree_data.animal <- function(){

 stop("Function not defined when no risk factor provided")


}

#############################################################################
##### Creation of STOCfree_data
##### herd level
##### 1 test
##### > 0 risk factor
#############################################################################
new_STOCfree_data.herd_1test_rf <- function(){


## Variable names - kept for later use by other functions
  var_names <- c(test_herd_col = test_herd_col,
                 test_date_col = test_date_col,
                 test_res_col = test_res_col,
                 risk_herd_col = risk_herd_col,
                 risk_date_col = risk_date_col)


}

#############################################################################
##### Creation of STOCfree_data
##### herd level
##### > 1 test
##### > 0 risk factor
#############################################################################
new_STOCfree_data.herd_ntests_rf <- function(test_data = test_data,
                                             test_herd_col = test_herd_col,
                                             test_date_col = test_date_col,
                                             test_res_col = test_res_col,
                                             test_name_col = test_name_col,
                                             test_level = test_level,
                                             risk_factor_data = risk_factor_data,
                                             risk_herd_col = risk_herd_col,
                                             risk_date_col = risk_date_col,
                                             risk_factor_col = risk_factor_col,
                                             risk_factor_type = risk_factor_type,
                                             all_tests = all_tests,
                                             month_last_id = month_last_id,
                                             herd_id_corresp = herd_id_corresp,
                                             all_months_list = all_months_list){

  ## Variable names - kept for later use by other functions
  var_names <- c(test_herd_col = test_herd_col,
                 test_date_col = test_date_col,
                 test_res_col = test_res_col,
                 test_name_col = test_name_col,
                 risk_herd_col = risk_herd_col,
                 risk_date_col = risk_date_col)

  ## Dataset with all test months
  test_data <- merge(test_data, all_tests,
                     by = c("herd_id", "month_id"),
                     all = TRUE)

  ## Adding row status_id to test_data
  ## these are the status ids that will be used by JAGS
  test_data <- test_data[order(test_data$herd_id, test_data$month_id),]
  test_data$status_id <- 1:nrow(test_data)

  ## dataset to be merged wit risk factor data
  rf_data_all_herds_months <- test_data[,
                                        c("status_id", "herd_id", "month_id")]

  ## going back to the data where test results are available
  test_data <- test_data[!is.na(test_data[, test_res_col]) | (is.na(test_data[, test_res_col]) & test_data$month_id == month_last_id),]

  ## herd test data
  herd_test_data <- by(test_data, list(test_data$herd_id), function(x){
    data.frame(
      herd_id = as.integer(unique(x$herd_id)),
      ind_i = as.integer(min(x$status_id)),
      ind_j = as.integer(min(x$status_id) + 1),
      ind_f = as.integer(max(x$status_id) - 1),
      ind_p = as.integer(max(x$status_id)),
      last_is_test = ifelse(is.na(x[x$month_id == month_last_id, test_res_col]), 0, 1)
    )
  })

  herd_test_data <- do.call("rbind", herd_test_data)

  ## test_data formatting
  test_data <- test_data[order(test_data$status_id),
                         c("status_id","herd_id", "month_id",
                           test_date_col, test_name_col, test_res_col)]
  rownames(test_data) <- 1:nrow(test_data)

  ## Priors for test performance
  test_names <- sort(unique(na.omit(test_data[[test_name_col]])))

  test_perf_prior <- data.frame(
    test = test_names,
    test_id = 1:length(test_names),
    Se_a = NA,
    Se_b = NA,
    Sp_a = NA,
    Sp_b = NA,
    stringsAsFactors = FALSE
  )

  ## test_id added to test_data
  ## used for indexing in the JAGS model
  test_data$test_id <- match(test_data[[test_name_col]], test_names)

  ## Priors for infection dynamics
  inf_dyn_priors <- c(
    pi1_a = NA,
    pi1_b = NA,
    tau2_a = NA,
    tau2_b = NA
  )

  if(is.null(risk_factor_data)){

    ## Subset of the risk factor data with herd id in test data
    sfd <-  list(
      var_names = var_names,
      herd_id_corresp = herd_id_corresp,
      test_data = test_data,
      herd_test_data = herd_test_data,
      test_perf_prior = test_perf_prior,
      risk_factors = character(),
      risk_factor_data = data.frame(),
      inf_dyn_priors = inf_dyn_priors
    )

    n_risk_factors <- 0

  } else {

    ## Risk factor list
    risk_factors <- data.frame(
      risk_factor = c("Intercept", risk_factor_col),
      type = c("intercept", risk_factor_type),
      modality = NA,
      ref = 0,
      stringsAsFactors = FALSE
    )

    risk_factors$type[grep("ca", tolower(risk_factors$type))] <- "categorical"

    ## identifcation of categorical risk factors
    rf_cat <- which(risk_factors$type == "categorical")

    rf_cat_list <- data.frame(
      risk_factor = NULL,
      type = NULL,
      modality = NULL,
      ref = NULL,
      stringsAsFactors = FALSE
    )

    if(length(rf_cat) > 0){

      for(i in length(rf_cat)){

        rf_i <- risk_factors$risk_factor[rf_cat[i]]
        rf_i_modalities <- na.omit(unique(unlist(risk_factor_data[, rf_i])))

        rf_cat_i <- data.frame(
          risk_factor = rep(rf_i, length(rf_i_modalities)),
          type = rep("categorical", length(rf_i_modalities)),
          modality = rf_i_modalities,
          ref = 0,
          stringsAsFactors = FALSE

        )

        rf_cat_list <- rbind(rf_cat_list, rf_cat_i)

      }

    }

    risk_factors <- risk_factors[risk_factors$type %in% c("intercept", "continuous"),]
    risk_factors <- rbind(risk_factors, rf_cat_list)

    risk_factors$n_obs <- rep(NA)

    for(i in 1:nrow(risk_factors)){

      if(risk_factors$type[i] %in% c("intercept", "continuous")){

        risk_factors$n_obs[i] <- nrow(risk_factor_data[!is.na(risk_factors$risk_factor[i]), ])

      } else {

        risk_factors$n_obs[i] <- nrow(risk_factor_data[risk_factor_data[, risk_factors$risk_factor[i]] == risk_factors$modality[i],])

      }

    }

    risk_factors$sd_prior <- risk_factors$mean_prior <- rep(NA)

    ## adding a reference modality to categorical risk factors
    rf_cat1 <- unique(risk_factors$risk_factor[risk_factors$type == "categorical"])

    for(j in 1:length(rf_cat1)){

      rf_j <- match(rf_cat1[j], risk_factors$risk_factor)[1]
      risk_factors$ref[rf_j] <- 1

    }

    ## risk factor data
    # selecting the columns that will be used
    risk_factor_data <- risk_factor_data[, c(risk_herd_col,
                                             risk_date_col,
                                             risk_factor_col)]
    ## Adding herd_id
    risk_factor_data <- merge(herd_id_corresp,
                              risk_factor_data,
                              by.x = test_herd_col,
                              by.y = risk_herd_col)
    ## Adding month_id
    risk_factor_data$date__1 <- format(as.Date(risk_factor_data[[risk_date_col]]), "%Y-%m")
    risk_factor_data <- merge(risk_factor_data, all_months_list, by = "date__1")

    ## merge with all tests
    risk_factor_data <- merge(rf_data_all_herds_months, risk_factor_data,
                              by = c("herd_id", "month_id"),
                              all.x = TRUE)

    risk_factor_data <- risk_factor_data[,
                                         c("status_id", "herd_id", "month_id", risk_factor_col)]

    ## Missing values are replaced by 0
    for(i in 1:length(risk_factor_col)){

      rf_i <- risk_factor_col[i]
      risk_factor_data[is.na(risk_factor_data[, rf_i]), rf_i] <- 0

    }

    for(i in 1:nrow(risk_factors)){

      rf_i <- risk_factors$risk_factor[i]
      type_i <- risk_factors$type[i]
      modalit_i <- risk_factors$modality[i]
      ref_i <- risk_factors$ref[i]

      if(i == 1) X <- data.frame(Intercept = rep(1, nrow(risk_factor_data)))
      if(type_i == "continuous"){


        X_i <- risk_factor_data[, rf_i]
        X <- cbind(X, X_i)
        colnames(X)[length(X)] <- rf_i

      }
      if(type_i == "categorical" & ref_i ==0){

        X_i <- ifelse(risk_factor_data[, rf_i] == modalit_i, 1, 0)
        X <- cbind(X, X_i)
        colnames(X)[length(X)] <- paste(rf_i, modalit_i, sep = "_")

      }

    }

    rf_data <- cbind(risk_factor_data[, 1:3], X)

    ## Subset of the risk factor data with herd id in test data
    sfd <-  list(
      var_names = var_names,
      herd_id_corresp = herd_id_corresp,
      test_data = test_data,
      herd_test_data = herd_test_data,
      test_perf_prior = test_perf_prior,
      risk_factors = risk_factors,
      risk_factor_data = rf_data,
      inf_dyn_priors = inf_dyn_priors
    )

    n_risk_factors <- length(which(sfd$risk_factors$ref == 0)) - 1

  }

  class(sfd) <- c("herd_ntests_rf", "STOCfree_data")

  attr(sfd, "level")  <- "herd"

  sfd

}


#############################################################################
##### Creation of STOCfree_data
##### animal level
##### 1 test
##### > 0 risk factor
#############################################################################
new_STOCfree_data.animal_rf <- function(){



}
