#' Generates an initial date from a reference date and a time lag in months
#'
#' @param date a date with format yyyy-mm-dd
#' @param time_lag an integer representing a number of months
#'
#' @return a date which is time_lag months before or after date depending on the sign of time_lag
#' @export
#'
#' @examples seq_dates <- seq(as.Date("2020-01-01"),
#'  as.Date("2020-12-01"),
#'  by = "1 month")
#'
#' data.frame(date = seq_dates,
#'  lagged_date = date_from_lag(
#'   date = sqt,
#'   time_lag = -2))
date_from_lag <- function(date = character(), time_lag = integer()){

  year  <- as.integer(format(as.Date(date), "%Y"))
  month <- as.integer(format(as.Date(date), "%m"))

  n_years <- (month + time_lag - 1) %/% 12

  date_month <- month + time_lag - n_years * 12

  as.Date(
    as.character(
      paste0(year + n_years, "-", sprintf("%02d", date_month), "-01")))

}

#' Creates a dataset with new infection events from herd level test results
#'
#' @param sfd a STOCfree_data object
#' @param time_of_inf a character string. When 'mid', the time of new infection is the median of the months between which the consecutive tests occurred.
#'
#' @return
#' @export
#'
#' @examples
make_nwinf_data <- function(sfd,
                            time_of_inf = c("mid", "first", "last")){


  month_first <- attr(sfd, "month first test")
  month_last  <- attr(sfd, "month last test")

  test_res_col <- sfd$var_names[["test_res_col"]]

  nwinf_data <- sfd$test_data
  nwinf_data <- nwinf_data[order(nwinf_data$status_id),]

  ## test_res_col name changed
  colnames(nwinf_data)[match(test_res_col, colnames(nwinf_data))] <- "test_res"

  nwinf_data$prev_herd_id    <- c(NA, nwinf_data$herd_id[-nrow(nwinf_data)])
  ## Previous month id
  nwinf_data$prev_month_id   <- c(NA, nwinf_data$month_id[-nrow(nwinf_data)])
  nwinf_data$prev_month_id   <- with(nwinf_data, ifelse(herd_id == prev_herd_id,
                                                        prev_month_id, NA))
  ## Previous test result
  nwinf_data$prev_test_res <- c(NA, nwinf_data$test_res[-nrow(nwinf_data)])
  nwinf_data$prev_test_res <- with(nwinf_data, ifelse(herd_id == prev_herd_id,
                                                      prev_test_res, NA))

  ## Only rows eligible for new infections are kept
  nwinf_data <- nwinf_data[!is.na(nwinf_data$test_res) &  !is.na(nwinf_data$prev_test_res) & nwinf_data$prev_test_res == 0,]

  ## New infection variable
  nwinf_data$nwinf <- with(nwinf_data, ifelse(test_res == 1, 1, 0))

  ## month_id is changed to the midpoint between previous and current test
  nwinf_data$mid_month_id <- with(nwinf_data,
                                  floor(prev_month_id + .5 *(month_id - prev_month_id)))

  ## Name of the column to keep for month_id
  ## by default, month_id is the month of the last test
  if(time_of_inf == "mid"){

    nwinf_data$month_id <- with(nwinf_data,
                                floor(prev_month_id + .5 *(month_id - prev_month_id)))

  }

  ## Name of the column to keep for month_id
  if(time_of_inf == "first"){

    nwinf_data$month_id <- nwinf_data$prev_month_id

  }

  nwinf_data <- nwinf_data[, c("status_id", "herd_id", "month_id",
                               "nwinf")]

  nwinf_data <- nwinf_data[order(nwinf_data$status_id),]
  rownames(nwinf_data) <- 1:nrow(nwinf_data)

  nwinf_list <- list(
    herd_id_corresp = sfd$herd_id_corresp,
    nwinf_data = nwinf_data)

  attr(nwinf_list, "month first test") <- month_first
  attr(nwinf_list, "month last test")  <- month_last
  attr(nwinf_list, "time of inf")      <- time_of_inf

  class(nwinf_list) <- "nwinf_data"

  nwinf_list

}


#' Adds a risk factor to a new infection dataset
#'
#' @param nwinf a dataset created with the make_nwinf_data() function
#' @param rf_data risk factor data. Ids for farms must be the same as in the original test results data.
#' @param rf_col
#' @param rf_date_col
#' @param lag1
#' @param lag2
#'
#' @return
#' @export
#'
#' @examples
add_risk_factor <- function(nwinf = nwinf_data(),
                            rf_data,
                            rf_col = character(),
                            rf_date_col = character(),
                            lag1 = 6,
                            lag2 = 12){

  nwinf_data <- nwinf$nwinf_data

  ## Month of first test -> month_id = 1
  month_first <- attr(nwinf, "month first test")
  ## Month of last test
  month_last <- attr(nwinf, "month last test")

  ## risk factor column in the risk factor data renamed
  colnames(rf_data)[match(rf_col, colnames(rf_data))] <- "risk_factor"

  ## List of months used in the study
  ## month_id = 0 for the first month in the STOCfree dataset
  rf_first_month <- date_from_lag(
    date = paste0(month_first, "-01"), time_lag = -lag2)

  all_months <- data.frame(
    date__1 = seq(as.Date(rf_first_month),
                  as.Date(paste0(month_last, "-01")), by = "1 month"),
    rf_month_id = rep(NA),
    stringsAsFactors = FALSE
  )

  all_months$rf_month_id <- 1:nrow(all_months) -
    which(all_months$date__1 == paste0(month_first, "-01"))

  ## name of column with old herd_id in risk factor data
  risk_herd_col <- colnames(nwinf$herd_id_corresp)[1]

  ## herd_id added to risk factor data
  rf_data <- merge(rf_data, nwinf$herd_id_corresp)

  ## month of risk factor occurrence
  rf_data$date__1 <- as.Date(paste0(
    format(rf_data[[rf_date_col]], "%Y-%m"), "-01"))

  rf_data <- merge(rf_data, all_months)

  ## columns of interest selected
  rf_data <- rf_data[, c("herd_id", rf_date_col, "rf_month_id", "risk_factor")]

  ## sequence of lags to study
  sq_lag <- expand.grid(month_id = sort(unique(nwinf_data$month_id)),
                        rf_month_id = lag1:lag2)

  nwinf_data <- merge(
    nwinf_data[, c("herd_id", "month_id", "nwinf")],
    sq_lag)

  nwinf_data$rf_month_id <- with(nwinf_data,
                                 month_id - rf_month_id)

  nwinf_data <- merge(nwinf_data, rf_data, all.x = TRUE)

  nwinf_data$risk_factor[is.na(nwinf_data$risk_factor)] <- 0

  ## data aggregation
  aggr_data <- dplyr::group_by(nwinf_data,
                               herd_id, month_id)
  aggr_data <- dplyr::summarise(aggr_data,
                                nwinf = unique(nwinf),
                                risk_factor = sum(risk_factor))

  ## Risk factor column renamed using risk factor name and lags
  colnames(aggr_data)[match("risk_factor", colnames(aggr_data))] <- paste(rf_col, lag1, lag2, sep = "_")

  nwinf$nwinf_data <- as.data.frame(aggr_data)

  nwinf

}


#' Modelling of all possible lagged relationships between risk factor occurrence and new infection
#'
#' @param sf_data a object of class STOCfree_data with test results
#' @param rf_data risk factor data
#' @param rf_herd_col name of the column containing the herd id. The ids should be the same as in the test data
#' @param rf_date_col name of the column the date of risk factor occurrence. The dates should be formatted as yyyy-mm-dd
#' @param rf_col name of the column with risk factor
#' @param lag1 an integer value representing the minimum number of months from which to consider for risk factor occurrence
#' @param lag2 an integer value representing the maximum number of months to consider for risk factor occurrence
#'
#' @details this function models a probability of new infection as a function of risk factor occurrence considered at various time intervals before a test is performed. New infection is defined, at the herd test level, as a positive test result recorded after a negative test result on the previous test. All herds with 2 test results of which the first one is negative are eligible for a new infection. The association with risk factors recorded between lag1 months and lag2 months before the second test result are modelled with logistic regression. For each combination of lag1 and lag2 the AIC of the logistic model is recorded. The best interval is the one with the lowest AIC.
#'
#'
#' @return a data.frame with a model AIC for each combination of lag1 (start of interval) and lag2 (end of interval)
#' @export
#'
#' @examples
logit_nwinf_lagged <- function(sf_data,
                             rf_data,
                             rf_date_col = character(),
                             rf_col = character(),
                             lag1 = 0,
                             lag2 = 36,
                             time_of_inf = c("mid", "first", "last")){

  ## Month of first test -> month_id = 1
  month_first <- attr(sf_data, "month first test")
  ## Month of last test
  month_last <- attr(sf_data, "month last test")
  # name of column with herd_id
  rf_herd_col <- colnames(sf_data$herd_id_corresp)[1]

  ## List of months used in the study
  ## month_id = 0 for the first month in the STOCfree dataset
  rf_first_month <- date_from_lag(
    date = paste0(month_first, "-01"), time_lag = -lag2)

  all_months <- data.frame(
    date__1 = seq(as.Date(rf_first_month),
                  as.Date(paste0(month_last, "-01")), by = "1 month"),
    month_id = rep(NA),
    stringsAsFactors = FALSE
  )

  all_months$month_id <- 1:nrow(all_months) -
    which(all_months$date__1 == paste0(month_first, "-01"))

  ## New infection data created from test results
  nwinf_data <- sf_data$test_data

  ## renaming columns with test results and test dates
  colnames(nwinf_data)[match(sf_data$var_names[c("test_res_col", "test_date_col")],
                             colnames(nwinf_data))] <- c("test_res", "test_date")

  nwinf_data <- nwinf_data[order(nwinf_data$status_id),]

  nwinf_data$prev_herd_id    <- c(NA, nwinf_data$herd_id[-nrow(nwinf_data)])
  ## Previous month id
  nwinf_data$prev_month_id   <- c(NA, nwinf_data$month_id[-nrow(nwinf_data)])
  nwinf_data$prev_month_id   <- with(nwinf_data, ifelse(herd_id == prev_herd_id,
                                                        prev_month_id, NA))
  ## Previous test result
  nwinf_data$prev_test_res <- c(NA, nwinf_data$test_res[-nrow(nwinf_data)])
  nwinf_data$prev_test_res <- with(nwinf_data, ifelse(herd_id == prev_herd_id,
                                                      prev_test_res, NA))

  ## Only rows eligible for new infections are kept
  nwinf_data <- nwinf_data[!is.na(nwinf_data$test_res) &  !is.na(nwinf_data$prev_test_res) & nwinf_data$prev_test_res == 0,]

  ## New infection variable
  nwinf_data$nwinf <- with(nwinf_data, ifelse(test_res == 1, 1, 0))

  ## month_id is changed to the midpoint between previous and current test
  nwinf_data$mid_month_id <- with(nwinf_data,
                                  floor(prev_month_id + .5 *(month_id - prev_month_id)))


  cln <- match(c("prev_month_id", "month_id", "mid_month_id"),
               colnames(nwinf_data))

  colnames(nwinf_data)[cln] <- c("first", "last", "mid")

  colnames(nwinf_data)[match(time_of_inf, colnames(nwinf_data))] <- "month_id"

  nwinf_data <- nwinf_data[, c("herd_id", "month_id", "nwinf")]

  ## Risk factor data - herd id added
  rf_data <- merge(sf_data$herd_id_corresp, rf_data, all.x = TRUE)
  rf_data$herd_id <- as.integer(rf_data$herd_id)

  ## month_id added
  rf_data$date__1 <- as.Date(paste0(format(as.Date(as.character(rf_data[[rf_date_col]])),
                                           "%Y-%m"), "-01"))
  rf_data <- merge(all_months, rf_data, by = "date__1")


  # ## sequence of lags to study
  sq_lag <- expand.grid(month_id = sort(unique(nwinf_data$month_id)),
                        time_lag = lag1:lag2)

  nwinf_model_data <- merge(
    nwinf_data[, c("herd_id", "month_id", "nwinf")],
    sq_lag)

  nwinf_model_data$rf_month_id <- with(nwinf_model_data,
                                       month_id - time_lag)

  rf_data <- rf_data[, c("herd_id", "month_id", rf_col)]
  colnames(rf_data)[2] <- "rf_month_id"
  colnames(rf_data)[3] <- "risk_factor"

  nwinf_model_data <- merge(nwinf_model_data, rf_data, all.x = TRUE)
  ## removing duplicates generated by the merge function
  nwinf_model_data <- nwinf_model_data[-which(duplicated(nwinf_model_data)),]
  ## Filling missing values with 0s
  nwinf_model_data$risk_factor[is.na(nwinf_model_data$risk_factor)] <- 0

  ## Modelling of lagged relationships

  # Grid of all possible lags to explore
  l1l2 <- expand.grid(
    lag1 = lag1:lag2,
    lag2 = lag1:lag2,
    l = rep(NA),
    AIC = rep(NA)
  )

  l1l2 <- l1l2[l1l2$lag2 >= l1l2$lag1,]
  l1l2$l <- with(l1l2, lag2 - lag1)

  ## length and sum of response variable - used for checking AIC values
  unique_y <- nwinf_model_data$nwinf[-which(duplicated(nwinf_model_data[, c("herd_id", "month_id")]))]
  lg_y  <- length(unique_y)
  sum_y <- sum(unique_y)

  for(i in 1:nrow(l1l2)){

    l1 <- l1l2$lag1[i]
    l2 <- l1l2$lag2[i]

    model_data <- dplyr::filter(nwinf_model_data,
                                time_lag >= l1 & time_lag <= l2)
    model_data <- dplyr::group_by(model_data,
                                  herd_id, month_id)
    model_data <- dplyr::summarise(model_data,
                                   nwinf = unique(nwinf),
                                   risk_factor = sum(risk_factor))

    lg_yi  <- length(model_data$nwinf)
    sum_yi <- sum(model_data$nwinf)

    if(lg_yi == lg_y & sum_yi == sum_y){

      l1l2$AIC[i] <- AIC(glm(nwinf ~ risk_factor,
                             data = model_data,
                             family = binomial(link = "logit")))
    } else {

      l1l2$AIC[i] <- NA

    }

  }

  rownames(l1l2) <- 1:nrow(l1l2)
  l1l2

}

#' Logistic regression for the probability of new infection
#'
#' @param nwinf a dataset created with the make_nwinf_data() and add_risk_factor() functions
#' @param risk_factors a character string with rhe risk factors to include in the regression
#'
#' @return results of a logistic regression model performed with the glm function
#' @export
#'
#' @examples
logit_nwinf <- function(nwinf = nwinf_data(),
                      risk_factors = character){

  data <- nwinf$nwinf_data

  formula <- as.formula(paste("nwinf ~ ",
                       paste(risk_factors, collapse = " + ")))

  log_reg <- glm(formula = formula,
                 data = data,
                 family = binomial(link = "logit"))

  log_reg

  }
