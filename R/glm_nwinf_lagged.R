#' Generates an initial date from a reference date and a time lag in months
#'
#' @param date a date with format yyyy-mm-dd
#' @param time_lag an integer representing a number of months
#'
#' @return a date which is time_lag months before or after date depending on the sign of time_lag
#' @export
#'
#' @examples date_from_lag(date = "2020-01-01", time_lag = -1)
date_from_lag <- function(date = character(), time_lag = integer()){

  year  <- as.integer(format(as.Date(date), "%Y"))
  month <- as.integer(format(as.Date(date), "%m"))

  n_years <- time_lag %/% 12

  date_month <- 1 + time_lag - n_years * 12

  as.Date(
    as.character(
      paste0(year + n_years, "-", sprintf("%02d", date_month), "-01")))

}


#' Modelling of lagged relationship between risk factor occurrence and new infection
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
glm_nwinf_lagged <- function(sf_data,
                             rf_data,
                             rf_herd_col = character(),
                             rf_date_col = character(),
                             rf_col = character(),
                             lag1 = 0,
                             lag2 = 36){

  ## Month of first test -> month_id = 1
  month_first <- attr(sf_data, "month first test")
  ## Month of last test
  month_last <- attr(sf_data, "month last test")

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

  # nwinf_data <- STOCfree_data$test_data
  nwinf_data <- sf_data$test_data

  ## renaming column with test results
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

  ## Herd level data for new infection
  rf_first_last <- by(nwinf_data,
                      list(nwinf_data$herd_id),
                      function(x){
                        data.frame(
                          herd_id = unique(x$herd_id),
                          month_first = min(x$month_id) - lag2,
                          month_last = max(x$month_id) - lag1
                        )
                      })

  rf_first_last  <- do.call("rbind", rf_first_last)

  ## For each herd, list of months for which risk factors need to be considered
  rf_herd_months <- apply(rf_first_last, 1, function(x){
    data.frame(
      herd_id = rep(x["herd_id"], as.integer(x["month_last"]) - as.integer(x["month_first"]) + 1),
      month_id = as.integer(x["month_first"]):as.integer(x["month_last"])
    )})

  rf_herd_months <- do.call("rbind", rf_herd_months)

  ## Risk factor data - herd id added
  rf_data <- merge(sf_data$herd_id_corresp, rf_data, all.x = TRUE)
  rf_data$herd_id <- as.integer(rf_data$herd_id)

  ## month_id added
  rf_data$date__1 <- as.Date(paste0(format(as.Date(as.character(rf_data[[rf_date_col]])),
                                           "%Y-%m"), "-01"))
  rf_data <- merge(all_months, rf_data, by = "date__1")
  # rf_data <- rf_data[order(rf_data$herd_id, rf_data$month_id),]

  ## data for all herds for all necessary months
  rf_data <- merge(rf_herd_months, rf_data, all.x = TRUE)

  ## Missing values replaced with 0
  rf_data[[rf_col]][is.na(rf_data[[rf_col]])] <- 0

  ## sequence of lags to study
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

  nwinf_model_data <- merge(nwinf_model_data, rf_data, all = TRUE)

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

    l1l2$AIC[i] <- AIC(glm(nwinf ~ risk_factor,
                           data = model_data,
                           family = binomial(link = "logit")))

  }

  l1l2

}
