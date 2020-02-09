#' Creates a data.frame with first date of month and corresponding month id
#'
#' @param start first day of the date sequence as a character string
#' @param end last day of the date sequence as a character string
#' @param orig origin of the sequence. The corresponding month is labelled 1 in the month sequence
#'
#' @return a data.frame
#' @export
#'
#' @examples
#' make_month_id(start = "2019-12-21", end = "2020-03-30", orig = "2020-02-20")
#'
make_month_id <- function(start = character(),
                          end = character(),
                          orig = start){

  date_start <- paste0(format(as.Date(start), "%Y-%m"), "-01")
  date_end   <- paste0(format(as.Date(end), "%Y-%m"), "-01")
  date_orig  <- paste0(format(as.Date(orig), "%Y-%m"), "-01")

  month_list <- data.frame(
    date__1 = seq(
      as.Date(date_start),
      as.Date(date_end),
      by = "1 month"),
    month_id = rep(NA)
  )

  id_month_orig <- which(month_list$date__1 == date_orig)

  month_list$month_id <- 1:nrow(month_list) - id_month_orig + 1

  month_list

}

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
