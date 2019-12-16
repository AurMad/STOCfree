#' Expansion of longitudinal datasets to obtain 1 row per month
#'
#' @param data a dataset in the data.frame or tibble formats with 1 row per test date
#' @param herd_colname name of the column containing the herd identifier
#' @param date_colname name of the column containing the test date
#' @param test_res_colname name of the column containing test results
#'
#' @examples
#' data(herdBTM)
#' expand_month(data = herdBTM, herd_colname = Farm, date_colname = DateOfTest, test_res_colname = TestResult)
#'
#' @importFrom magrittr %>%
#'
#' @export
expand_month <- function(data, herd_colname, date_colname, test_res_colname){

  herd <- dplyr::enquo(herd_colname)
  date <- dplyr::enquo(date_colname)
  test <- dplyr::enquo(test_res_colname)

  data <- herd_renumber(data = data, herd_colname = !! herd)

  data <- data %>%
    dplyr::mutate(date__1 = as.Date(paste0(format(as.Date(!! date), "%Y-%m"), "-01")))

  herd_dates <- data %>%
    dplyr::select(herd_id, date__1)

  ## First and last months in the dataset
  date_min <- as.Date(min(herd_dates$date__1))
  date_max <- as.Date(max(herd_dates$date__1))
  ## Dataset with all the months listed and numbered
  all_months_list <- dplyr::tibble(
    date__1 = seq(date_min, date_max, by = "1 month"))  %>%
    dplyr::mutate(
      month_id = 1:dplyr::n()
    )

  month_max <- max(all_months_list$month_id)

  ## Month ID added to the herdDates dataset
  herd_dates <- dplyr::left_join(herd_dates, all_months_list)

  ## First and last months of test for each herd
  herd_first_last <- herd_dates %>%
    dplyr::group_by(herd_id) %>%
    dplyr::summarise(
      month_min = min(month_id),
      month_max = max(month_id)
    )

  ## Sequence of all possible test month IDs for all herds
  nm <-  colnames(herd_first_last)[1]

  all_tests <- apply(herd_first_last, 1, function(x){
    data.frame(
      nameTmp = rep(x["herd_id"], month_max - x["month_min"] + 1),
      month_id = x["month_min"]:month_max
    )
  })

  all_tests <- do.call("rbind", all_tests)
  colnames(all_tests)[1] <- colnames(herd_first_last)[1]

  all_tests <- dplyr::as_tibble(all_tests)

  ## Month added to the sequence of all possible test month IDs for all herds
  all_tests <- dplyr::left_join(all_tests, all_months_list) %>%
    dplyr::select(1, 3, 2)

  all_tests <- dplyr::left_join(all_tests, data)

  all_tests <- all_tests %>%
    dplyr::arrange(herd_id, month_id) %>%
    dplyr::mutate(row_id = row_number()) %>%
    dplyr::select(row_id, tidyselect::everything())

  all_tests
  herd <- group_by_herd(test_data = all_tests,
                        herd_colname = herd_id,
                        month_colname = month_id,
                        row_id_colname = row_id,
                        test_res_colname = !! test)

  JAGS_data <- list(
    test_data = all_tests,
    herd_data = herd
  )

}
