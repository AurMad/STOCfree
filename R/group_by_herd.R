#' Generates a dataset with herd level indices to be used by JAGS
#'
#' @param data a dataset with each row representing a test result in a herd on a given month
#' @param herd_colname name of the column with herd identifier
#' @param row_id_colname name of the column with row identifier
#' @param month_colname name of the column with month of test
#' @param test_res_colname name of the column with test results
#'
#' @return
#' @export
#'
#'
group_by_herd <- function(data, herd_colname, month_colname, row_id_colname, test_res_colname){

  herd     <- dplyr::enquo(herd_colname)
  month    <- dplyr::enquo(month_colname)
  row_id   <- dplyr::enquo(row_id_colname)
  test_res <- dplyr::enquo(test_res_colname)

  data <- dplyr::rename(data,
                  herd_id = !!herd,
                  month = !!month,
                  row_id = !!row_id,
                  test_res = !!test_res)

    herd_indices <- dplyr::group_by(data, herd_id)
    herd_indices <-   dplyr::summarise( herd_indices,
       ind_i = row_id[month ==  min(month)],
       ind_j = row_id[month ==  min(month)] + 1,
       ind_f = row_id[month ==  max(month)] - 1,
       ind_p = row_id[month ==  max(month)],
       last_is_test = ifelse(is.na(test_res[month ==  max(month)]), 0, 1)
    )

    herd_indices
}
