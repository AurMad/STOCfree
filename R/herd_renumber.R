#' Numbering of herds from 1 to number of herds
#'
#' @param data a dataset in the data.frame or tibble formats
#' @param herd_colname name of the column containing the herd identifier
#'
#' @return
#'
#' @examples
#' herd_renumber(data = herdBTM, herd_colname = Farm)
#'
#' @importFrom magrittr %>%
#'
#' @export
herd_renumber <- function(data, herd_colname){

  herd <- dplyr::enquo(herd_colname)

  data <- data %>%
    dplyr::mutate(herd_id = match(!! herd, unique(!! herd)))%>%
    dplyr::select(herd_id, tidyselect::everything())

  data <- dplyr::as_tibble(data)

  data

}
