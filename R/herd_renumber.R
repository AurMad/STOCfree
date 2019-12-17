#' Numbering of herds from 1 to number of herds
#'
#' @param data a dataset in the data.frame or tibble formats
#' @param herd_colname name of the column containing the herd identifier
#'
#' @examples
#' herd_renumber(data = herdBTM, herd_colname = "Farm")
#'
#' @export
#'
herd_renumber <- function(data, herd_colname){

  herd <- dplyr::enquo(herd_colname)

  ## Creation of a new herd id
  data <- dplyr::mutate(data, herd_id = match(!! herd, unique(!! herd)))
  ## re ordering columns
  data <- dplyr::select(data, herd_id, tidyselect::everything())

  data

}

