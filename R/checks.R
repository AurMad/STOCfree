#' Extract data from a STOCfree_data object based on status ids
#'
#' @param sfd a STOCfree_data object
#' @param status_id status id usually returned by error message
#'
#' @details This function is designed to help debugging error messages when compiling the JAGS model. Error messages can refer to Status whih is the internal variable refering to a month in a herd. The function will extract herds associated with specific statuses.
#'
#' @return
#' @export
#'
extract_herd_from_status_id <- function(sfd,
                                        status_id = integer()){

  ## herd level test data
  herd_test <- sfd$herd_test_data

  ## list of herds corresponding to the status_id provided
  herds <- rep(NA, length(status_id))

  for(i in 1:length(status_id)){

    herds[i] <- herd_test$herd_id[herd_test$ind_i <= status_id[i] & herd_test$ind_p >= status_id[i]]

  }

  sfd$herd_id_corresp <- sfd$herd_id_corresp[sfd$herd_id_corresp$herd_id %in% herds,]
  sfd$test_data <- sfd$test_data[sfd$test_data$herd_id %in% herds,]
  sfd$herd_test_data <- sfd$herd_test_data[sfd$herd_test_data$herd_id %in% herds,]

  herd_subset <- list(
    herd_id_corresp = sfd$herd_id_corresp[sfd$herd_id_corresp$herd_id %in% herds,],
    test_data = sfd$test_data[sfd$test_data$herd_id %in% herds,],
    herd_test_data = sfd$herd_test_data[sfd$herd_test_data$herd_id %in% herds,]
  )

  herd_subset
}
