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
#' @return A list of class STOCfree_data. A specific sub-class is defined based on the input.
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
                          risk_factor_type = c("continuous", "categorical"),
                          lag1 = 0,
                          lag2 = 0,
                          FUN = sum){

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
  ## List of all months in the dataset numbered
  all_months_list <- make_month_id(start = paste0(month_first, "-01"),
                                   end = paste0(month_last, "-01"))

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
  if(test_level == "animal") test_data <- make_animal_test(test_data, test_res_col)

  ## dataset to be merged with risk factor data
  ## numbers the status_id in the dataset
  rfd <- test_data[, c("herd_id", "month_id")]
  rfd <- unique(rfd)
  rfd <- rfd[order(rfd$herd_id, rfd$month_id),]
  rfd$status_id <- 1:nrow(rfd)
  rfd <- rfd[, c("status_id", "herd_id", "month_id")]

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
    ref = 0,
    stringsAsFactors = FALSE
  )

  ## going back to the data where test results are available
  if(test_level == "animal"){

    test_data <- test_data[!is.na(test_data[, "n_pos"]) | (is.na(test_data[, "n_pos"]) & test_data$month_id == month_last_id),]

  }

  if(test_level == "herd"){

    test_data <- test_data[!is.na(test_data[, test_res_col]) | (is.na(test_data[, test_res_col]) & test_data$month_id == month_last_id),]

  }

  ## Selecting the columns to keep from the test_data
  ## Selecting the columns to keep from the test_data
  if(test_level == "animal"){

    cln <- c("status_id", "herd_id", "month_id", "test_id", "n_tested", "n_pos")

    test_data <- test_data[, cln]

    ## herd test data
    herd_test_data <- by(test_data, list(test_data$herd_id), function(x){
      data.frame(
        herd_id = as.integer(unique(x$herd_id)),
        ind_i = as.integer(min(x$status_id)),
        ind_j = as.integer(min(x$status_id) + 1),
        ind_f = as.integer(max(x$status_id) - 1),
        ind_p = as.integer(max(x$status_id)),
        last_is_test = ifelse(is.na(x[x$month_id == month_last_id, "n_pos"]), 0, 1)
      )
    })
    herd_test_data <- do.call("rbind", herd_test_data)


  }

  if(test_level == "herd"){

    test_data$test_res <- test_data[[test_res_col]]
    cln   <- c("status_id", "herd_id", "month_id", "test_id", "test_res")

    test_data <- test_data[, cln]

    ## herd test data
    herd_test_data <- by(test_data, list(test_data$herd_id), function(x){
      data.frame(
        herd_id = as.integer(unique(x$herd_id)),
        ind_i = as.integer(min(x$status_id)),
        ind_j = as.integer(min(x$status_id) + 1),
        ind_f = as.integer(max(x$status_id) - 1),
        ind_p = as.integer(max(x$status_id)),
        last_is_test = ifelse(is.na(x[x$month_id == month_last_id, "test_res"]), 0, 1)
      )
    })
    herd_test_data <- do.call("rbind", herd_test_data)

  }

  ## Priors for infection dynamics
  ## all possible possibilities are listed and the unecessary ones removed
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

  if(test_level == "herd"){

    inf_dyn_priors <- inf_dyn_priors[-grep("pi_within", names(inf_dyn_priors))]

  }

  if(n_risk_factors > 0){

    inf_dyn_priors <- inf_dyn_priors[-grep("tau1", names(inf_dyn_priors))]

  }

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
    risk_factor_data = rfd,
    inf_dyn_priors = inf_dyn_priors
  )

  attr(sfd, "level")  <- test_level
  attr(sfd, "number of herds")  <- n_herds
  attr(sfd, "number of tests")  <- n_tests
  attr(sfd, "month first test") <- month_first
  attr(sfd, "month last test")  <- month_last
  attr(sfd, "number of risk factors")  <- n_risk_factors

  class(sfd) <- c(paste0(test_level,
                         paste0("_", ifelse(n_tests == 1, "1test", "ntests")),
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




#' Add a risk factor to STOCfree_data
#'
#' @param sfd STOC free data
#' @param risk_factor_data a data.frame with a risk factor
#' @param risk_herd_col
#' @param risk_date_col
#' @param risk_factor_col
#' @param risk_factor_type
#' @param lag1
#' @param lag2
#' @param FUN function used between lag1 and lag2. The default is a sum.
#'
#' @return
#' @export
#'
#' @examples
sf_add_risk_factor <- function(sfd,
                                  risk_factor_data,
                                  risk_herd_col = NULL,
                                  risk_date_col = NULL,
                                  risk_factor_col = NULL,
                                  risk_factor_type = c("continuous", "categorical"),
                                  lag1 = 0,
                                  lag2 = 0,
                                  FUN = sum){

  if(length(lag1) == 1 & length(risk_factor_col) > 1){ lag1 <- rep(lag1, length(risk_factor_col)) }
  if(length(lag2) == 1 & length(risk_factor_col) > 1){ lag2 <- rep(lag2, length(risk_factor_col)) }

  ## First month in the test data
  month_first <- attr(sfd, "month first test")

  ## Final risk factor data
  rfd <- sfd$risk_factor_data

  ## Risk factor data to incorporate in the final data
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

  risk_factor_data <- merge(risk_factor_data, months_list)

  colnames(risk_factor_data)[match("month_id", colnames(risk_factor_data))] <- "lagged_month"

  for(i in 1:length(risk_factor_col)){

  ## final risk factor data
  risk_factor_data_i <- risk_factor_data[
    order(risk_factor_data$herd_id,
          risk_factor_data$lagged_month),
    c("herd_id", "lagged_month", risk_factor_col[i])]

  colnames(risk_factor_data_i)[3] <- "risk_factor"

  ## lagged variables
  rfd_i <- merge(rfd[, c("status_id", "herd_id", "month_id")],
                 data.frame(
                 lagged_month = lag1[i]:lag2[i]
               ))

  rfd_i$lagged_month <- with(rfd_i, month_id - lagged_month)

  ## merging with risk factor data
  risk_factor_data_i <- merge(rfd_i,
                            risk_factor_data_i,
                            all.x = TRUE)

  risk_factor_data_i$risk_factor[is.na(risk_factor_data_i$risk_factor)] <- 0

  ## aggregation by month id
  rf_aggr <- dplyr::group_by(risk_factor_data_i, status_id, herd_id, month_id)
  rf_aggr <- dplyr::summarise(rf_aggr, risk_factor = FUN(risk_factor))

  if(risk_factor_type[i] == "categorical"){

    rf_aggr$risk_factor <- with(rf_aggr,
                                ifelse(risk_factor >= 1, 1, 0))

  }

  ## Final risk factor data
  rfd <- merge(
    rfd, rf_aggr,
    all.x = TRUE)

  if(lag1[i] == 0 & lag2[i] == 0){

    rf_name <- risk_factor_col[i]

  } else {

    rf_name <- paste(risk_factor_col[i], lag1[i], lag2[i], sep = "_")

    }


  if(!is.na(match(rf_name, colnames(rfd)))){

    sfd <- sf_remove_risk_factor(sfd,
                          risk_factor = rf_name)

  }
  colnames(rfd)[match("risk_factor", colnames(rfd))] <- rf_name

  ## updating risk factor list
  sfd$risk_factors <- rbind(sfd$risk_factors,
                            data.frame(
                              risk_factor = rf_name,
                              type = risk_factor_type[i],
                              modality = ifelse(risk_factor_type[i] == "categorical", 1, NA),
                              ref = 0,
                              stringsAsFactors = FALSE))

  }

  rfd <- rfd[order(rfd$status),]
  rownames(rfd) <- 1:nrow(rfd)

  sfd$risk_factor_data <- rfd

  ## If there are risk factors, there are no priors needed for tau1
  if(length(grep("tau1", names(sfd$inf_dyn_priors))) > 0){

    sfd$inf_dyn_priors <- sfd$inf_dyn_priors[-grep("tau1", names(sfd$inf_dyn_priors))]

  }

  ## Number of risk factors
  n_risk_factors <- length(sfd$risk_factor_data) - 4
  attr(sfd, "number of risk factors")  <- n_risk_factors

  test_level <- attr(sfd, "level")
  n_tests <- attr(sfd, "number of tests")

  class(sfd) <- c(paste0(test_level,
                         paste0("_", ifelse(n_tests == 1, "1test", "ntests")),
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
#' @examples
sf_remove_risk_factor <- function(sfd,
                               risk_factor = character()){

  cln_to_rem <- match(risk_factor, colnames(sfd$risk_factor_data))
  sfd$risk_factor_data <- sfd$risk_factor_data[, -cln_to_rem]

  rf_to_rem <- match(risk_factor, sfd$risk_factors["risk_factor"])
  sfd$risk_factors <- sfd$risk_factors[-rf_to_rem, ]

 }

#' Aggregate animal level data at the herd-month level
#'
#' @param test_data
#'
#' @return
#' @export
#'
#' @examples
make_animal_test <- function(test_data, test_res_col){

  test_data <- by(test_data,
                list(test_data$herd_id, test_data$month_id),
                function(x){
                  data.frame(
                    herd_id = unique(x$herd_id),
                    month_id = unique(x$month_id),
                    test_id = unique(x$test_id),
                    n_tested = length(na.omit(x[[test_res_col]])),
                    n_pos = sum(x[[test_res_col]])
                  )
                })

  test_data <- do.call("rbind", test_data)
  test_data <- test_data[order(test_data$herd_id, test_data$month_id),]
  rownames(test_data) <- 1:nrow(test_data)

  test_data

  }




