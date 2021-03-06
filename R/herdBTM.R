#' Bulk tank milk ELISA test results for 100 herds
#'
#' A dataset containing results of ELISA test results performed approximately
#' every 6 months between 2010 and 2017 in 100 herds
#'
#' @docType data
#'
#' @usage data(herdBTM)
#'
#' @format A data frame with 1100 rows and 6 variables:
#' \describe{
#'   \item{Farm}{farm ID}
#'   \item{DateOfTest}{date of test}
#'   \item{ODR}{result of test as optical density ratio}
#'   \item{Test}{type of test used. Either bulk tank milk ELISA or confirmatory test}
#'   \item{TestResult}{test result. 0 for negative, 1 for positive}
#'   \item{LocalSeroPrev}{risk factor 2, proportion of herds positive to the test in the municipality to which the herd belongs on the previous date of test}
#' }
"herdBTM"
