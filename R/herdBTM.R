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
#'   \item{TestResult}{test result. 0 for negative, 1 for positive}
#'   \item{ln_nOrig6_12}{risk factor 1, natural logarithm of the number of origins of the cattle purchased between 6 and 12 months before the date of test}
#'   \item{LocalSeroPrev}{risk factor 2, proportion of herds positive to the test in the municipality to which the herd belongs on the previous date of test}
#' }
"herdBTM"
