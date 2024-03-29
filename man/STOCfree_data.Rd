% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sf_data.R
\name{STOCfree_data}
\alias{STOCfree_data}
\title{STOCfree_data class}
\usage{
STOCfree_data(
  test_data = data.frame(),
  test_herd_col = NULL,
  test_date_col = NULL,
  test_res_col = NULL,
  test_name_col = NULL,
  status_dynamics_scale = "logit",
  risk_factor_data = NULL,
  risk_herd_col = NULL,
  risk_date_col = NULL,
  risk_factor_col = NULL,
  risk_factor_type = c("continuous", "categorical"),
  lag1 = 0,
  lag2 = 0,
  time_interval = c("month", "year", "week"),
  FUN = sum
)
}
\arguments{
\item{test_data}{a data.frame containing test results}

\item{test_herd_col}{name of the column with herd / farm identifiers}

\item{test_date_col}{name of the column with date of test}

\item{test_res_col}{name of the column with test results. Test results should be codes as 0 for negative results and 1 for positive results.}

\item{test_name_col}{when several tests are used, name of the column containing the test names}

\item{status_dynamics_scale}{scale on which priors for status dynamics are defined. The default is 'proba' and uses Beta priors. If 'logit' is used as an argument, the prior distributions will be normal distributions on the logit scale.}

\item{risk_factor_data}{a data.frame containing the risk factors}

\item{risk_herd_col}{name of the column with herd / farm identifiers}

\item{risk_date_col}{name of the column with date when the risk factors apply}

\item{risk_factor_col}{name of the column(s) with risk factor values}

\item{risk_factor_type}{risk factor type. Must be either continuous or categorical}

\item{lag1}{start of the time interval for risk factor aggregation}

\item{lag2}{end of the time interval for risk factor aggregation}

\item{time_interval}{function used for risk factor aggregation. By default sum() is used, the values are added.}

\item{FUN}{function used when aggregating the data. By default sum() is used, the values are added.}
}
\value{
A list of class STOCfree_data. A specific sub-class is defined based on the input.
}
\description{
STOCfree_data class
}
