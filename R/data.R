#' Simple simulated data
#'
#' A simulated dataset to demonstrate how to model different types of measurement error and missing data using the 'inlami' package.
#'
#' @format ## `simple_data`
#' A data frame with 1000 rows and 4 columns:
#' \describe{
#'   \item{y}{Response variable}
#'   \item{x}{Covariate measured with error, both Berkson and classical error and missing observations}
#'   \item{x_true}{Correct version of the covariate with error}
#'   \item{z}{Error free covariate, correlated with x}
#' }
#' @source The dataset is simulated.
"simple_data"

#' Framingham heart study data
#'
#' A data set with observations of heart disease status systolic blood pressure (SBP) and smoking status.
#'
#' @format ## `framingham`
#' A data frame with 641 rows and 4 columns:
#' \describe{
#'   \item{disease}{A binary response, 1 if heart disease, 0 otherwise}
#'   \item{sbp1}{log(SBP − 50) at examination 1 (centered)}
#'   \item{sbp2}{log(SBP − 50) at examination 2 (centered)}
#'   \item{smoking}{Smoking status, 1 if smoking, 0 otherwise.}
#' }
#' @source Carrol book?
"framingham"
