#' Simple simulated data
#'
#' A simulated dataset to demonstrate how to model different types of measurement error and missing data using the 'inlamemi' package.
#'
#' @format ## `simple_data`
#' A data frame with 1000 rows and 4 columns:
#' \describe{
#'   \item{y}{Response variable}
#'   \item{x}{Covariate measured with error, both Berkson and classical error and missing observations}
#'   \item{x_true}{Correct version of the covariate with error}
#'   \item{z}{Error free covariate, correlated with x.}
#' }
#' @source The dataset is simulated.
"simple_data"
