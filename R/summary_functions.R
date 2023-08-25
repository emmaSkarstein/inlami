#' Simplify the "raw" model summary for printing and plotting
#'
#' @param inlami_model the model returned from the fit_inlami function.
#'
#' @return A list of four data frames, containing the summaries for different components of the model. These are the coefficients of the model of interest, the coefficient of the variable with error, the coefficients of the imputation model, and the hyperparameters.
#' @export
simplify_inlami_model_summary <- function(inlami_model){
  # Extract the fixed effects and the hyperparameters from the inla summary
  fixed <- inlami_model$summary.fixed
  hyper <- inlami_model$summary.hyperpar

  # Identify the moi coefs and the imputation coefs, based on whether they start with alpha or beta.
  moi_coef <- dplyr::filter(fixed, grepl("beta.", rownames(fixed)))
  imp_coef <- dplyr::filter(fixed, grepl("alpha.", rownames(fixed)))

  error_coef <- dplyr::filter(hyper, grepl("Beta for ", rownames(hyper)))

  error_var_name <- sub(".*Beta for ", "", rownames(error_coef))
  rownames(error_coef) <- error_var_name

  other_model_hyperpar <- dplyr::filter(hyper, !grepl("Beta for ", rownames(hyper)))

  return(list(moi_coef = moi_coef,
              error_coef = error_coef,
              imp_coef = imp_coef,
              other_model_hyperpar = other_model_hyperpar))
}

#' Summary method for inlami
#'
#' Takes a fitted `inlami` object produced by
#' `fit_inlami` and produces a summary from it.
#'
#' @aliases summary.inlami print.summary.inlami
#' @param object model of class `inlami`.
#' @param ... other arguments
#'
#' @return `summary.inlami` returns an object of class `summary.inlami`, a list of components to print.
#' @method summary inlami
#' @rdname summary
#' @export
#'
#' @examples
#' simple_moi <- y ~ x + z
#' simple_imp <- x ~ z
#'
#' # Prior for beta.x
#' prior.beta <- c(0, 1/1000) # N(0, 10^3)
#'
#' # Priors for y, measurement error and true x-value precision
#' prior.prec.y <- c(0.5, 0.5) # Gamma(0.5, 0.5)
#' prior.prec.u_b <- c(10, 9) # Gamma(0.5, 0.5)
#' prior.prec.u_c <- c(10, 9) # Gamma(0.5, 0.5)
#' prior.prec.r <- c(0.5, 0.5) # Gamma(0.5, 0.5)
#'
#' # Initial values
#' initial.prec.y <- 1
#' initial.prec.u_b <- 1
#' initial.prec.u_c <- 1
#' initial.prec.r <- 1
#'
#' # Fit the model
#' simple_model <- fit_inlami(data = simple_data,
#'                          formula_moi = simple_moi,
#'                          formula_imp = simple_imp,
#'                          family_moi = "gaussian",
#'                          error_type = c("berkson", "classical"),
#'                          prior.prec.y = prior.prec.y,
#'                          prior.prec.u_b = prior.prec.u_b,
#'                          prior.prec.u_c = prior.prec.u_c,
#'                          prior.prec.r = prior.prec.r,
#'                          initial.prec.y = initial.prec.y,
#'                          initial.prec.u_b = initial.prec.u_b,
#'                          initial.prec.u_c = initial.prec.u_c,
#'                          initial.prec.r = initial.prec.r)
#'
#' summary(simple_model)
summary.inlami <- function(object, ...){
  simple_summary <- simplify_inlami_model_summary(object)

  class(simple_summary) <- "summary.inlami"

  return(simple_summary)
}

#' Print method for summary.inlami
#'
#' @param x object of class summary.inlami.
#' @param ... other arguments
#'
#' @method print summary.inlami
#' @rdname summary
#' @export
#'
print.summary.inlami <- function(x, ...){
  # inlami_summary <- summary(inlami_model)

  cat("Fixed effects for model of interest: \n")
  print(x$moi_coef)
  cat("\n")

  cat("Coefficient for error prone variable: \n")
  print(x$error_coef)
  cat("\n")

  cat("Fixed effects for imputation model: \n")
  print(x$imp_coef)
  cat("\n")

  cat(paste0("Model hyperparameters (apart from ",
             rownames(x$error_coef), "): \n"))
  print(x$other_model_hyperpar)
  cat("\n")
}


