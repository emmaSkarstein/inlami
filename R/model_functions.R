#' Extract and group variables from formulas
#'
#' Helper function that takes in the formulas for the model of interest and the
#' imputation model, and groups them into responses, covariates, covariate with
#' error and covariate(s) without error, for both sub-models.
#'
#' @param formula_moi an object of class "formula", describing the main model to be fitted.
#' @param formula_imp an object of class "formula", describing the imputation model for the mismeasured and/or missing observations.
#'
#' @return A list containing the names of the different variables of the model. The names of the elements in the list are "response_moi" (the response for the moi), "covariates_moi" (all covariates in the moi), "error_variable" (the name of the variable with error or missing data), "covariates_error_free" (the moi covariates without error), "response_imp" (imputation model response), "covariates_imp" (imputation model covariates).
#' @export
#'
#' @examples
#' f_moi <- y ~ x + z
#' f_imp <- x ~ z
#' extract_variables_from_formula(formula_moi = f_moi, formula_imp = f_imp)
extract_variables_from_formula <- function(formula_moi,
                                           formula_imp){
  if(!(methods::is(formula_moi, "formula") & methods::is(formula_imp, "formula"))){
    stop("One of the input objects is not of class 'formula'")
  }
  # Model of interest variables ------------------------------------------------
  # Extract name of response:
  response_moi <- all.vars(formula_moi)[1]

  # Covariates in moi:
  covariates_moi <- all.vars(formula_moi)[-1]

  # Identify error prone variable:
  error_variable <- all.vars(formula_imp)[1]

  # Extract index of error variable in formula:
  error_var_index <- which(covariates_moi == error_variable)

  # Error-free covariates:
  covariates_error_free <- covariates_moi[-error_var_index]

  # Imputation model variables -------------------------------------------------
  response_imp <- all.vars(formula_imp)[1]
  covariates_imp <- all.vars(formula_imp)[-1]

  return(list(response_moi = response_moi,
              covariates_moi = covariates_moi,
              error_variable = error_variable,
              covariates_error_free = covariates_error_free,
              response_imp = response_imp,
              covariates_imp = covariates_imp))
}


#' Make formula for measurement error and missing data model
#'
#' @param formula_moi an object of class "formula", describing the main model to be fitted.
#' @param formula_imp an object of class "formula", describing the imputation model for the mismeasured and/or missing observations.
#' @param error_type Type of error (one of "classical", "berkson", "missing")
#'
#' @return An object of class "formula".
#' @export
#'
#' @examples
#' f_moi <- y ~ x + z
#' f_imp <- x ~ z
#' make_inlami_formula(formula_moi = f_moi, formula_imp = f_imp, error_type = "classical")
make_inlami_formula <- function(formula_moi,
                              formula_imp,
                              error_type = "classical"){
  # Extract variables from formulas (error variable can be extracted from imputation model response)
  # How to differentiate between moi-vars and imp-vars? Bind with beta. and alpha.?
  # Add the copy-stuff (this will vary depending on error type)

  # Weaknesses:
  #

  # Extract and group all variables from formulas:
  vars <- extract_variables_from_formula(formula_moi = formula_moi,
                                         formula_imp = formula_imp)

  # Define covariates in output formula as all variables in moi formula
  #  except the response and error prone covariate:
  covariates_error_free_string <- paste0("beta.", vars$covariates_error_free)

  # Covariates for imputation model:
  covariates_imp_string <- paste0("alpha.", vars$covariates_imp)

  # The copy term to ensure the mismeasured variable is copied correctly through the models
  copy_term1 <- paste0("f(", paste0("beta.", vars$response_imp),
                       ", copy = 'id.x', hyper = list(beta = list(param = c(0, 1/1000), fixed = FALSE)))")
  # TODO: What does the "param = c(0, 1/1000)" above actually control? Should this be given as an argument?
  # I think it's the prior for beta.x?

  copy_term2 <- "f(id.x, weight.x, model='iid', values = 1:n, hyper = list(prec = list(initial = -15, fixed=TRUE)))"

  r_term <- "f(id.r, weight.r, model='iid', values = 1:n, hyper = list(prec = list(initial = -15, fixed = TRUE)))"

  formula <- paste("Y ~ -1",
                   "beta.0",
                   paste(covariates_error_free_string, collapse = " + "),
                   "alpha.0",
                   paste(covariates_imp_string, collapse = " + "),
                   copy_term1,
                   copy_term2,
                   r_term,
                   sep = " + ")

  return(stats::as.formula(formula))
}


#' Make matrices for joint model specification in INLA
#'
#' @param data A data frame with response and covariate variables for the main model and the imputation model.
#' @param formula_moi an object of class "formula", describing the main model to be fitted.
#' @param formula_imp an object of class "formula", describing the imputation model for the mismeasured and/or missing observations.
#' @param error_type Type of error (one or more of "classical", "berkson", "missing")
#'
#' @return A data frame with all the data needed for the model, structured in the right way for INLA.
#' @export
#'
#' @examples
#' f_moi <- y ~ x + z
#' f_imp <- x ~ z
#' make_inlami_matrices(data = simple_data,
#'               formula_moi = f_moi,
#'               formula_imp = f_imp,
#'               error_type = "classical")
make_inlami_matrices <- function(data,
                          formula_moi,
                          formula_imp,
                          error_type = "classical"){
  # Weaknesses:
  #  limited to the one classical, one berkson, one missing (same var) case, so model size must be 2+2

  # Set up sizes
  n <- nrow(data)
  n_models <- 2 + 2 # Num of error mods plus MOI and imp. mod.
  Y_ncol <- n_models
  Y_nrow <- n_models*n
  #cov_nrows <- n_models*n

  # Extract and group all variables from formulas:
  vars <- extract_variables_from_formula(formula_moi = formula_moi,
                                         formula_imp = formula_imp)

  if(vars$error_variable == "error_variable"){
    stop("Please name the variable with error something other than 'error_variable', as this name is used internally in the code and will lead to errors.")
  }

  # For all non-error covariates, the procedure is the same. So only the error-prone covariate needs to be treated separately.
  # Similar for the imputation model.

  # Set up the response matrix:
  Y <- matrix(NA, Y_nrow, Y_ncol)

  Y[1:n, 1] <- as.matrix(data[vars$response_moi])         # Regression model of interest response
  Y[n+(1:n), 2] <- rep(0, n)                              # Berkson error model response
  Y[2*n+(1:n), 3] <- as.matrix(data[vars$error_variable]) # Classical error model response
  Y[3*n+(1:n), 4] <- rep(0, n)                            # Imputation model response

  # Set up intercept vector for MOI:
  beta.0 <- c(rep(1, n), rep(NA, 3*n))

  # Set up vector for the variable with error/missingness, and corresponding
  # vectors to enable imputation:
  beta.error_variable <- paste0("beta.", vars$error_variable)
  assign(beta.error_variable, c(1:n, rep(NA, 3*n)))

  id.x <- c(rep(NA, n), 1:n, rep(NA, n), rep(NA, n))
  weight.x <- c(rep(NA, n), rep(-1, n), rep(NA, n), rep(NA, n))

  id.r <- c(rep(NA, n), 1:n, 1:n, 1:n)
  weight.r <- c(rep(NA, n), rep(1, n), rep(1, n), rep(-1, n))

  # Set up vectors for all error-free covariates with the correct names:
  cov_moi_names <- c()

  for(variable in vars$covariates_error_free){
    var_name <- paste0("beta.", variable)
    cov_moi_names <- c(cov_moi_names, var_name)
    assign(var_name, c(as.matrix(data[variable]), rep(NA, 3*n)))
  }

  # Intercept for imputation model:
  alpha.0 <- c(rep(NA, 3*n), rep(1, n))

  # Set up vectors for covariates in the imputation model
  cov_imp_names <- c()

  for(variable in vars$covariates_imp){
    var_name <- paste0("alpha.", variable)
    cov_imp_names <- c(cov_imp_names, var_name)
    assign(var_name, c(rep(NA, 3*n), as.matrix(data[variable])))
  }

  # How do I save the objects that have names based on user input into a list?
  # I generate a string that contains the code I want to run, and then I use
  # eval(parse(.)) to actually run that code. Since I do have the names of the
  # variables as strings, such a string will be easy to generate.

  # The names of the objects I defined explicitly above:
  fixed_objects <- "Y = Y, beta.0 = beta.0, id.x = id.x, weight.x = weight.x, id.r = id.r, weight.r = weight.r, alpha.0 = alpha.0, "

  # The names of the variables that change based on the covariate names:
  variable_objects <- paste0(c(beta.error_variable, cov_moi_names, cov_imp_names), " = ",
                             c(beta.error_variable, cov_moi_names, cov_imp_names), collapse = ", ")

  dd <- NULL # Assign NULL to dd, just to avoid notes when running checks.

  # A string containing the code needed to define a list of all the objects:
  define_dataframe <- paste0("dd <- list(", fixed_objects, variable_objects, ")")

  # Evaluate the string of code from above:
  eval(parse(text = define_dataframe))

  return(dd)
}

#' Make data stacks for joint model specification in INLA
#'
#' @param data A data frame with response and covariate variables for the main model and the imputation model.
#' @param formula_moi an object of class "formula", describing the main model to be fitted.
#' @param formula_imp an object of class "formula", describing the imputation model for the mismeasured and/or missing observations.
#' @param error_type Type of error (one or more of "classical", "berkson", "missing")
#'
#' @return An object of class inla.data.stack with data structured according to specified formulas and error models.
#' @export
#'
#' @examples
#' f_moi <- y ~ x + z
#' f_imp <- x ~ z
#' make_inlami_stacks(data = simple_data,
#'               formula_moi = f_moi,
#'               formula_imp = f_imp,
#'               error_type = "classical")
make_inlami_stacks <- function(data,
                               formula_moi,
                               formula_imp,
                               error_type = "classical"){

  # Extract and group all variables from formulas:
  vars <- extract_variables_from_formula(formula_moi = formula_moi,
                                         formula_imp = formula_imp)

  if(vars$error_variable == "error_variable"){
    stop("Please name the variable with error something other than 'error_variable', as this name is used internally in the code and will lead to errors.")
  }

  n <- nrow(data)


  # Model of interest stack  ---------------------------------------------------

  # Set up intercept vector for MOI:
  beta.0 <- rep(1, n)

  # Set up vector for variable with error
  beta.error_variable <- paste0("beta.", vars$error_variable)
  assign(beta.error_variable, 1:n)

  # Set up vectors for all error-free covariates with the correct names:
  cov_moi_names <- c()

  for(variable in vars$covariates_error_free){
    var_name <- paste0("beta.", variable)
    cov_moi_names <- c(cov_moi_names, var_name)
    assign(var_name, as.matrix(data[variable]))
  }

  # The names of the variables that change based on the covariate names:
  moi_effects <- paste0(c("beta.0", beta.error_variable, cov_moi_names), " = ",
                             c("beta.0", beta.error_variable, cov_moi_names), collapse = ", ")

  moi_effects_list <- NULL # Assign NULL to moi_effects_list, just to avoid notes when running checks.

  # A string containing the code needed to define a list of all the objects:
  moi_code <- paste0("moi_effects_list <- list(", moi_effects, ")")

  # Evaluate the string of code from above:
  eval(parse(text = moi_code))

  # Response
  if("berkson" %in% error_type){
    response_moi <- list(Y = cbind(as.matrix(data[vars$response_moi]), NA, NA, NA))
  }else{
    response_moi <- list(Y = cbind(as.matrix(data[vars$response_moi]), NA, NA))
  }

  # y = x + z...
  stk_moi <- INLA::inla.stack(data = response_moi,
                        A = list(1),
                        effects = list(moi_effects_list),
                        tag = "moi")


  # Berkson stack --------------------------------------------------------------

  if("berkson" %in% error_type){
    # 0 = -x_true + r + u_b
    stk_b <- INLA::inla.stack(data = list(Y = cbind(NA, rep(0, n), NA, NA)),
                              A = list(1),
                              effects = list(
                                list(id.x = 1:n,
                                     weight.x = -1,
                                     id.r = 1:n,
                                     weight.r = 1)),
                              tag = "berkson")
  }


  # Classical stack  -----------------------------------------------------------

  # Latent variable r if Berkson ME, otherwise x
  if("berkson" %in% error_type){
    classical_effects_list <- list(id.r = 1:n, weight.r = 1)
  }else{
    classical_effects_list <- list(id.x = 1:n, weight.x = 1)
  }

  # Response
  if("berkson" %in% error_type){
    response_classical <- list(Y = cbind(NA, NA, as.matrix(data[vars$error_variable]), NA))
  }else{
    response_classical <- list(Y = cbind(NA, as.matrix(data[vars$error_variable]), NA))
  }

  # x = r + u_c
  stk_c <- INLA::inla.stack(data = response_classical,
                      A = list(1),
                      effects = list(classical_effects_list),
                      tag = "classical")

  # Imputation stack  ----------------------------------------------------------

  # Intercept for imputation model:
  alpha.0 <- rep(1, n)

  # Set up vectors for covariates in the imputation model
  cov_imp_names <- c()

  for(variable in vars$covariates_imp){
    var_name <- paste0("alpha.", variable)
    cov_imp_names <- c(cov_imp_names, var_name)
    assign(var_name, as.matrix(data[variable]))
  }

  # The names of the variables that change based on the covariate names:
  imp_effects <- paste0(c("alpha.0", cov_imp_names), " = ",
                        c("alpha.0", cov_imp_names), collapse = ", ")

  imp_effects_list <- NULL # Assign NULL to imp_effects_list, just to avoid notes when running checks.

  # A string containing the code needed to define a list of all the objects:
  imp_code <- paste0("imp_effects_list <- list(", imp_effects, ")")

  # Evaluate the string of code from above:
  eval(parse(text = imp_code))

  # Latent variable r if Berkson ME, otherwise x
  if("berkson" %in% error_type){
    imp_effects_list$id.r <- 1:n
    imp_effects_list$weight.r <- rep(-1, n)
  }else{
    imp_effects_list$id.x <- 1:n
    imp_effects_list$weight.x <- rep(-1, n)
  }


  # Response
  if("berkson" %in% error_type){
    response_imputation <- list(Y = cbind(NA, NA, NA, rep(0, n)))
  }else{
    response_imputation <- list(Y = cbind(NA, NA, rep(0, n)))
  }

  # r = z + ...
  stk_imp <- INLA::inla.stack(data = response_imputation,
                        A = list(1),
                        effects = list(imp_effects_list),
                        tag = "imputation")

  # Join stacks ----

  if("berkson" %in% error_type){
    stk_full <- INLA::inla.stack(stk_moi, stk_b, stk_c, stk_imp)
  }else{
    stk_full <- INLA::inla.stack(stk_moi, stk_c, stk_imp)
  }

  return(stk_full)
}

#' Construct scaling vector to scale the precision of non-mismeasured observations
#'
#' @param data a data frame containing the variables in the model.
#' @param error_type one of "classical", "berkson" or "missing".
#' @param classical_error_scaling can be specified if the classical measurement error varies across observations. Must be a vector of the same length as the data.
#'
#' @return A vector reflecting the scaling factor for the residual terms in each model level.
#' @export
#'
#' @examples
#' make_inlami_scaling_vector(simple_data,
#'                            error_type = c("classical", "berkson"))
make_inlami_scaling_vector <- function(data,
                                error_type,
                                classical_error_scaling = NULL){
  # Scale the classical error precision of the perfectly measured values to a very high value (10^12).
  # If we have only missing data/perfectly observed data, then the precision can be scaled for all (because it makes no difference for the missing values)
  # but if we have measurement error (possibly varying), then the value of the precision is more meaningful.

  n <- nrow(data)

  # Default case:
  scale_classical <- rep(1, n)
  scale_berkson <- rep(1, n)

  if(!"classical" %in% error_type){
    # If there is no classical measurement error, then the classical error
    # precision should be scaled very high to indicate that there is no error.
    # Even if there is missingness, this scaling should be done, as the missing
    # values will anyway be imputed.
    scale_classical <- rep(10^12, n)
  }

  if(!"berkson" %in% error_type){
    # If there is no Berkson error, we scale this to a large value to "skip"
    # that level of the model.
    scale_berkson <- rep(10^12, n)
  }

  if(!is.null(classical_error_scaling)){
    if(length(classical_error_scaling) != n){
      stop(paste0("The length of classical_error_scaling (",
                  length(classical_error_scaling),
                  ") is not the same as the number of observations (",
                  n,
                  ")."))
    }

    # This can be used if the error varies for each observation somehow.
    scale_classical <- classical_error_scaling
  }

  precision_scaling <- c(rep(1, n), scale_berkson, scale_classical, rep(1, n))
  return(precision_scaling)
}

#' Fit model for measurement error and missing data in INLA
#'
#' A wrapper function around "INLA::inla()", providing the necessary structure to fit the hierarchical measurement error model that adjusts coefficient estimates to account for biases due to measurement error and missing data.
#'
#' @param formula_moi an object of class "formula", describing the main model to be fitted.
#' @param formula_imp an object of class "formula", describing the imputation model for the mismeasured and/or missing observations.
#' @param family_moi a string indicating the likelihood family for the model of interest (the main model).
#' @param data a data frame containing the variables in the model.
#' @param error_type type of error (one or more of "classical", "berkson", "missing")
#' @param classical_error_scaling can be specified if the classical measurement error varies across observations. Must be a vector of the same length as the data.
#' @param prior.prec.y a string containing the parameters for the prior for the precision of the residual term for the model of interest.
#' @param prior.prec.u_b a string containing the parameters for the prior for the precision of the error term for the Berkson error model.
#' @param prior.prec.u_c a string containing the parameters for the prior for the precision of the error term for the classical error model.
#' @param prior.prec.r a string containing the parameters for the precision of the latent variable r, which is the variable being described in the imputation model.
#' @param initial.prec.y the initial value for the precision of the residual term for the model of interest.
#' @param initial.prec.u_b the initial value for the precision of the residual term for the Berkson error term.
#' @param initial.prec.u_c the initial value for the precision of the residual term for the classical error term.
#' @param initial.prec.r the initial value for the precision of the residual term for the latent variable r.
#' @param control.family control.family for use in inla (can be provided directly instead of passing the "prior.prec...." and "initial.prec..." arguments.
#' @param control.predictor control.predictor for use in inla.
#' @param ... other arguments to pass to `inla`.
#'
#' @return An object of class \code{inlami}.
#' @export
#'
#' @examples
#' simple_moi <- y ~ x + z
#' simple_imp <- x ~ z
#'
#' # Fit the model
#' simple_model <- fit_inlami(data = simple_data,
#'                            formula_moi = simple_moi,
#'                            formula_imp = simple_imp,
#'                            family_moi = "gaussian",
#'                            error_type = c("berkson", "classical"),
#'                            prior.prec.y = c(0.5, 0.5),
#'                            prior.prec.u_b = c(10, 9),
#'                            prior.prec.u_c = c(10, 9),
#'                            prior.prec.r = c(0.5, 0.5),
#'                            initial.prec.y = 1,
#'                            initial.prec.u_b = 1,
#'                            initial.prec.u_c = 1,
#'                            initial.prec.r = 1)
fit_inlami <- function(formula_moi,
                     formula_imp,
                     family_moi,
                     data,
                     error_type = "classical",
                     classical_error_scaling = NULL,
                     prior.prec.y = NULL,
                     prior.prec.u_b = NULL,
                     prior.prec.u_c = NULL,
                     prior.prec.r = NULL,
                     initial.prec.y = NULL,
                     initial.prec.u_b = NULL,
                     initial.prec.u_c = NULL,
                     initial.prec.r = NULL,
                     control.family = NULL,
                     control.predictor = NULL,
                     ...){

  # Make formula from sub-models -----------------------------------------------
  formula_full <- make_inlami_formula(formula_moi = formula_moi,
                               formula_imp = formula_imp,
                               error_type = error_type)

  # Make scaling vector --------------------------------------------------------
  scaling_vec <- make_inlami_scaling_vector(
   data = data,
   error_type = error_type,
   classical_error_scaling = classical_error_scaling)

  # Make the stacks for the joint model --------------------------------------
  # data_matrices <- make_inlami_matrices(data = data,
  #                                       formula_moi = formula_moi,
  #                                       formula_imp = formula_imp,
  #                                       error_type = error_type)
  data_stack <- make_inlami_stacks(data = data,
                                   formula_moi = formula_moi,
                                   formula_imp = formula_imp,
                                   error_type = error_type)

  # Append scaling vector and data size to the matrices to pass to inla() ------
  n <- nrow(data) # Define number of observations
  data_for_inla <- INLA::inla.stack.data(data_stack, n = n, scaling_vec = scaling_vec)

  #data_with_other_stuff <- data_matrices
  #data_with_other_stuff$n <- n
  #data_with_other_stuff$scaling_vec <- scaling_vec

  # Make "control.family" ------------------------------------------------------
  if(is.null(control.family)){
    if(family_moi == "binomial"){
      control.family.y <- list(hyper = list())
    }
    if(family_moi == "gaussian"){
      control.family.y <- list(hyper = list(prec = list(initial = log(initial.prec.y),
                                                        param = prior.prec.y,
                                                        fixed = FALSE)))
    }
    control.family.u_b <- list(hyper = list(prec = list(initial = log(initial.prec.u_b),
                                                        param = prior.prec.u_b,
                                                        fixed = FALSE)))
    control.family.u_c <- list(hyper = list(prec = list(initial = log(initial.prec.u_c),
                                                        param = prior.prec.u_c,
                                                        fixed = FALSE)))
    control.family.r <- list(hyper = list(prec = list(initial = log(initial.prec.r),
                                                      param = prior.prec.r,
                                                      fixed = FALSE)))
    control.family <- list(
      control.family.y,
      control.family.u_b,
      control.family.u_c,
      control.family.r
    )
  }

  # Specify "control.predictor" ------------------------------------------------
  if(is.null(control.predictor)){
    control.predictor <- list(compute = TRUE)
  }

  # Run everything in the "inla"-function --------------------------------------
  inlami_model <- INLA::inla(
    formula = formula_full,
    family = c(family_moi, "gaussian", "gaussian", "gaussian"),
    data = data_for_inla,
    scale = scaling_vec,
    control.family = control.family,
    control.predictor = control.predictor,
    ...
  )

  # Set call to be full formula
  inlami_model$call <- formula_full

  # Add interesting arguments to output
  inlami_model$.args$formula_moi <- formula_moi
  inlami_model$.args$formula_imp <- formula_imp
  inlami_model$.args$error_type <- error_type

    # Set class ------------------------------------------------------------------
  class(inlami_model) <- c("inlami", class(inlami_model))


  return(inlami_model)
}
