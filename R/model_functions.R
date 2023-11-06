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

  # Model of interest variables --
  moi_components <- as.list(formula_moi)

  ## Extract name of response:
  response_moi <- all.vars(moi_components[[2]])

  ## Check for random effects in moi:
  reff_moi_index <- which(grepl("f\\((.*)\\)", moi_components[[3]]))
  reff_vars_moi <- ""
  reff_moi <- ""
  if(length(reff_moi_index) > 0){
    reff_vars_moi <- setdiff(all.names(moi_components[[3]][reff_moi_index]), "f")
    reff_moi <- as.character(moi_components[[3]][reff_moi_index])
  }

  ## Covariates in moi (select variables that are NOT random effects)
  covariates_moi <- setdiff(all.vars(moi_components[[3]]), reff_vars_moi)

  ## Identify error prone variable:
  error_variable <- rlang::as_string(as.list(formula_imp)[[2]])

  ## Extract index of error variable in formula:
  error_var_index <- which(covariates_moi == error_variable)

  ## Error-free covariates:
  covariates_error_free <- covariates_moi[-error_var_index]

  # Imputation model variables --
  imp_components <- as.list(formula_imp)
  response_imp <- error_variable

  ## Check for random effects in imp:
  reff_imp_index <- which(grepl("f\\((.*)\\)", imp_components[[3]]))
  reff_imp <- ""
  reff_vars_imp <- ""
  if(length(reff_imp_index) > 0){
    reff_vars_imp <- setdiff(all.names(imp_components[[3]][reff_imp_index]), "f")
    reff_imp <- as.character(imp_components[[3]][reff_imp_index])
  }

  ## Imputation model covariates:
  covariates_imp <- setdiff(all.vars(imp_components[[3]]), reff_vars_imp)

  return(list(response_moi = response_moi,
              covariates_moi = covariates_moi,
              random_effects_moi = reff_moi,
              random_effect_variables_moi = reff_vars_moi,
              error_variable = error_variable,
              covariates_error_free = covariates_error_free,
              response_imp = response_imp,
              covariates_imp = covariates_imp,
              random_effects_imp = reff_imp,
              random_effect_variables_imp = reff_vars_imp))
}


#' Make formula for measurement error and missing data model
#'
#' @param formula_moi an object of class "formula", describing the main model to be fitted.
#' @param formula_imp an object of class "formula", describing the imputation model for the mismeasured and/or missing observations.
#' @param family_moi a string indicating the likelihood family for the model of interest (the main model).
#' @param error_type type of error (one of "classical", "berkson", "missing")
#' @param prior.beta.error parameters for the prior for the coefficient of the error prone variable. TODO: Which distribution? Gamma?
#'
#' @return An object of class "formula".
#' @export
#'
#' @examples
#' f_moi <- y ~ x + z
#' f_imp <- x ~ z
#' make_inlami_formula(formula_moi = f_moi,
#'                     formula_imp = f_imp,
#'                     error_type = "classical",
#'                     prior.beta.error = c(0, 1/1000)
#'                     )
make_inlami_formula <- function(formula_moi,
                                formula_imp,
                                family_moi = "gaussian",
                                error_type = "classical",
                                prior.beta.error){

  # Weaknesses:
  # Fixed size of id.x and id.r values (1:n)

  # Extract and group all variables from formulas:
  vars <- extract_variables_from_formula(formula_moi = formula_moi,
                                         formula_imp = formula_imp)

  # Define covariates in output formula as all variables in moi formula
  #  except the response and error prone covariate:
  covariates_error_free_string <- paste0("beta.", vars$covariates_error_free)

  # Covariates for imputation model:
  if(length(vars$covariates_imp) > 0){
    covariates_imp_string <- paste0("alpha.", vars$covariates_imp)
  }else{
    covariates_imp_string <- ""
  }


  # The copy term to ensure the mismeasured variable is copied correctly through the models
  #prior.beta.error <- c(0, 1/1000)
  copy_term1 <- paste0("f(", paste0("beta.", vars$response_imp),
                       ", copy = 'id.x', hyper = list(beta = list(param = ",
                       deparse(prior.beta.error), ", fixed = FALSE)))")
  # TODO: What does the "param = c(0, 1/1000)" above actually control? Should this be given as an argument?
  # I think it's the prior for beta.x?

  copy_term2 <- "f(id.x, weight.x, model='iid', values = 1:n, hyper = list(prec = list(initial = -15, fixed=TRUE)))"

  # Setting up the response formula
  if(!family_moi %in% inla_survival_families()){
    y_moi <- "y_moi, "
  }else{
    y_moi <- "inla.surv(time = y_time, event = y_event), "
  }
  if("berkson" %in% error_type){
    y_berkson <- "y_berkson, "
  }else{
    y_berkson <- ""
  }

  response_list <- paste0("list(", y_moi, y_berkson, "y_classical, ", "y_imp)")

  formula_RHS <- paste(
    # Remove default intercept
    "-1",
    # MOI covariates
    "beta.0",
    paste(covariates_error_free_string, collapse = " + "),
    # MOI random effects
    paste(vars$random_effects_moi, collapse = " + "),
    # Imputation model covariates
    "alpha.0",
    paste(covariates_imp_string, collapse = " + "),
    # Imputation model random effects
    paste(vars$random_effects_imp, collapse = " + "),
    # Copy terms
    copy_term1,
    copy_term2,
    # Between each term add a "+"
    sep = " + ")

  # Remove duplicate "+"'s
  formula_RHS <- gsub(pattern = "((\\+)( *)){2,}", replacement = "+ ",
                      formula_RHS)

  # Optionally add Berkson layer
  if("berkson" %in% error_type){
    r_term <- "f(id.r, weight.r, model='iid', values = 1:n, hyper = list(prec = list(initial = -15, fixed = TRUE)))"
    formula_RHS <- paste(formula_RHS, "+", r_term)
  }


  formula <- paste(response_list, "~", formula_RHS)

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
#' @param formula_moi an object of class "formula", describing the main model to be fitted.
#' @param formula_imp an object of class "formula", describing the imputation model for the mismeasured and/or missing observations.
#' @param family_moi a string indicating the likelihood family for the model of interest (the main model).
#' @param data A data frame with response and covariate variables for the main model and the imputation model.
#' @param error_type Type of error (one or more of "classical", "berkson", "missing")
#' @param repeated_observations Does the variable with measurement error and/or missingness have repeated observations? If so, set this to "TRUE". In that case, when specifying the formula, use the name of the variable without any numbers, but when specifying the data, make sure that the repeated measurements end in a number, i.e "sbp1" and "sbp2".
#'
#' @return An object of class inla.stack with data structured according to specified formulas and error models.
#' @export
#'
#' @examples
#' make_inlami_stacks(formula_moi = y ~ x + z,
#'                    formula_imp = x ~ z,
#'                    data = simple_data,
#'                    error_type = "classical")
make_inlami_stacks <- function(formula_moi,
                               formula_imp,
                               family_moi = "gaussian",
                               data,
                               error_type = "classical",
                               repeated_observations = FALSE){

  # Extract and group all variables from formulas:
  vars <- extract_variables_from_formula(formula_moi = formula_moi,
                                         formula_imp = formula_imp)

  # Some checks ----
  if(vars$error_variable == "error_variable"){
    stop("Please name the variable with error something other than 'error_variable', as this name is used internally in the code and will lead to errors.")
  }

  # Check that variables are named correctly in case of repeated measurements
  if(repeated_observations){
    expected_data_variable_names <- paste0(vars$error_variable, c(1,2))
    if(!all(expected_data_variable_names %in% names(data))){
      stop("It looks like you have repeated measurements, and that the formula or data may not be specified correctly. When fitting the model with repeated observations for the variable with missingness or measurement error, make sure that you specify the variable with its general name in the formula, but provide numbered versions in the data, such that each repeat is labelled with a number after the variable name. For instance, assume we have repeated measurements of systolic blood pressure (sbp). Then we would write 'sbp' in the formula, but in the data they would be labelled 'sbp1' and 'sbp2'.")
    }
  }

  # Check if it looks like there may be repeated measurements
  repeated_obs_warning <- FALSE
  if(sum(grepl(paste0("\\b", vars$error_variable, "(\\d)\\b"),
               names(data)))>1 &&
     !repeated_observations){
    repeated_obs_warning <- TRUE
    warning("It looks like you may have repeated measurements of the variable with measurement error or missingness. If this is the case, specify the argument 'repeated_observations = TRUE' to ensure correct analysis.")
  }

  # Check that all variables in formula exist in data
  formula_vars <- union(all.vars(formula_moi), all.vars(formula_imp))
  data_vars <- names(data)
  if(repeated_observations | repeated_obs_warning){
    # If reason to believe repeated measurements, ignore the error variable in the comparison
    formula_vars <- setdiff(formula_vars, vars$error_variable)
  }
  diff_vars <- setdiff(formula_vars, data_vars) # vars in formula but not in data
  if(length(diff_vars)>0){
    stop(paste0("It looks like the following variable(s) are in the formula but not the data: ", toString(diff_vars)))
  }


  # Defining sizes -------------------------------------------------------------
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

  # Add (optional) random effect variable
  if(nchar(vars$random_effects_moi) > 0){
    # "list", "-" and "c" will be recognized as variables... so remove them
    reff_moi_names <- setdiff(vars$random_effect_variables_moi, c("list", "-", "c"))
    assign(reff_moi_names, as.matrix(data[reff_moi_names]))

    reff_list_moi <- paste0(reff_moi_names, " = ", reff_moi_names)
    moi_effects <- paste0(moi_effects, ", ", reff_list_moi)
  }

  moi_effects_list <- NULL # Assign NULL to moi_effects_list, just to avoid notes when running checks.

  # A string containing the code needed to define a list of all the objects:
  moi_code <- paste0("moi_effects_list <- list(", moi_effects, ")")

  # Evaluate the string of code from above:
  eval(parse(text = moi_code))

  # Response
  ## if survival model, structure the response accordingly
  if(family_moi %in% inla_survival_families()){
    time_var <- vars$response_moi[1]
    event_var <- vars$response_moi[2]
    response_moi <- list(y_time = as.matrix(data[time_var]),
                         y_event = as.matrix(data[event_var]))
  }else{ ## otherwise response is just one vector
    response_moi <- list(y_moi = as.matrix(data[vars$response_moi]))
  }
  # y = x + z...
  stk_moi <- INLA::inla.stack(data = response_moi,
                        A = list(1),
                        effects = list(moi_effects_list),
                        tag = "moi")


  # Berkson stack --------------------------------------------------------------

  if("berkson" %in% error_type){
    # 0 = -x_true + r + u_b
    response_berkson <- list(y_berkson = rep(0, n))

    stk_b <- INLA::inla.stack(data = response_berkson,
                              A = list(1),
                              effects = list(
                                list(id.x = 1:n,
                                     weight.x = -1,
                                     id.r = 1:n,
                                     weight.r = 1)),
                              tag = "berkson")
  }


  # Classical stack  -----------------------------------------------------------

  # Check if repeated measurements
  error_var_list <- names(data)[grepl(paste0("\\b", vars$error_variable, "(\\d|)\\b"), names(data))]
  nr_repeats <- length(error_var_list)
  classical_id <- rep(1:n, nr_repeats)

  # Latent variable r if Berkson ME, otherwise x
  if("berkson" %in% error_type){
    classical_effects_list <- list(id.r = classical_id, weight.r = 1)
  }else{
    classical_effects_list <- list(id.x = classical_id, weight.x = 1)
  }

  # Response
  classical_response <- as.matrix(utils::stack(data[error_var_list])[1])

  response_classical <- list(y_classical = classical_response)

  # x = r + u_c
  stk_c <- INLA::inla.stack(data = response_classical,
                      A = list(1),
                      effects = list(classical_effects_list),
                      compress = FALSE,
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

  # Add (optional) random effect variable
  if(nchar(vars$random_effects_imp) > 0){
    reff_imp_names <- setdiff(vars$random_effect_variables_imp, c("list", "-", "c"))
    assign(reff_imp_names, as.matrix(data[reff_imp_names]))

    reff_list_imp <- paste0(reff_imp_names, " = ", reff_imp_names)
    imp_effects <- paste0(imp_effects, ", ", reff_list_imp)
  }

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
  response_imputation <- list(y_imp = rep(0, n))

  # r = z + ...
  stk_imp <- INLA::inla.stack(data = response_imputation,
                        A = list(1),
                        effects = list(imp_effects_list),
                        tag = "imputation")

  # Join stacks ----

  if("berkson" %in% error_type){
    stk_full <- INLA::inla.stack(stk_moi, stk_b, stk_c, stk_imp, compress = FALSE)
  }else{
    stk_full <- INLA::inla.stack(stk_moi, stk_c, stk_imp, compress = FALSE)
  }

  return(stk_full)
}

#' Construct scaling vector to scale the precision of non-mismeasured observations
#'
#' @param inlami_stack an object of class inlami.data.stack containing data structured
#' @param error_type one of "classical", "berkson" or "missing".
#' @param classical_error_scaling can be specified if the classical measurement error varies across observations. Must be a vector of the same length as the data.
#'
#' @return A vector reflecting the scaling factor for the residual terms in each model level.
#' @export
#'
#' @examples
#' stacks <- make_inlami_stacks(data = simple_data,
#'                              formula_moi = y ~ x + z,
#'                              formula_imp = x ~ z,
#'                              error_type = c("classical", "berkson"))
#' make_inlami_scaling_vector(stacks,
#'                            error_type = c("classical", "berkson"))
make_inlami_scaling_vector <- function(inlami_stack,
                                       error_type,
                                       classical_error_scaling = NULL){
  # Scale the classical error precision of the perfectly measured values to a very high value (10^12).
  # If we have only missing data/perfectly observed data, then the precision can be scaled for all (because it makes no difference for the missing values)
  # but if we have measurement error (possibly varying), then the value of the precision is more meaningful.

  index_moi <- inlami_stack$data$index$moi
  index_berkson <- inlami_stack$data$index$berkson
  index_classical <- inlami_stack$data$index$classical
  index_imputation <- inlami_stack$data$index$imputation

  n_moi <- length(index_moi)
  n_berkson <- length(index_berkson)
  n_classical <- length(index_classical)
  n_imputation <- length(index_imputation)

  # Default case:
  scale_berkson <- rep(1, n_berkson)
  scale_classical <- rep(1, n_classical)

  if(!"classical" %in% error_type){
    # If there is no classical measurement error, then the classical error
    # precision should be scaled very high to indicate that there is no error.
    # Even if there is missingness, this scaling should be done, as the missing
    # values will anyway be imputed.
    scale_classical <- rep(10^12, n_classical)
  }

  if(!"berkson" %in% error_type){
    # If there is no Berkson error, we scale this to a large value to "skip"
    # that level of the model.
    scale_berkson <- rep(10^12, n_berkson)
  }

  if(!is.null(classical_error_scaling)){
    if(length(classical_error_scaling) != n_classical){
      stop(paste0("The length of classical_error_scaling (",
                  length(classical_error_scaling),
                  ") is not the same as the number of observations (",
                  n_classical,
                  ")."))
    }

    # This can be used if the error varies for each observation somehow.
    scale_classical <- classical_error_scaling
  }

  precision_scaling <- c(rep(1, n_moi), scale_berkson, scale_classical, rep(1, n_imputation))
  return(precision_scaling)
}

#' Fit model for measurement error and missing data in INLA
#'
#' A wrapper function around "INLA::inla()", providing the necessary structure to fit the hierarchical measurement error model that adjusts coefficient estimates to account for biases due to measurement error and missing data.
#'
#' @param formula_moi an object of class "formula", describing the main model to be fitted.
#' @param formula_imp an object of class "formula", describing the imputation model for the mismeasured and/or missing observations.
#' @param family_moi a string indicating the likelihood family for the model of interest (the main model).
#' @param data an object of class data.frame or list containing the variables in the model.
#' @param error_type type of error (one or more of "classical", "berkson", "missing")
#' @param repeated_observations Does the variable with measurement error and/or missingness have repeated observations? If so, set this to "TRUE". In that case, when specifying the formula, use the name of the variable without any numbers, but when specifying the data, make sure that the repeated measurements end in a number, i.e "sbp1" and "sbp2".
#' @param classical_error_scaling can be specified if the classical measurement error varies across observations. Must be a vector of the same length as the data.
#' @param prior.prec.moi a string containing the parameters for the prior for the precision of the residual term for the model of interest.
#' @param prior.prec.berkson a string containing the parameters for the prior for the precision of the error term for the Berkson error model.
#' @param prior.prec.classical a string containing the parameters for the prior for the precision of the error term for the classical error model.
#' @param prior.prec.imp a string containing the parameters for the precision of the latent variable r, which is the variable being described in the imputation model.
#' @param prior.beta.error parameters for the prior for the coefficient of the error prone variable. TODO: Which distribution? Gamma?
#' @param initial.prec.moi the initial value for the precision of the residual term for the model of interest.
#' @param initial.prec.berkson the initial value for the precision of the residual term for the Berkson error term.
#' @param initial.prec.classical the initial value for the precision of the residual term for the classical error term.
#' @param initial.prec.imp the initial value for the precision of the residual term for the latent variable r.
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
#'                            prior.prec.moi = c(10, 9),
#'                            prior.prec.berkson = c(10, 9),
#'                            prior.prec.classical = c(10, 9),
#'                            prior.prec.imp = c(10, 9),
#'                            prior.beta.error = c(0, 1/1000),
#'                            initial.prec.moi = 1,
#'                            initial.prec.berkson = 1,
#'                            initial.prec.classical = 1,
#'                            initial.prec.imp = 1)
fit_inlami <- function(formula_moi,
                     formula_imp,
                     family_moi,
                     data,
                     error_type = "classical",
                     repeated_observations = FALSE,
                     classical_error_scaling = NULL,
                     prior.prec.moi = NULL,
                     prior.prec.berkson = NULL,
                     prior.prec.classical = NULL,
                     prior.prec.imp = NULL,
                     prior.beta.error = NULL,
                     initial.prec.moi = NULL,
                     initial.prec.berkson = NULL,
                     initial.prec.classical = NULL,
                     initial.prec.imp = NULL,
                     control.family = NULL,
                     control.predictor = NULL,
                     ...){

  # Some checks ----------------------------------------------------------------


  # Make formula from sub-models -----------------------------------------------
  formula_full <- make_inlami_formula(formula_moi = formula_moi,
                               formula_imp = formula_imp,
                               family_moi = family_moi,
                               error_type = error_type,
                               prior.beta.error = prior.beta.error)

  # Define some stuff ----------------------------------------------------------
  model_levels_without_moi <- union(error_type, c("classical", "imp")) # Always need classical
  model_levels_without_moi_and_missing <- setdiff(model_levels_without_moi, "missing")

  # Make the stacks for the joint model ----------------------------------------
  data_stack <- make_inlami_stacks(formula_moi = formula_moi,
                                   formula_imp = formula_imp,
                                   family_moi = family_moi,
                                   data = data,
                                   error_type = error_type,
                                   repeated_observations = repeated_observations)

  # Make scaling vector --------------------------------------------------------
  scaling_vec <- make_inlami_scaling_vector(
    inlami_stack = data_stack,
    error_type = error_type,
    classical_error_scaling = classical_error_scaling)

  # Append scaling vector and data size to the matrices to pass to inla() ------
  n <- length(data_stack$data$index$moi) # Define number of observations
  data_for_inla <- INLA::inla.stack.data(data_stack,
                                         n = n,
                                         scaling_vec = scaling_vec,
                                         compress = FALSE)


  # Make "control.family" ------------------------------------------------------
  # If control.family is specified, just use that directly. If not:
  if(is.null(control.family)){
    # Check if any defaults are needed:
    prec.list <- list(moi = prior.prec.moi, berkson = prior.prec.berkson,
                      classical = prior.prec.classical, imp = prior.prec.imp)
    # Not so easy to find indexes of NULL elements in list.
    # Workaround is to find elements that have length 0.
    which_prec_null <- names(prec.list[lengths(prec.list)==0])

    # Find out which of these are missing:
    which_need_defaults <- intersect(which_prec_null, model_levels_without_moi)

    # Set defaults for necessary priors that have not been specified
    for(component in which_need_defaults){
      #warning(paste0("Using default prior Gamma(10, 9) for the precision of the ",
      #               component, " model, and the initial value is set to 1 for the same term. This should be given a better value by specifying 'prior.prec.",
      #               component, "' and 'initial.prec.", component, "'.\n"))
      prior_code <- paste0("prior.prec.", component, " <- c(10, 9); ",
                           "initial.prec.", component, " <- 1")
      eval(parse(text = prior_code))
    }

    # Define control.family.moi
    if(family_moi == "binomial"){
      control.family.moi <- list(hyper = list())
    }
    if(family_moi == "gaussian"){
      # Check if default is needed
      if(any(is.null(initial.prec.moi), is.null(prior.prec.moi))){
        warning("Using default prior Gamma(10, 9) for the precision of the model of interest, and the initial value is set to 1 for the same term. This should be given a better value by specifying 'prior.prec.moi' and 'initial.prec.moi'.\n")
        prior.prec.moi <- c(10, 9)
        initial.prec.moi <- 1
      }
      control.family.moi <- list(hyper = list(prec = list(initial = log(initial.prec.moi),
                                                        param = prior.prec.moi,
                                                        fixed = FALSE)))
    }
    # What if poisson or survival?

    control.family <- list(control.family.moi)

    # Optionally define control.family.berkson
    if("berkson" %in% error_type){
      control.family.berkson <- list(hyper = list(
        prec = list(initial = log(initial.prec.berkson),
                    param = prior.prec.berkson,
                    fixed = FALSE)))
      control.family <- append(control.family, list(control.family.berkson))
    }

    # Define control.family.classical and control.family.imp
    control.family.classical <- list(hyper = list(
      prec = list(initial = log(initial.prec.classical),
                  param = prior.prec.classical,
                  fixed = FALSE)))
    control.family <- append(control.family, list(control.family.classical))

    control.family.imp <- list(hyper = list(
      prec = list(initial = log(initial.prec.imp),
                  param = prior.prec.imp,
                  fixed = FALSE)))
    control.family <- append(control.family, list(control.family.imp))
  }

  # Specify "control.predictor" ------------------------------------------------
  if(is.null(control.predictor)){
    control.predictor <- list(compute = TRUE)
  }

  # Run everything in the "inla"-function --------------------------------------
  inlami_model <- INLA::inla(
    formula = formula_full,
    family = c(family_moi, rep("gaussian", length(model_levels_without_moi_and_missing))),
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
  inlami_model$stack_data <- data_stack

    # Set class ------------------------------------------------------------------
  class(inlami_model) <- c("inlami", class(inlami_model))


  return(inlami_model)
}
