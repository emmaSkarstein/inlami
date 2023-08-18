test_that("extract_variables_from_formula works", {
  simple_moi <- y ~ x + z
  simple_imp <- x ~ z
  expected <- extract_variables_from_formula(formula_moi = simple_moi, formula_imp = simple_imp)

  # Catching errors -----------------------------------------------------------
  expect_error(extract_variables_from_formula(formula_moi = "string", formula_imp = simple_imp))

})

test_that("make_inlami_matrices works", {
  simple_moi <- y ~ x + z
  simple_imp <- x ~ z
  matrices <- make_inlami_matrices(simple_data, simple_moi, simple_imp)
  # Check if the length of the list is 10
  expect_equal(length(matrices), 10)
  # Check if the matrices (list elements) have the correct names
  expected_names <- c("Y", "beta.0", "id.x", "weight.x", "id.r", "weight.r",
                      "alpha.0", "beta.x", "beta.z", "alpha.z")
  expect_equal(names(matrices), expected_names)

  # Catching errors -----------------------------------------------------------

  # Check if the error variable is called "error_variable"
  formula_moi <- mpg ~ error_variable + cyl + disp
  formula_imp <- error_variable ~ cyl
  expect_error(make_inlami_matrices(data = mtcars,
                             formula_moi = formula_moi,
                             formula_imp = formula_imp))


})

test_that("fit_inlami works", {
  simple_moi <- y ~ x + z
  simple_imp <- x ~ z

  # Priors for y, measurement error and true x-value precision
  prior.prec.y <- c(0.5, 0.5) # Gamma(0.5, 0.5)
  prior.prec.u_b <- c(10, 9) # Gamma(0.5, 0.5)
  prior.prec.u_c <- c(10, 9) # Gamma(0.5, 0.5)
  prior.prec.r <- c(0.5, 0.5) # Gamma(0.5, 0.5)

  # Initial values
  initial.prec.y <- 1
  initial.prec.u_b <- 1
  initial.prec.u_c <- 1
  initial.prec.r <- 1

  # Fit the model
  #simple_model <- fit_inlami(data = simple_data,
  #                         formula_moi = simple_moi,
  #                         formula_imp = simple_imp,
  #                         family_moi = "gaussian",
  #                         error_type = c("berkson", "classical"),
  #                         prior.prec.y = prior.prec.y,
  #                         prior.prec.u_b = prior.prec.u_b,
  #                         prior.prec.u_c = prior.prec.u_c,
  #                         prior.prec.r = prior.prec.r,
  #                         initial.prec.y = initial.prec.y,
  #                         initial.prec.u_b = initial.prec.u_b,
  #                         initial.prec.u_c = initial.prec.u_c,
  #                         initial.prec.r = initial.prec.r)

  # Test if the number of hyperparameters is 5
  #expect_equal(nrow(simple_model$summary.hyperpar), 5)

  # Catching errors -----------------------------------------------------------

})
