test_that("extract_variables_from_formula works", {
  # Check that all elements are characters
  expected <- extract_variables_from_formula(formula_moi = y ~ x + z,
                                             formula_imp = x ~ z)
  expect_true(all(sapply(expected, class) == "character"))

  # Check that it correctly identifies interaction between error variable and another covariate.
  expected2 <- extract_variables_from_formula(
    formula_moi = y ~ x + x:z + x:m + z:s,
    formula_imp = x ~ z)
  expect_equal(expected2[["error_interaction_variables"]], c("z", "m"))

  # Make sure it also works if only term is interaction term
  expected3 <- extract_variables_from_formula(formula_moi <- y ~ smoking:sbp,
                                              formula_imp = sbp ~ smoking)
  expect_equal(expected3[["error_interaction_variables"]], "smoking")

  # Catching errors -----------------------------------------------------------
  expect_error(extract_variables_from_formula(formula_moi = "string", formula_imp = simple_imp))

})

test_that("make_inlami_stacks works", {
  # Some checks with simulated data, stacks with berkson level
  simple_moi <- y ~ x + z
  simple_imp <- x ~ z
  simple_stacks_cb <- make_inlami_stacks(formula_moi = simple_moi,
                                         formula_imp = simple_imp,
                                         data = simple_data,
                                         error_type = c("classical", "berkson"))

  # Check rows and columns of "data"
  expect_equal(nrow(simple_stacks_cb$data$data), 4*1000)
  expect_equal(ncol(simple_stacks_cb$data$data), 4)

  # Check rows and columns of "effects"
  expect_equal(nrow(simple_stacks_cb$effects$data), 4*1000)
  expect_equal(ncol(simple_stacks_cb$effects$data), 9)

  # Check the names
  expected_names_effects <- c("beta.0", "id.x", "weight.x", "id.r", "weight.r",
                      "alpha.0", "beta.x", "beta.z", "alpha.z")
  actual_names <- unlist(simple_stacks_cb$effects$names)
  names(actual_names) <- NULL
  expect_equal(sort(actual_names), sort(expected_names_effects))

  # Some checks with framingham data
  framingham_stacks <- make_inlami_stacks(data = framingham,
                                          formula_moi = disease ~ sbp + smoking,
                                          formula_imp = sbp ~ smoking,
                                          error_type = "classical",
                                          repeated_observations = TRUE)
  # Check rows and columns of "data"
  expect_equal(nrow(framingham_stacks$data$data), 4*641)
  expect_equal(ncol(framingham_stacks$data$data), 3)

  # Check rows and columns of "effects"
  expect_equal(nrow(framingham_stacks$effects$data), 4*641)
  expect_equal(ncol(framingham_stacks$effects$data), 7)

  # Catching errors -----------------------------------------------------------

  # Check if the error variable is called "error_variable"
  formula_moi <- mpg ~ error_variable + cyl + disp
  formula_imp <- error_variable ~ cyl
  expect_error(make_inlami_stacks(data = mtcars,
                             formula_moi = formula_moi,
                             formula_imp = formula_imp))

  # Check if repeated observations are specified
  expect_error(make_inlami_stacks(data = framingham,
                                 formula_moi = disease ~ sbp + smoking,
                                 formula_imp = sbp ~ smoking,
                                 error_type = "classical"))

  # Check if variables in formula are in data (but not sensitive to repeated measurements)
  expect_error(make_inlami_stacks(
    formula_moi = inla.surv(t, d) ~ sbp + aged + smoking + sex + diabetes,
    formula_imp = sbp ~ age + smoke + sex + diabetes,
    data = nhanes_survival,
    error_type = c("classical", "missing"),
    repeated_observations = TRUE)
  )

})

test_that("make_inlami_scaling_vector works", {
  # Check that the scaling vector is the correct length
  expect_error(
    make_inlami_scaling_vector(simple_data,
                               error_variable = "x",
                               error_type = c("classical", "berkson"),
                               classical_error_scaling = rep(1, 500)))
})

test_that("make_inlami_control.family works", {
  cont.fam <- make_inlami_control.family(
    formula_moi = simple_moi,
    formula_imp = simple_imp,
    family_moi = "gaussian",
    error_type = c("berkson", "classical"),
    prior.prec.moi = c(10, 9),
    prior.prec.berkson = c(10, 9),
    prior.prec.classical = c(10, 9),
    prior.prec.imp = c(10, 9),
    prior.beta.error = c(0, 1/1000),
    initial.prec.moi = 1,
    initial.prec.berkson = 1,
    initial.prec.classical = 1,
    initial.prec.imp = 1)
  expect_equal(length(cont.fam), 4)
  expect_equal(cont.fam[[1]]$hyper$prec$initial, log(1))
  expect_equal(cont.fam[[1]]$hyper$prec$param, c(10, 9))
})

test_that("fit_inlami works", {
  simple_moi <- y ~ x + z
  simple_imp <- x ~ z

  # Fit the model
  simple_model <- fit_inlami(data = simple_data,
                             formula_moi = simple_moi,
                             formula_imp = simple_imp,
                             family_moi = "gaussian",
                             error_type = c("berkson", "classical"),
                             prior.prec.moi = c(10, 9),
                             prior.prec.berkson = c(10, 9),
                             prior.prec.classical = c(10, 9),
                             prior.prec.imp = c(10, 9),
                             prior.beta.error = c(0, 1/1000),
                             initial.prec.moi = 1,
                             initial.prec.berkson = 1,
                             initial.prec.classical = 1,
                             initial.prec.imp = 1)

  # Test if the number of hyperparameters is 5
  expect_equal(nrow(simple_model$summary.hyperpar), 5)

  # Check that the model hass the correct class
  expect_s3_class(simple_model, c("inlami", "inla"))

  # Check the estimate for beta.x and the precision of MOI-level
  beta.x <- simple_model$summary.hyperpar["Beta for beta.x", "mean"]
  expect_equal(beta.x, 2, tolerance = 0.3)
  prec.moi <- simple_model$summary.hyperpar["Precision for the Gaussian observations", "mean"]
  expect_equal(prec.moi, 1, tolerance = 0.2)

  # Catching errors -----------------------------------------------------------

})
