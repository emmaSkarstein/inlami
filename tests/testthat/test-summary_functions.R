test_that("simplify_inlami_model_summary works", {
  # CHANGE THIS!!
  #expect_equal(2 * 2, 4)
})

test_that("symmary.inlami works", {
  simple_moi <- y ~ x + z
  simple_imp <- x ~ z

  # Fit the model
  simple_model <- fit_inlami(data = simple_data,
                             formula_moi = simple_moi,
                             formula_imp = simple_imp,
                             family_moi = "gaussian",
                             error_type = c("berkson", "classical"),
                             prior.prec.y = c(0.5, 0.5),
                             prior.prec.u_b = c(10, 9),
                             prior.prec.u_c = c(10, 9),
                             prior.prec.r = c(0.5, 0.5),
                             initial.prec.y = 1,
                             initial.prec.u_b = 1,
                             initial.prec.u_c = 1,
                             initial.prec.r = 1)

  inlami_summary <- summary(simple_model)

  # Check summary class
  expect_s3_class(inlami_summary, "summary.inlami")

  # Catching errors -----------------------------------------------------------

})
