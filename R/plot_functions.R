
#' Plot model summary
#'
#' @param inlami_model the model returned from the fit_inlami function.
#' @param plot_moi should the posterior mean for the coefficients of the model of interest be plotted? Defaults to TRUE.
#' @param plot_imp should the posterior mean for the coefficients of the imputation model be plotted? Defaults to TRUE.
#' @param plot_intercepts should the posterior mean for the intercept(s) be plotted?
#'
#' @return An object of class "ggplot2" that plots the posterior mean and 95 \% credible interval for each coefficient in the model. The coefficients are colored to indicate if they belong to the main or imputation model, and the variable with error is also highlighted.
#' @export
#'
#' @importFrom ggplot2 ggplot aes vars
#' @importFrom rlang .data
#' @examples
#' #' #' simple_moi <- y ~ x + z
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
#' #simple_model <- fit_inlami(data = simple_data,
#' #                         formula_moi = simple_moi,
#' #                         formula_imp = simple_imp,
#' #                         family_moi = "gaussian",
#' #                         error_type = c("berkson", "classical"),
#' #                         prior.prec.y = prior.prec.y,
#' #                         prior.prec.u_b = prior.prec.u_b,
#' #                         prior.prec.u_c = prior.prec.u_c,
#' #                         prior.prec.r = prior.prec.r,
#' #                         initial.prec.y = initial.prec.y,
#' #                         initial.prec.u_b = initial.prec.u_b,
#' #                         initial.prec.u_c = initial.prec.u_c,
#' #                         initial.prec.r = initial.prec.r)
#'
#' #plot_inlami(simple_model)
plot_inlami <- function(inlami_model, plot_moi = TRUE, plot_imp = TRUE, plot_intercepts = TRUE){
  # Plot the posterior mean and 0.975 and 0.025 quantiles for all the coefficients
  # Highlight the error prone variable?
  # Return ggplot2 object so that it can be futher modified by user

  # What other arguments would be useful?
  # - Which covariates to plot
  # - Font? (check if this can be changed for the ggplot object after being returned)
  # - colors? (check if this can be changed for the ggplot object after being returned)

  simple_summary <- simplify_inlami_model_summary(inlami_model)

  joint_summary_df <- dplyr::bind_rows(moi_coef = simple_summary$moi_coef[,1:6],
                                       error_coef = simple_summary$error_coef[,1:6],
                                       imp_coef = simple_summary$imp_coef[,1:6],
                                       .id = "coefficient_type") |>
    dplyr::rename(quant_0.025 = "0.025quant", quant_0.975 = "0.975quant") |>
    dplyr::mutate(coefficient_type = dplyr::recode(.data$coefficient_type,
                                                   error_coef = "Variable w. error",
                                                   imp_coef = "Imputation model",
                                                   moi_coef = "Model of interest"))

  joint_summary_df$coefficient_name <- rownames(joint_summary_df)

  if(!plot_moi && !plot_imp){
    stop("This will plot nothing. Please set either 'plot_moi' or 'plot_imp' to TRUE.")
  }

  if(!plot_moi){
    # If plot_moi == FALSE, only plot the imputation model coefficients.
    joint_summary_df <- dplyr::filter(joint_summary_df, .data$coefficient_type == "Imputation model")
  }
  if(!plot_imp){
    # If plot_imp == FALSE, do not plot the imputation model coefficients.
    joint_summary_df <- dplyr::filter(joint_summary_df, .data$coefficient_type != "Imputation model")
  }
  if(!plot_intercepts){
    joint_summary_df <- dplyr::filter(joint_summary_df, !(.data$coefficient_name %in% c("beta.0", "alpha.0")))
  }

  color_pal <- c("#7AA6CA", "#DDAA33", "#EE849D")

  coef_plot <- ggplot(joint_summary_df, aes(y = .data$coefficient_name)) +
    # Error lines
    ggplot2::geom_linerange(aes(xmin = .data$quant_0.025,
                                xmax = .data$quant_0.975,
                                color = .data$coefficient_type),
                            linewidth = 1) +
    # Point for mean
    ggplot2::geom_point(aes(x = .data$mean,
                            color = .data$coefficient_type),
                        size = 3) +
    # Numeric text at mean
    ggplot2::geom_text(aes(x = .data$mean,
                           y = .data$coefficient_name,
                           label = round(.data$mean, digits = 3)),
                       vjust = -1, size = 3.5) +
    # Colors
    ggplot2::scale_color_manual(values = color_pal) +
    # Lables
    ggplot2::labs(x = "Posterior mean",
                  y = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = "bottom")

  coef_plot
}
