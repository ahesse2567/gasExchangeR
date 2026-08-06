#' Find the gas exchange threshold (GET) (~1st ventilatory threshold (VT1)) using the 1st derivative crossing point
#'
#' Use polynomial regression to find the gas exchange threshold (GET) (practically equivalent to the first ventilatory threshold (VT1)). This method by Wisen and Wohlfart (2004) is similar to the V-slope method by Beaver et al. (1986) in that it finds when the rate of CO2 production outpaces the rate of O2 production. In contrast to piecewise regression in the original V-slope method, this uses the first derivative of VO2 and VCO2 vs. time (differentiated analytically). The threshold is when dVCO2/dt surpasses dVO2/dt for the final time.
#'
#' By default a fixed sixth-degree polynomial is fit to each of VO2 and VCO2 vs. time. \strong{Note:} this degree-6 default differs from the degree-5 default used by the other polynomial and spline methods in this package (e.g. \code{d2_poly_reg_maxima}, \code{d2_reg_spline_maxima}). We use 6 here because the original Wisen and Wohlfart (2004) technique fit sixth-degree polynomials. Set \code{search_degree = TRUE} to instead select the degree from the data (between \code{degree} and \code{max_degree}) using the corrected AIC (\code{"aicc"}, default) or a sequential likelihood-ratio test (\code{"lrt"}).
#'
#' @param .data Gas exchange data frame or tibble.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param degree The polynomial degree fit to VO2 and VCO2 vs. time. Default \code{6}, matching Wisen and Wohlfart (2004) (note this differs from the degree-5 default of the other methods). When \code{search_degree = TRUE}, this is the lowest degree searched. Must be at least 5.
#' @param search_degree Search for the best polynomial degree rather than using the fixed \code{degree}? Default \code{FALSE}.
#' @param max_degree Highest polynomial degree considered when \code{search_degree = TRUE}. Default \code{8}.
#' @param criterion Model-selection rule when \code{search_degree = TRUE}: \code{"aicc"} (default) or \code{"lrt"}.
#' @param vo2 Name of the \code{vo2} variable
#' @param vco2 Name of the \code{vco2} variable
#' @param ve Name of the \code{ve} variable
#' @param time Name of the \code{time} variable
#' @param alpha_linearity Significance value used by the \code{"lrt"} degree search.
#' @param bp Is this the first (\code{vt1}) or the second (\code{vt2}) breakpoint? This method is specifically designed for finding VT1 (GET)
#' @param pos_change Do you expect the change in slope to be positive (default) or negative? If a two-line regression explains significantly reduces the sum square error but the change in slope does not match the expected underlying physiology, the breakpoint will be classified as indeterminate.
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param front_trim_vt1 Seconds of data to remove from the beginning of the test before fitting (VT1). Default \code{60}.
#' @param front_trim_vt2 Seconds of data to remove from the beginning of the test before fitting (VT2). Default \code{60}.
#' @param ci Should the output include confidence interval data? Default is `FALSE`.
#' @param conf_level Confidence level to use if calculating confidence intervals.
#' @param plots Should this function generate plots? Set to `FALSE` to save time.
#'
#' @returns A list include breakpoint data, best-fit and related functions, and a plot.
#'
#' @export
#'
#' @references
#' Beaver, W. L., Wasserman, K. A. R. L. M. A. N., & Whipp, B. J. (1986). A new method for detecting anaerobic threshold by gas exchange. Journal of applied physiology, 60(6), 2020-2027.
#' Wisen, A. G., & Wohlfart, B. (2004). A refined technique for determining the respiratory gas exchange responses to anaerobic metabolism during progressive exercise–repeatability in a group of healthy men. Clinical physiology and functional imaging, 24(1), 1-9.
#'
#' @examples
#'
#' cpet_bbb <- utils::read.csv(
#'   system.file("extdata", "anton_vo2max_clean.csv", package = "gasExchangeR")
#' )
#'
#' cpet_15b_roll <- avg_exercise_test(
#'   cpet_bbb,
#'   method = "breath", calc_type = "rolling", time_col = "time", window = 15
#' )
#'
#' # GET/VT1 as the final crossing of dVCO2/dt above dVO2/dt. The crossing is
#' # located in TIME, so .x must be "time" (thr_y is the fitted VO2, so .y = "vo2").
#' d1_crossing(
#'   cpet_15b_roll,
#'   .x = "time",
#'   .y = "vo2",
#'   bp = "vt1",
#'   vo2 = "vo2",
#'   vco2 = "vco2",
#'   ve = "ve",
#'   time = "time",
#'   plots = FALSE
#' )
#'
d1_crossing <- function(.data,
                        .x,
                        .y,
                        bp = "vt1",
                        ...,
                        degree = 6,
                        search_degree = FALSE,
                        max_degree = 8,
                        criterion = c("aicc", "lrt"),
                        vo2 = "vo2",
                        vco2 = "vco2",
                        ve = "ve",
                        time = "time",
                        alpha_linearity = 0.05, # change to just alpha?
                        pos_change = TRUE,
                        ordering = c("by_x", "time"),
                        front_trim_vt1 = 60,
                        front_trim_vt2 = 60,
                        ci = FALSE,
                        conf_level = 0.95,
                        plots = TRUE) {
  # check for crucial missing data
  rlang::check_required(.data)
  rlang::check_required(.x)
  rlang::check_required(.y)
  rlang::check_required(bp)

  bp <- rlang::arg_match(bp, values = c("vt1", "vt2"))
  ordering <- rlang::arg_match(ordering)
  criterion <- rlang::arg_match(criterion)
  .data <- order_cpet_df(.data,
    .x = .x, time = time,
    ordering = ordering
  )

  plot_df <- .data

  # remove the noisy start of the test before fitting (as the other threshold
  # functions do). Without this, spurious early derivative crossings appear.
  front_trim <- set_front_trim(
    bp = bp,
    front_trim_vt1 = front_trim_vt1,
    front_trim_vt2 = front_trim_vt2
  )
  .data <- .data %>%
    dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))

  # polynomial fit for vo2 (shared poly-fitting loop)
  lm_poly_vo2 <- fit_poly_reg(
    .data = .data, .x = time, .y = vo2,
    degree = degree, search_degree = search_degree,
    max_degree = max_degree, criterion = criterion,
    alpha_linearity = alpha_linearity
  )

  # polynomial fit for vco2
  lm_poly_vco2 <- fit_poly_reg(
    .data = .data, .x = time, .y = vco2,
    degree = degree, search_degree = search_degree,
    max_degree = max_degree, criterion = criterion,
    alpha_linearity = alpha_linearity
  )

  # if creating the linear model fails, return a quick summary
  if (is.null(lm_poly_vo2) |
    is.null(lm_poly_vco2) |
    any(
      is.na(lm_poly_vo2$coefficients),
      is.na(lm_poly_vco2$coefficients)
    )) {
    bp_dat <- return_indeterminant_findings(
      .data = plot_df, # change to plot_data later?
      bp = bp,
      algorithm = as.character(match.call()[[1]]),
      .x = .x,
      .y = .y,
      est_ci = "estimate"
    )

    return(list(breakpoint_data = bp_dat))
  }


  # 1st derivative for vo2
  poly_expr_vo2 <- expr_from_coefs(lm_poly_vo2$coefficients)
  deriv1_vo2 <- Deriv::Deriv(poly_expr_vo2, x = "x", nderiv = 1) # slope

  # 1st derivative for vco2
  poly_expr_vco2 <- expr_from_coefs(lm_poly_vco2$coefficients)
  deriv1_vco2 <- Deriv::Deriv(poly_expr_vco2, x = "x", nderiv = 1) # slope

  # turn derivative expressions into functions to use with uniroot.all()
  deriv1_vo2_func <- expr_to_func(deriv1_vo2)
  deriv1_vco2_func <- expr_to_func(deriv1_vco2)

  # calculate derivative of difference in 1st derivative functions
  # this lets us know if vco2 was surpassing vo2 or the other way around
  deriv_diff_deriv1 <- paste0(
    deriv1_vo2 %>%
      Ryacas::y_fn("Simplify") %>%
      Ryacas::yac_str(),
    " - (",
    deriv1_vco2 %>%
      Ryacas::y_fn("Simplify") %>%
      Ryacas::yac_str(),
    ")"
  ) %>%
    parse(text = .) %>%
    Deriv::Deriv(x = "x", nderiv = 1)

  roots <- rootSolve::uniroot.all(
    function(x) deriv1_vo2_func(x) - deriv1_vco2_func(x),
    interval = c(min(.data[[.x]]), max(.data[[.x]]))
  )

  # filter by roots within range of x values
  roots <- roots[roots >= min(.data[[.x]]) &
    roots <= max(.data[[.x]])]

  # filter by which roots have a negative derivative. This indicates that
  # co2 is rising above o2. Finding max finds the last time this occurs.
  final_crossing <- roots[eval(deriv_diff_deriv1,
    envir = list(x = roots)
  ) < 0] %>%
    max()

  # find the fitted VO2 value and use that to find threshold values.
  # You could use VCO2 or somehow integrate both VO2 & VCO2, but that might be a future update

  if (length(final_crossing) > 0) {
    y_hat_threshold <- eval(poly_expr_vo2,
      envir = list(x = final_crossing)
    )
    # get values at threshold
    bp_dat <- find_threshold_vals(
      .data = .data, thr_x = final_crossing,
      thr_y = y_hat_threshold, .x = .x,
      .y = .y, ...
    )

    bp_dat <- bp_dat %>%
      dplyr::mutate(
        bp = bp,
        algorithm = "d1_crossing",
        x_var = .x,
        y_var = .y,
      )
  } else {
    bp_dat <- tibble::tibble()
  }

  if (nrow(bp_dat) == 0) { # no breakpoint found
    bp_dat <- bp_dat %>%
      dplyr::add_row() %>%
      dplyr::mutate(determinant_bp = FALSE)
  } else { # breakpoint found
    bp_dat <- bp_dat %>%
      dplyr::mutate(determinant_bp = TRUE)
  }
  bp_dat <- bp_dat %>%
    dplyr::mutate(
      algorithm = "d1_crossing",
      x_var = .x,
      y_var = .y,
      bp = bp
    ) %>%
    dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp)

  if (plots) {
    bp_plot <- ggplot2::ggplot(
      data = .data,
      ggplot2::aes(
        x = .data[[.x]],
        y = .data[[.y]]
      )
    ) +
      ggplot2::geom_point(ggplot2::aes(y = vo2, color = "vo2"),
        alpha = 0.5
      ) +
      ggplot2::geom_point(ggplot2::aes(y = vco2, color = "vco2"),
        alpha = 0.5
      ) +
      ggplot2::geom_line(ggplot2::aes(y = eval(
        poly_expr_vo2,
        envir = list(x = .data[[.x]])
      ))) +
      ggplot2::geom_line(ggplot2::aes(y = eval(
        poly_expr_vco2,
        envir = list(x = .data[[.x]])
      ))) +
      ggplot2::geom_vline(xintercept = bp_dat[[.x]]) +
      ggplot2::scale_color_manual(values = c(
        "vo2" = "red", "vco2" = "blue"
      )) +
      ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) +
      ggplot2::theme_minimal()
  } else {
    bp_plot <- NULL
  }

  return(list(
    breakpoint_data = bp_dat,
    lm_poly_vo2 = lm_poly_vo2,
    deriv1_vo2 = deriv1_vo2,
    lm_poly_vco2 = lm_poly_vco2,
    deriv1_vco2 = deriv1_vco2,
    bp_plot = bp_plot
  ))
}




## plotting code to figure out later

# max_diff_deriv1 <- max(eval(deriv_diff_deriv1, envir = list(x = df_avg[["time"]])))
#
# sum_range_deriv1 <- eval(deriv_diff_deriv1,
#                          envir = list(x = df_avg[["time"]])) %>%
#     range() %>%
#     abs() %>%
#     sum()
#
# scale_factor <- max(c(.data[[vo2]], .data[[vco2]])) /
#     sum_range_deriv2
#
#
# bp_plot +
#     scale_y_continuous(name = paste(vo2, "&", vco2),
#                        sec.axis = sec_axis(
#                            ~ . / scale_factor*2 - sum_range_deriv2,
#                            name = paste("1st Derivative of",
#                                         vo2, "%", vco2))) +
#     geom_line(ggplot2::aes(y = eval(deriv_diff_deriv1,
#                            envir = list(x = df_avg[["time"]]))),
#               linetype = "dashed")

# this doesn't give the negative values we'll need though
