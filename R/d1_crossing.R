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
#'   degree = 6,
#'   search_degree = FALSE,
#'   max_degree = 8,
#'   criterion = "aicc",
#'   vo2 = "vo2",
#'   vco2 = "vco2",
#'   ve = "ve",
#'   time = "time",
#'   ci = FALSE,
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

  # find VO2/VCO2 polynomial fits and the final derivative crossing (shared
  # with the bootstrap statistic so this logic exists in exactly one place)
  crossing_res <- find_d1_crossing(
    .data, time = time, vo2 = vo2, vco2 = vco2,
    degree = degree, search_degree = search_degree,
    max_degree = max_degree, criterion = criterion,
    alpha_linearity = alpha_linearity
  )
  lm_poly_vo2 <- crossing_res$lm_poly_vo2
  lm_poly_vco2 <- crossing_res$lm_poly_vco2

  # if creating the linear model fails, return a quick summary
  if (is.null(lm_poly_vo2) ||
    is.null(lm_poly_vco2) ||
    anyNA(lm_poly_vo2$coefficients) ||
    anyNA(lm_poly_vco2$coefficients)) {
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

  poly_expr_vo2 <- crossing_res$poly_expr_vo2
  poly_expr_vco2 <- crossing_res$poly_expr_vco2
  deriv1_vo2 <- crossing_res$deriv1_vo2
  deriv1_vco2 <- crossing_res$deriv1_vco2
  final_crossing <- crossing_res$crossing_x

  # find the fitted VO2 value and use that to find threshold values.
  # You could use VCO2 or somehow integrate both VO2 & VCO2, but that might be a future update

  if (!is.na(final_crossing)) {
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
        est_ci = "estimate"
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
      bp = bp,
      est_ci = "estimate"
    ) %>%
    dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp)

  # bootstrapping a CI around an indeterminate estimate is not meaningful, so
  # only run it when a threshold was actually found on the original data.
  if (ci && isTRUE(bp_dat$determinant_bp[1])) {
    boot_res <- bootstrap_threshold_ci(
      .data = .data,
      statistic = boot_stat_d1_crossing,
      R = 999,
      conf_level = conf_level,
      time = time, vo2 = vo2, vco2 = vco2,
      degree = degree, alpha_linearity = alpha_linearity
    )

    if (!is.null(boot_res)) {
      y_lower <- eval(poly_expr_vo2, envir = list(x = boot_res$lower))
      lower_ci_res <- find_threshold_vals(
        .data = .data, thr_x = boot_res$lower, thr_y = y_lower,
        .x = .x, .y = .y, ...
      ) %>%
        dplyr::mutate(
          bp = bp, algorithm = "d1_crossing",
          x_var = .x, y_var = .y, est_ci = "lower_ci"
        )

      y_upper <- eval(poly_expr_vo2, envir = list(x = boot_res$upper))
      upper_ci_res <- find_threshold_vals(
        .data = .data, thr_x = boot_res$upper, thr_y = y_upper,
        .x = .x, .y = .y, ...
      ) %>%
        dplyr::mutate(
          bp = bp, algorithm = "d1_crossing",
          x_var = .x, y_var = .y, est_ci = "upper_ci"
        )

      bp_dat <- dplyr::bind_rows(lower_ci_res, bp_dat, upper_ci_res) %>%
        dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp)
      if (any(bp_dat[["determinant_bp"]], na.rm = TRUE)) {
        bp_dat[["determinant_bp"]] <- TRUE
      }
    }
  }

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

#' Fit VO2 and VCO2 vs. time and locate the final 1st-derivative crossing.
#'
#' Shared by [d1_crossing()] and its bootstrap statistic ([boot_stat_d1_crossing()])
#' so this logic exists in exactly one place. Returns `crossing_x = NA_real_`
#' (not `-Inf`) when the models fail to fit or no valid crossing exists; a
#' previous version of this logic took `max()` of a possibly-empty vector,
#' which silently returns `-Inf` rather than signaling "not found."
#'
#' Root-finding uses exact polynomial arithmetic (see `poly_deriv_coefs()`,
#' `poly_eval_coefs()` in utils.R) rather than symbolic differentiation
#' (Deriv::Deriv) plus Ryacas simplification and rootSolve::uniroot.all().
#' This is both faster (the symbolic + Ryacas route is the dominant cost when
#' this function is called ~1000 times for a bootstrap CI) and, since the
#' derivative of a polynomial is itself an exact polynomial, MORE reliable:
#' polyroot() finds every root of that polynomial exactly (to numerical
#' precision), whereas uniroot.all() is a numerical scan that can miss closely
#' spaced roots depending on its sampling grid. The symbolic expressions are
#' still built once (cheaply, via Deriv::Deriv with no Ryacas call) purely so
#' the return value keeps its existing shape for plotting/inspection.
#'
#' @keywords internal
#' @noRd
find_d1_crossing <- function(.data, time, vo2, vco2, degree,
                             search_degree = FALSE, max_degree = 8,
                             criterion = c("aicc", "lrt"),
                             alpha_linearity = 0.05,
                             build_exprs = TRUE) {
  lm_poly_vo2 <- fit_poly_reg(
    .data = .data, .x = time, .y = vo2,
    degree = degree, search_degree = search_degree,
    max_degree = max_degree, criterion = criterion,
    alpha_linearity = alpha_linearity
  )
  lm_poly_vco2 <- fit_poly_reg(
    .data = .data, .x = time, .y = vco2,
    degree = degree, search_degree = search_degree,
    max_degree = max_degree, criterion = criterion,
    alpha_linearity = alpha_linearity
  )

  failed <- is.null(lm_poly_vo2) || is.null(lm_poly_vco2) ||
    anyNA(lm_poly_vo2$coefficients) || anyNA(lm_poly_vco2$coefficients)
  if (failed) {
    return(list(
      crossing_x = NA_real_,
      lm_poly_vo2 = lm_poly_vo2, lm_poly_vco2 = lm_poly_vco2,
      poly_expr_vo2 = NULL, poly_expr_vco2 = NULL,
      deriv1_vo2 = NULL, deriv1_vco2 = NULL
    ))
  }

  coefs_vo2 <- stats::coef(lm_poly_vo2)
  coefs_vco2 <- stats::coef(lm_poly_vco2)
  deriv1_vo2_coefs <- poly_deriv_coefs(coefs_vo2)
  deriv1_vco2_coefs <- poly_deriv_coefs(coefs_vco2)

  # the difference of the two 1st derivatives is itself an exact polynomial;
  # its roots ARE the crossing candidates.
  diff_coefs <- poly_subtract_coefs(deriv1_vo2_coefs, deriv1_vco2_coefs)
  roots <- find_real(polyroot(diff_coefs))

  rng <- range(.data[[time]])
  roots <- roots[roots >= rng[1] & roots <= rng[2]]

  # keep roots where CO2 is rising above O2 (negative derivative-of-difference,
  # i.e. the difference polynomial is decreasing there); the last such
  # crossing is the threshold. Check length BEFORE max() -- max(numeric(0))
  # is -Inf, not an empty/NA result.
  diff_deriv_coefs <- poly_deriv_coefs(diff_coefs)
  neg_roots <- roots[poly_eval_coefs(diff_deriv_coefs, roots) < 0]
  crossing_x <- if (length(neg_roots) > 0) max(neg_roots) else NA_real_

  # build symbolic expressions once (cheap; no Ryacas) purely for the return
  # value's existing shape -- not needed for the root-finding above, and
  # skipped by the bootstrap statistic, where it would run ~1000x per call.
  if (build_exprs) {
    poly_expr_vo2 <- expr_from_coefs(coefs_vo2)
    deriv1_vo2 <- Deriv::Deriv(poly_expr_vo2, x = "x", nderiv = 1)
    poly_expr_vco2 <- expr_from_coefs(coefs_vco2)
    deriv1_vco2 <- Deriv::Deriv(poly_expr_vco2, x = "x", nderiv = 1)
  } else {
    poly_expr_vo2 <- poly_expr_vco2 <- deriv1_vo2 <- deriv1_vco2 <- NULL
  }

  list(
    crossing_x = crossing_x,
    lm_poly_vo2 = lm_poly_vo2, lm_poly_vco2 = lm_poly_vco2,
    poly_expr_vo2 = poly_expr_vo2, poly_expr_vco2 = poly_expr_vco2,
    deriv1_vo2 = deriv1_vo2, deriv1_vco2 = deriv1_vco2
  )
}

#' Bootstrap statistic for [d1_crossing()]'s confidence interval.
#'
#' Refits both polynomials on the resampled data at the fixed `degree`
#' selected from the full data (no re-searching per replicate) and returns
#' the crossing x-location, or `NA_real_` if this resample has no valid
#' crossing.
#' @keywords internal
#' @noRd
boot_stat_d1_crossing <- function(.data, i, time, vo2, vco2, degree,
                                  alpha_linearity) {
  find_d1_crossing(
    .data[i, ],
    time = time, vo2 = vo2, vco2 = vco2,
    degree = degree, search_degree = FALSE,
    alpha_linearity = alpha_linearity,
    build_exprs = FALSE # skip building symbolic expressions per replicate
  )$crossing_x
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
