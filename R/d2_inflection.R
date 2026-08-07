#' Second derivative inflection
#'
#' Find a breakpoint in ventilatory data based on the 2nd derivative inflection point after Santos et al. (2004). In their paper, they use pulmonary ventilation (VE) vs. VO2 to find the first ventilatory (VT1) / gas exchange threshold (GET) and VE vs. VCO2 to find the respiratory compensation point (RCP) / second ventilatory threshold (VT2). It is generally suitable to also use VE vs. time to find VT2. One can truncate the test and rerun VE vs. time to find VT1.
#'
#' @param .data Gas exchange data frame or tibble.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param bp Is this the first (\code{vt1}) or the second (\code{vt2}) breakpoint?
#' @param degree The polynomial degree fit to the curve. Default \code{5}, matching Santos and Giannella-Neto (2004). When \code{search_degree = TRUE}, this is the lowest degree searched. Must be at least 5.
#' @param search_degree Search for the best polynomial degree rather than using the fixed \code{degree}? Default \code{FALSE}.
#' @param max_degree Highest polynomial degree considered when \code{search_degree = TRUE}. Default \code{8}.
#' @param criterion Model-selection rule when \code{search_degree = TRUE}: \code{"aicc"} (default) or \code{"lrt"}.
#' @param vo2 Name of the \code{vo2} variable
#' @param vco2 Name of the \code{vco2} variable
#' @param ve Name of the \code{ve} variable
#' @param time Name of the \code{time} variable
#' @param alpha_linearity Significance value used by the \code{"lrt"} degree search.
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param front_trim_vt1 Seconds of data to remove from the beginning of the test before fitting (VT1). Default \code{60}.
#' @param front_trim_vt2 Seconds of data to remove from the beginning of the test before fitting (VT2). Default \code{60}.
#' @param ci Should the output include confidence interval data? Default is `FALSE`. When `TRUE`, a nonparametric (case-resampling) bootstrap with 999 replicates is used to estimate a percentile confidence interval, refitting at the same fixed \code{degree} used for the estimate.
#' @param conf_level Confidence level to use if calculating confidence intervals.
#' @param plots Should this function generate plots? Set to `FALSE` to save time.
#'
#' @returns A list containing breakpoint data, best-fit and related functions, and a plot.
#' @export
#'
#' @details To mimic the Santos et al. (2004) paper, this fits a 5th-order polynomial regression to the VE vs. VO2 or VCO2 data by default. Set \code{search_degree = TRUE} to instead select the degree from the data (between \code{degree} and \code{max_degree}) using the corrected AIC (\code{"aicc"}, default) or a sequential likelihood-ratio test (\code{"lrt"}). After establishing the polynomial, we find the real roots of the 2nd derivative within the given range of the `.x` variable to find the inflection points. We then filter by concavity changes that go from up to down because this matches the expected shape of the pulmonary ventilation curves.
#'
#' We are unsure how to handle multiple inflection points at this time, but this is unlikely to be an issue if the polynomial degree is kept equal to 5.
#'
#' @references
#' Santos, E. L., & Giannella-Neto, A. (2004). Comparison of computerized methods for detecting the ventilatory thresholds. European journal of applied physiology, 93, 315-324.
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
#' # RCP from the inflection point of VE vs. VCO2
#' d2_inflection(
#'   cpet_15b_roll,
#'   .x = "vco2",
#'   .y = "ve",
#'   bp = "vt2",
#'   degree = 5,
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
d2_inflection <- function(.data,
                          .x,
                          .y,
                          bp,
                          ...,
                          degree = 5,
                          search_degree = FALSE,
                          max_degree = 8,
                          criterion = c("aicc", "lrt"),
                          vo2 = "vo2",
                          vco2 = "vco2",
                          ve = "ve",
                          time = "time",
                          alpha_linearity = 0.05,
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

  .data <- order_cpet_df(.data, .x = .x, time = time, ordering = ordering)

  # remove the noisy start of the test before fitting (as the other threshold
  # functions do), to avoid spurious early inflection points.
  front_trim <- set_front_trim(
    bp = bp,
    front_trim_vt1 = front_trim_vt1,
    front_trim_vt2 = front_trim_vt2
  )
  .data <- .data %>%
    dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))

  # fit the polynomial and locate the chosen inflection point (shared with the
  # bootstrap statistic so this logic exists in exactly one place)
  infl_res <- find_d2_inflection(
    .data, .x = .x, .y = .y,
    degree = degree, search_degree = search_degree,
    max_degree = max_degree, criterion = criterion,
    alpha_linearity = alpha_linearity
  )
  lm_poly <- infl_res$lm_poly

  # return quick summary if generating the model fails
  if (is.null(lm_poly) || anyNA(stats::coef(lm_poly))) {
    bp_dat <- return_indeterminant_findings(
      .data = .data, bp = bp,
      algorithm = "d2_inflection",
      .x = .x, .y = .y, est_ci = "estimate"
    )
    return(list(breakpoint_data = bp_dat))
  }

  poly_expr <- infl_res$poly_expr
  deriv1 <- infl_res$deriv1
  deriv2 <- infl_res$deriv2
  inflection_points_x <- infl_res$inflection_x

  if (!is.na(inflection_points_x)) {
    inflection_points_y <- eval(poly_expr, envir = list(x = inflection_points_x))
    bp_dat <- find_threshold_vals(
      .data = .data, thr_x = inflection_points_x,
      thr_y = inflection_points_y, .x = .x,
      .y = .y, ...
    )
  } else {
    bp_dat <- tibble::tibble()
  }

  if (nrow(bp_dat) == 0) { # no breakpoint found
    bp_dat <- bp_dat %>%
      dplyr::add_row() %>%
      dplyr::mutate(determinant_bp = FALSE)
    # add character or factor columns that are all the same value (e.g. ids)
    non_numeric_df <- .data %>%
      dplyr::select(tidyselect::where(function(x) {
        is.character(x) | is.factor(x) &
          all(x == x[1])
      })) %>%
      dplyr::slice(1)
    bp_dat <- dplyr::bind_cols(bp_dat, non_numeric_df)
  } else { # breakpoint found
    bp_dat <- bp_dat %>%
      dplyr::mutate(determinant_bp = TRUE)
  }
  bp_dat <- bp_dat %>%
    dplyr::mutate(
      bp = bp,
      algorithm = "d2_inflection",
      x_var = .x,
      y_var = .y,
      est_ci = "estimate"
    ) %>%
    dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp)

  # bootstrapping a CI around an indeterminate estimate is not meaningful, so
  # only run it when a threshold was actually found on the original data.
  if (ci && isTRUE(bp_dat$determinant_bp[1])) {
    boot_res <- bootstrap_threshold_ci(
      .data = .data,
      statistic = boot_stat_d2_inflection,
      R = 999,
      conf_level = conf_level,
      .x = .x, .y = .y,
      degree = degree, alpha_linearity = alpha_linearity
    )

    if (!is.null(boot_res)) {
      y_lower <- eval(poly_expr, envir = list(x = boot_res$lower))
      lower_ci_res <- find_threshold_vals(
        .data = .data, thr_x = boot_res$lower, thr_y = y_lower,
        .x = .x, .y = .y, ...
      ) %>%
        dplyr::mutate(
          bp = bp, algorithm = "d2_inflection",
          x_var = .x, y_var = .y, est_ci = "lower_ci"
        )

      y_upper <- eval(poly_expr, envir = list(x = boot_res$upper))
      upper_ci_res <- find_threshold_vals(
        .data = .data, thr_x = boot_res$upper, thr_y = y_upper,
        .x = .x, .y = .y, ...
      ) %>%
        dplyr::mutate(
          bp = bp, algorithm = "d2_inflection",
          x_var = .x, y_var = .y, est_ci = "upper_ci"
        )

      bp_dat <- dplyr::bind_rows(lower_ci_res, bp_dat, upper_ci_res) %>%
        dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp)
      if (any(bp_dat[["determinant_bp"]], na.rm = TRUE)) {
        bp_dat[["determinant_bp"]] <- TRUE
      }
    }
  }

  # make plot for output
  if (plots) {
    bp_plot <- ggplot2::ggplot(
      data = .data,
      ggplot2::aes(x = .data[[.x]], y = .data[[.y]])
    ) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_line(ggplot2::aes(y = eval(poly_expr,
        envir = list(x = .data[[.x]])
      ))) +
      ggplot2::geom_vline(xintercept = bp_dat[[.x]]) +
      ggplot2::theme_minimal()
  } else {
    bp_plot <- NULL
  }

  return(list(
    breakpoint_data = bp_dat,
    lm_poly_reg = lm_poly,
    deriv1_expr = deriv1,
    deriv2_expr = deriv2,
    bp_plot = bp_plot
  ))
}

#' Fit the polynomial and locate the chosen "up to down" inflection point.
#'
#' Shared by [d2_inflection()] and its bootstrap statistic
#' ([boot_stat_d2_inflection()]) so this logic exists in exactly one place.
#' Returns `inflection_x = NA_real_` when the model fails to fit or no
#' qualifying inflection point exists.
#'
#' Root-finding uses exact polynomial arithmetic (see `poly_deriv_coefs()`,
#' `poly_eval_coefs()` in utils.R) rather than symbolic differentiation
#' (Deriv::Deriv) plus a Ryacas-simplify/string-parse/polyroot() route. This
#' is both faster (Ryacas, an external computer algebra system, is the
#' dominant cost when this function is called ~1000 times for a bootstrap CI)
#' and simpler: the 2nd derivative's coefficients come directly from the
#' original coefficients via two closed-form differentiation steps, with no
#' string round-trip through a symbolic simplifier. The symbolic expressions
#' are still built once (cheaply, via Deriv::Deriv with no Ryacas call)
#' purely so the return value keeps its existing shape for
#' plotting/inspection.
#'
#' @keywords internal
#' @noRd
find_d2_inflection <- function(.data, .x, .y, degree,
                               search_degree = FALSE, max_degree = 8,
                               criterion = c("aicc", "lrt"),
                               alpha_linearity = 0.05,
                               build_exprs = TRUE) {
  lm_poly <- fit_poly_reg(
    .data = .data, .x = .x, .y = .y,
    degree = degree, search_degree = search_degree,
    max_degree = max_degree, criterion = criterion,
    alpha_linearity = alpha_linearity
  )

  if (is.null(lm_poly) || anyNA(stats::coef(lm_poly))) {
    return(list(
      inflection_x = NA_real_, lm_poly = lm_poly,
      poly_expr = NULL, deriv1 = NULL, deriv2 = NULL
    ))
  }

  coefs <- stats::coef(lm_poly)
  deriv1_coefs <- poly_deriv_coefs(coefs)
  deriv2_coefs <- poly_deriv_coefs(deriv1_coefs)

  # roots of the 2nd derivative (an exact polynomial) are the inflection
  # candidates; polyroot() finds them directly, no symbolic engine needed.
  inflection_points_x <- find_real(polyroot(deriv2_coefs))

  # filter by roots within range of x values
  rng <- range(.data[[.x]])
  inflection_points_x <- inflection_points_x[
    inflection_points_x >= rng[1] & inflection_points_x <= rng[2]
  ]

  # keep only "concave up to concave down" points: the 2nd derivative goes
  # from + to - across the point (tol matches the previous default).
  tol <- 0.1
  left_sign <- sign(poly_eval_coefs(deriv2_coefs, inflection_points_x - tol))
  right_sign <- sign(poly_eval_coefs(deriv2_coefs, inflection_points_x + tol))
  inflection_points_x <- inflection_points_x[left_sign > 0 & right_sign < 0]

  # if there's more than one up to down inflection point, choose whichever
  # point has the higher slope. which.max() on an empty vector correctly
  # returns integer(0), so no valid inflection point yields NA below.
  inflection_points_x <- inflection_points_x[which.max(
    poly_eval_coefs(deriv1_coefs, inflection_points_x)
  )]

  if (length(inflection_points_x) == 0) {
    inflection_points_x <- NA_real_
  }

  # build symbolic expressions once (cheap; no Ryacas) purely for the return
  # value's existing shape -- not needed for the root-finding above, and
  # skipped by the bootstrap statistic, where it would run ~1000x per call.
  if (build_exprs) {
    poly_expr <- expr_from_coefs(coefs)
    deriv1 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 1)
    deriv2 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 2)
  } else {
    poly_expr <- deriv1 <- deriv2 <- NULL
  }

  list(
    inflection_x = inflection_points_x,
    lm_poly = lm_poly, poly_expr = poly_expr,
    deriv1 = deriv1, deriv2 = deriv2
  )
}

#' Bootstrap statistic for [d2_inflection()]'s confidence interval.
#'
#' Refits the polynomial on the resampled data at the fixed `degree` selected
#' from the full data (no re-searching per replicate) and returns the
#' inflection x-location, or `NA_real_` if this resample has none.
#' @keywords internal
#' @noRd
boot_stat_d2_inflection <- function(.data, i, .x, .y, degree, alpha_linearity) {
  find_d2_inflection(
    .data[i, ],
    .x = .x, .y = .y,
    degree = degree, search_degree = FALSE,
    alpha_linearity = alpha_linearity,
    build_exprs = FALSE # skip building symbolic expressions per replicate
  )$inflection_x
}

