# d2_reg_spline_maxima(): fit a B-spline regression, take its 2nd derivative, and
# find the highest maxima. Notable implementation choices:
#   1. The 2nd derivative is computed ANALYTICALLY from the B-spline basis
#      (splineDesign(derivs = 2)) rather than by interpolating the fitted values
#      with a second cubic spline. Exact at any degree, no boundary artifact.
#   2. Knot selection mirrors d2_poly_reg_maxima(): a fixed, small number of
#      interior knots by default (n_knots = 2), with an opt-in search (AICc by
#      default, LRT optional).
# The maxima-selection logic lives in select_deriv2_max_idx(), shared by the
# main function and the bootstrap statistic.

#' Analytic derivative of a fitted B-spline regression.
#'
#' Returns a closure with the same `f(x, deriv)` interface as
#' [stats::splinefun()], but evaluates the exact derivative of the fitted
#' regression spline from its B-spline basis instead of re-interpolating the
#' predictions.
#' @keywords internal
#' @noRd
reg_spline_deriv <- function(lm_spline, .x) {
  bs_obj <- stats::model.frame(lm_spline)[[2]] # the bs() basis, carries knots
  ord <- attr(bs_obj, "degree") + 1
  full_knots <- c(
    rep(attr(bs_obj, "Boundary.knots")[1], ord),
    attr(bs_obj, "knots"),
    rep(attr(bs_obj, "Boundary.knots")[2], ord)
  )
  cf <- stats::coef(lm_spline)
  function(x, deriv = 0) {
    basis <- splines::splineDesign(full_knots, x, ord,
      derivs = deriv, outer.ok = TRUE
    )
    # bs() drops the first basis column (absorbed by the intercept), so drop it
    # here too; the derivative of the intercept term is 0.
    vals <- as.vector(basis[, -1] %*% cf[-1])
    if (deriv == 0) vals <- vals + cf[1]
    vals
  }
}

#' Fit a B-spline regression, optionally searching the number of interior knots.
#'
#' @param n_knots Interior knots to use (fixed) or, when searching, the lowest
#'   number considered. B-spline df is `n_knots + degree`.
#' @param search_knots Search for the best number of knots?
#' @param max_knots Upper bound of the knot search.
#' @param criterion `"aicc"` (default) or `"lrt"`.
#' @keywords internal
#' @noRd
loop_d2_reg_spline <- function(.data, .x, .y,
                                          degree = 5,
                                          n_knots = 2,
                                          search_knots = FALSE,
                                          max_knots = 6,
                                          criterion = c("aicc", "lrt"),
                                          alpha_linearity = 0.05) {
  if (nrow(.data) == 0) {
    return(NULL)
  }
  criterion <- rlang::arg_match(criterion)
  if (n_knots < 0) {
    rlang::abort("`n_knots` must be a non-negative integer.")
  }

  fit_spline <- function(k) {
    paste0(
      .y, " ~ 1 + splines::bs(", .x,
      ", df = ", k + degree, ", degree = ", degree, ")"
    ) %>%
      stats::lm(data = .data)
  }

  # default: a single fixed-knot fit (no search)
  if (!search_knots) {
    return(fit_spline(n_knots))
  }

  knots_seq <- seq.int(n_knots, max_knots)

  if (criterion == "aicc") {
    fits <- lapply(knots_seq, fit_spline)
    scores <- vapply(fits, aicc, numeric(1))
    return(fits[[which.min(scores)]])
  }

  # criterion == "lrt": sequential likelihood-ratio test up the knot counts,
  # stopping at the first non-significant improvement (original behavior).
  lm_spline <- fit_spline(knots_seq[1])
  for (k in knots_seq[-1]) {
    lm_next <- fit_spline(k)
    lrt <- stats::anova(lm_spline, lm_next)
    if (is.na(lrt$`Pr(>F)`[2]) || lrt$`Pr(>F)`[2] >= alpha_linearity) {
      break
    }
    lm_spline <- lm_next
  }
  lm_spline
}

#' Locate the chosen 2nd-derivative maximum among candidate sign changes.
#'
#' Shared by the main function and the bootstrap statistic. `spline_func` must
#' accept `(x, deriv)`. Returns the selected index into `equi_spaced_x`
#' (length 0 or 1). This is the original maxima-selection logic, unchanged.
#' @keywords internal
#' @noRd
select_deriv2_max_idx <- function(spline_func, equi_spaced_x, pos_change) {
  if (pos_change) {
    sign_change_idx <- slope_sign_changes(
      y = spline_func(equi_spaced_x, deriv = 2),
      change = "pos_to_neg"
    )
    if (length(sign_change_idx) > 1) {
      neg_slope_idx <- which(
        spline_func(equi_spaced_x[sign_change_idx], deriv = 1) < 0
      )
      while (length(neg_slope_idx) > 0 & length(sign_change_idx) > 1) {
        sign_change_idx <- sign_change_idx[-neg_slope_idx[1]]
        neg_slope_idx <- which(
          spline_func(equi_spaced_x[sign_change_idx], deriv = 1) < 0
        )
      }
    }
  } else {
    sign_change_idx <- slope_sign_changes(
      y = spline_func(equi_spaced_x, deriv = 2),
      change = "neg_to_pos"
    )
    if (length(sign_change_idx) > 1) {
      pos_slope_idx <- which(
        spline_func(equi_spaced_x[sign_change_idx], deriv = 1) > 0
      )
      while (length(pos_slope_idx) > 0 & length(sign_change_idx) > 1) {
        sign_change_idx <- sign_change_idx[-pos_slope_idx[1]]
        pos_slope_idx <- which(
          spline_func(equi_spaced_x[sign_change_idx], deriv = 1) > 0
        )
      }
    }
  }

  # if several remain, take the largest (pos_change) or smallest 2nd-derivative
  if (length(sign_change_idx) > 1) {
    if (pos_change) {
      sign_change_idx <- sign_change_idx[which.max(
        spline_func(equi_spaced_x[sign_change_idx], deriv = 2)
      )]
    } else {
      sign_change_idx <- sign_change_idx[which.min(
        spline_func(equi_spaced_x[sign_change_idx], deriv = 2)
      )]
    }
  }
  sign_change_idx
}

#' @keywords internal
#' @noRd
boot_stat_d2_reg_spline <- function(.data, i, .x, .y, degree,
                                               n_knots, pos_change) {
  lm_spline <- loop_d2_reg_spline(
    .data = .data[i, ], .x = .x, .y = .y,
    degree = degree, n_knots = n_knots, search_knots = FALSE
  )
  if (is.null(lm_spline) || anyNA(stats::coef(lm_spline))) {
    return(NA_real_)
  }
  spline_func <- reg_spline_deriv(lm_spline, .x)
  equi_spaced_x <- seq(
    from = min(.data[i, ][[.x]]), to = max(.data[i, ][[.x]]),
    length.out = round(diff(range(.data[i, ][[.x]])))
  )
  idx <- select_deriv2_max_idx(spline_func, equi_spaced_x, pos_change)
  if (length(idx) == 0) {
    return(NA_real_)
  }
  equi_spaced_x[idx]
}

#' Find a breakpoint using the maxima in the 2nd derivative of a regression spline
#'
#' Fit a B-spline regression to the data, take its second derivative, and find the
#' highest maxima of that second derivative. Traditionally used with ventilatory
#' equivalents (VE/VO2 or VE/VCO2), PetO2, or PetCO2 vs. time or vs. VO2.
#'
#' @details
#' Several works use a similar algorithm, though it is often unclear whether they
#' used "smoothing splines" or "regression splines" (which are not the same). By
#' default this function fits a fixed B-spline regression with \code{n_knots = 2}
#' interior knots. A small, fixed number of knots keeps the fitted curve from
#' developing spurious wiggles that would create extra second-derivative maxima.
#' If you would rather have the number of knots chosen from the data, set
#' \code{search_knots = TRUE}; the search selects a knot count between
#' \code{n_knots} and \code{max_knots} using the corrected AIC (AICc) by default,
#' or a sequential likelihood-ratio test (\code{criterion = "lrt"}).
#'
#' The second derivative is computed analytically from the B-spline basis
#' (\code{splines::splineDesign(derivs = 2)}), which is exact at any degree and
#' avoids the boundary artifacts of interpolating the fitted values with a second
#' spline.
#'
#' @param .data Gas exchange data frame or tibble.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param bp Is this the first (\code{vt1}) or the second (\code{vt2}) breakpoint?
#' @param degree The B-spline (polynomial) degree. Most previous studies use 3 or 5. Default \code{5}.
#' @param n_knots Number of interior knots to use (fixed) or, with \code{search_knots = TRUE}, the lowest number searched. B-spline df is \code{n_knots + degree}. Default \code{2}.
#' @param search_knots Search for the best number of interior knots rather than using the fixed \code{n_knots}? Default \code{FALSE}.
#' @param max_knots Highest number of interior knots considered when \code{search_knots = TRUE}. Default \code{6}.
#' @param criterion Model-selection rule when \code{search_knots = TRUE}: \code{"aicc"} (default) or \code{"lrt"}.
#' @param vo2 Name of the \code{vo2} variable
#' @param vco2 Name of the \code{vco2} variable
#' @param ve Name of the \code{ve} variable
#' @param time Name of the \code{time} variable
#' @param alpha_linearity Significance value used by the \code{"lrt"} knot search.
#' @param pos_change Do you expect the slope to be increasing at the breakpoint? This helps with filtering maxima. Set \code{FALSE} for PetCO2.
#' @param ordering Prior to fitting, should the data be reordered by the x-axis variable or by time? Default reorders by the x variable, using time to break ties.
#' @param front_trim_vt1 Seconds of data to remove from the beginning of the test before fitting (VT1).
#' @param front_trim_vt2 Seconds of data to remove from the beginning of the test before fitting (VT2).
#' @param ci Should the output include confidence interval data? Default \code{FALSE}.
#' @param conf_level Confidence level to use if calculating confidence intervals.
#' @param plots Should this function generate plots? Set to \code{FALSE} to save time.
#' @param ... Passed through so this function plays nicely when \code{breakpoint()} supplies extra arguments.
#'
#' @returns A list whose \code{breakpoint_data} is a slice of the original data
#'   frame at the threshold, with an \code{algorithm} column, plus the fitted
#'   regression spline and its analytic derivative function.
#'
#' @references
#' Leo, J. A., Sabapathy, S., Simmonds, M. J., & Cross, T. J. (2017). The Respiratory Compensation Point is Not a Valid Surrogate for Critical Power. Medicine and Science in Sports and Exercise, 49(7), 1452-1460.
#' Santos, E. L., & Giannella-Neto, A. (2004). Comparison of computerized methods for detecting the ventilatory thresholds. European Journal of Applied Physiology, 93(3), 315-324.
#' Sherrill, D. L., Anderson, S. J., & Swanson, G. (1990). Using smoothing splines for detecting ventilatory thresholds. Medicine and Science in Sports and Exercise, 22(5), 684-689.
#' Wade, T. D., Anderson, S. J., Bondy, J., Ramadevi, V. A., Jones, R. H., & Swanson, G. D. (1988). Using smoothing splines to make inferences about the shape of gas-exchange curves. Computers and Biomedical Research, 21(1), 16-26.
#'
#' @export
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
#' cpet_15b_roll$ve_vco2 <- cpet_15b_roll$ve / cpet_15b_roll$vco2
#'
#' # VT2/RC from the 2nd-derivative maximum of VE/VCO2 vs. time
#' d2_reg_spline_maxima(
#'   cpet_15b_roll,
#'   .x = "time",
#'   .y = "ve_vco2",
#'   bp = "vt2",
#'   vo2 = "vo2",
#'   vco2 = "vco2",
#'   ve = "ve",
#'   time = "time",
#'   plots = FALSE
#' )
d2_reg_spline_maxima <- function(.data,
                                            .x,
                                            .y,
                                            bp,
                                            degree = 5,
                                            n_knots = 2,
                                            search_knots = FALSE,
                                            max_knots = 6,
                                            criterion = c("aicc", "lrt"),
                                            vo2 = "vo2",
                                            vco2 = "vco2",
                                            ve = "ve",
                                            time = "time",
                                            alpha_linearity = 0.05,
                                            pos_change = TRUE,
                                            ordering = c("by_x", "time"),
                                            front_trim_vt1 = 60,
                                            front_trim_vt2 = 60,
                                            ci = FALSE,
                                            conf_level = 0.95,
                                            plots = TRUE,
                                            ...) {
  rlang::check_required(.data)
  rlang::check_required(.x)
  rlang::check_required(.y)
  rlang::check_required(bp)

  bp <- rlang::arg_match(bp, values = c("vt1", "vt2"))
  ordering <- rlang::arg_match(ordering)
  criterion <- rlang::arg_match(criterion)

  .data <- order_cpet_df(.data, .x = .x, time = time, ordering = ordering)
  plot_df <- .data

  front_trim <- set_front_trim(
    bp = bp,
    front_trim_vt1 = front_trim_vt1,
    front_trim_vt2 = front_trim_vt2
  )
  .data <- .data %>%
    dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))

  lm_spline <- loop_d2_reg_spline(
    .data = .data, .x = .x, .y = .y,
    degree = degree, n_knots = n_knots,
    search_knots = search_knots, max_knots = max_knots,
    criterion = criterion, alpha_linearity = alpha_linearity
  )

  if (is.null(lm_spline) || anyNA(stats::coef(lm_spline))) {
    bp_dat <- return_indeterminant_findings(
      .data = plot_df, bp = bp,
      algorithm = "d2_reg_spline_maxima",
      .x = .x, .y = .y, est_ci = "estimate"
    )
    return(list(breakpoint_data = bp_dat))
  }

  # analytic derivative closure (same f(x, deriv) interface as splinefun)
  spline_func <- reg_spline_deriv(lm_spline, .x)

  equi_spaced_x <- seq(
    from = min(.data[[.x]]), to = max(.data[[.x]]),
    length.out = round(diff(range(.data[[.x]])))
  )
  pred <- stats::predict(lm_spline, newdata = tibble::tibble("{.x}" := equi_spaced_x))

  sign_change_idx <- select_deriv2_max_idx(spline_func, equi_spaced_x, pos_change)

  if (length(sign_change_idx) > 0) {
    x_threshold <- equi_spaced_x[sign_change_idx]
    y_hat_threshold <- stats::predict(
      lm_spline, tibble::tibble("{.x}" := x_threshold)
    )
    bp_dat <- find_threshold_vals(
      .data = .data, thr_x = x_threshold, thr_y = y_hat_threshold,
      .x = .x, .y = .y, ...
    )
  } else {
    bp_dat <- tibble::tibble()
  }

  if (nrow(bp_dat) == 0) {
    bp_dat <- bp_dat %>%
      dplyr::add_row() %>%
      dplyr::mutate(determinant_bp = FALSE)
    non_numeric_df <- .data %>%
      dplyr::select(tidyselect::where(
        function(x) is.character(x) | is.factor(x) & all(x == x[1])
      )) %>%
      dplyr::slice(1)
    bp_dat <- dplyr::bind_cols(bp_dat, non_numeric_df)
  } else {
    bp_dat <- bp_dat %>% dplyr::mutate(determinant_bp = TRUE)
  }

  bp_dat <- bp_dat %>%
    dplyr::mutate(
      algorithm = "d2_reg_spline_maxima",
      x_var = .x, y_var = .y, bp = bp, est_ci = "estimate"
    )

  # bootstrapping a CI around an indeterminate estimate is not meaningful, so
  # only run it when a threshold was actually found on the original data.
  if (ci && isTRUE(bp_dat$determinant_bp[1])) {
    poly_df <- length(stats::coef(lm_spline)) - 1 # fitted bs df
    fitted_knots <- poly_df - degree
    boot_res <- bootstrap_threshold_ci(
      .data = .data,
      statistic = boot_stat_d2_reg_spline,
      R = 999, conf_level = conf_level,
      .x = .x, .y = .y, degree = degree,
      n_knots = fitted_knots, pos_change = pos_change
    )

    if (!is.null(boot_res)) {
      ci_lower_x <- boot_res$lower
      ci_lower_y <- stats::predict(lm_spline, tibble::tibble("{.x}" := ci_lower_x))
      lower_ci_res <- find_threshold_vals(
        .data = .data, thr_x = ci_lower_x, thr_y = ci_lower_y,
        .x = .x, .y = .y, ...
      ) %>%
        dplyr::mutate(
          bp = bp, algorithm = "d2_reg_spline_maxima",
          x_var = .x, y_var = .y, est_ci = "lower_ci"
        )

      ci_upper_x <- boot_res$upper
      ci_upper_y <- stats::predict(lm_spline, tibble::tibble("{.x}" := ci_upper_x))
      upper_ci_res <- find_threshold_vals(
        .data = .data, thr_x = ci_upper_x, thr_y = ci_upper_y,
        .x = .x, .y = .y, ...
      ) %>%
        dplyr::mutate(
          bp = bp, algorithm = "d2_reg_spline_maxima",
          x_var = .x, y_var = .y, est_ci = "upper_ci"
        )

      bp_dat <- dplyr::bind_rows(lower_ci_res, bp_dat, upper_ci_res)
      if (any(bp_dat[["determinant_bp"]], na.rm = TRUE)) {
        bp_dat[["determinant_bp"]] <- TRUE
      }
    }
  }

  bp_dat <- bp_dat %>%
    dplyr::relocate(bp, algorithm, x_var, y_var, est_ci, determinant_bp)

  if (plots) {
    bp_plot <- ggplot2::ggplot(
      data = plot_df,
      ggplot2::aes(x = .data[[.x]], y = .data[[.y]])
    ) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_line(
        data = tibble::tibble(x = equi_spaced_x, y = pred),
        ggplot2::aes(x = equi_spaced_x, y = pred)
      ) +
      ggplot2::geom_vline(xintercept = bp_dat[[.x]]) +
      ggplot2::theme_minimal()
  } else {
    bp_plot <- NULL
  }

  list(
    breakpoint_data = bp_dat,
    lm_reg_spline = lm_spline,
    spline_deriv_func = spline_func,
    bp_plot = bp_plot
  )
}
