#' Finding a breakpoint using the Jones-Molitoris algorithm.
#'
#' @param .data Gas exchange data.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.
#' @param bp Is this algorithm being used to find vt1 or vt2?
#' @param vo2 The name of the vo2 column in \code{.data}
#' @param vco2 The name of the vco2 column in \code{.data}
#' @param ve The name of the ve column in \code{.data}
#' @param time The name of the time column in \code{.data}
#' @param pos_change Do you expect the change in slope to be positive (default) or negative? If a two-line regression explains significantly reduces the sum square error but the change in slope does not match the expected underlying physiology, the breakpoint will be classified as indeterminate.
#' @param front_trim_vt1 How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param front_trim_vt2 How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param pos_slope_after_bp Should the slope after the breakpoint be positive? Default is `TRUE`. This catches cases when the percent change in slope is positive, but the second slope is still negative. Change to `FALSE` when PetCO2 is the y-axis variable.
#' @param ci Should the output include confidence interval data? Default is `FALSE`.
#' @param conf_level Confidence level to use if calculating confidence intervals.
#' @param plots Should this function generate plots? Set to `FALSE` to save time.
#'
#' @returns A list including slice of the original data frame at the threshold index with new columns `algorithm`, `determinant_bp`, `pct_slope_change`, `f_stat`, and `p_val_f.` The list also includes the fitted values, the left and right sides of the piecewise regression, and a simple linear regression.
#'
#' @importFrom rlang :=
#'
#' @export
#'
#' @references
#' Jones, R. H., & Molitoris, B. A. (1984). A statistical method for determining the breakpoint of two lines. Analytical Biochemistry, 141(1), 287–290. https://doi.org/10.1016/0003-2697(84)90458-5
#'
#' @examples
#'
#' # Load breath-by-breath graded exercise testing file
#' cpet_bbb <- utils::read.csv(
#'   system.file("extdata", "anton_vo2max_clean.csv", package = "gasExchangeR")
#' )
#'
#' # 15-breath (rolling breath) average
#' cpet_15b_roll <- avg_exercise_test(
#'   cpet_bbb,
#'   method = "breath",
#'   calc_type = "rolling",
#'   time_col = "time",
#'   window = 15,
#'   align = "center",
#'   mos = "mean"
#' )
#'
#' # Find a specific breakpoint using the JM algorithm. Note: using breakpoint
#' # algorithms directly does NOT perform best practices such as removing the
#' # first 60 seconds of exercise.
#' jm_res <- jm(.data = cpet_15b_roll,
#'   .x = "vco2",
#'   .y = "ve",
#'   bp = "vt2"
#' )
#' jm_res
#'

jm <- function(.data,
               .x,
               .y,
               bp,
               ...,
               vo2 = "vo2",
               vco2 = "vco2",
               ve = "ve",
               time = "time",
               ordering = c("by_x", "time"),
               alpha_linearity = 0.05,
               front_trim_vt1 = 60,
               front_trim_vt2 = 60,
               pos_change = TRUE,
               pos_slope_after_bp = TRUE,
               ci = FALSE,
               conf_level = 0.95,
               plots = TRUE) {

  stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
  bp <- match.arg(bp, choices = c("vt1", "vt2"), several.ok = FALSE)

  .data <- order_cpet_df(.data,
                         .x = .x,
                         time = time,
                         ordering = ordering)
  plot_df <- .data

  # determime how much data to remove
  front_trim <- set_front_trim(
    bp = bp,
    front_trim_vt1 = front_trim_vt1,
    front_trim_vt2 = front_trim_vt2
  )

  # remove data from the start of the test
  .data <- .data %>%
    dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))
  loop_res <- loop_jm(
    .data = .data,
    .x = .x,
    .y = .y,
    alpha_linearity = alpha_linearity,
    conf_level = conf_level
  )

  if (!is.null(loop_res)) {
    best_idx <- get_best_piecewise_idx(
      loop_res,
      range(.data[[.x]]),
      alpha_linearity = alpha_linearity,
      pos_change = pos_change,
      pos_slope_after_bp = pos_slope_after_bp
    )
  } else {
    best_idx <- numeric()
  }

  # return quick summary if generating models fails
  if (is.null(loop_res) | length(best_idx) == 0) {
    bp_dat <- return_indeterminant_findings(
      .data = plot_df,
      bp = bp,
      algorithm = as.character(match.call()[[1]]),
      .x = .x,
      .y = .y,
      est_ci = "estimate"
    )

    return(list(breakpoint_data = bp_dat))
  }

  estimate_res <- get_jm_res(
    .data = .data,
    bp_idx = best_idx,
    .x = .x,
    .y = .y,
    bp = bp,
    alpha_linearity = alpha_linearity,
    pos_change = pos_change,
    pos_slope_after_bp = pos_slope_after_bp,
    est_ci = "estimate"
  )
  if (ci) {
    # in the future, create a function to make DYRer code
    ci_lower_idx <- loop_res %>%
      dplyr::filter(inside_ci) %>%
      dplyr::filter(int_point_x == min(int_point_x)) %>%
      dplyr::filter(p == min(p)) %>% # for breaking ties
      dplyr::select(idx) %>%
      dplyr::pull()

    lower_ci_res <- get_jm_res(
      .data = .data,
      bp_idx = ci_lower_idx,
      .x = .x,
      .y = .y,
      bp = bp,
      alpha_linearity = alpha_linearity,
      pos_change = pos_change,
      pos_slope_after_bp = pos_slope_after_bp,
      est_ci = "lower"
    )

    ci_upper_idx <- loop_res %>%
      dplyr::filter(inside_ci) %>%
      dplyr::filter(int_point_x == max(int_point_x)) %>%
      dplyr::filter(p == min(p)) %>% # for breaking ties
      dplyr::select(idx) %>%
      dplyr::pull()

    upper_ci_res <- get_jm_res(
      .data = .data,
      bp_idx = ci_upper_idx,
      .x = .x,
      .y = .y,
      bp = bp,
      alpha_linearity = alpha_linearity,
      pos_change = pos_change,
      pos_slope_after_bp = pos_slope_after_bp,
      est_ci = "upper"
    )

    # combine estimate and both CI breakpoint res into one tibble
    estimate_res$bp_dat <- dplyr::bind_rows(
      lower_ci_res$bp_dat,
      estimate_res$bp_dat,
      upper_ci_res$bp_dat
    )
  }

  if (plots) {
    plot_x <- plot_df[[.x]]
    plot_df <- tibble::tibble("{.x}" := plot_x) %>%
      dplyr::mutate(
        y_hat_left = stats::predict(
          estimate_res$lm_left,
          tibble::tibble("{.x}" := plot_x)
        ),
        y_hat_right = stats::predict(
          estimate_res$lm_right,
          tibble::tibble("{.x}" := plot_x)
        ) +
          estimate_res$b0_plus_b1x0
      )

    bp_plot <- ggplot2::ggplot(
      data = .data,
      ggplot2::aes(x = .data[[.x]], y = .data[[.y]])
    ) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_line(
        data = plot_df,
        ggplot2::aes(x = get(.x), y = y_hat_left)
      ) +
      ggplot2::geom_line(
        data = plot_df,
        ggplot2::aes(x = get(.x), y = y_hat_right)
      ) +
      ggplot2::geom_vline(xintercept = estimate_res$bp_dat[[.x]]) +
      ggplot2::theme_minimal()
  } else {
    bp_plot <- NULL
  }

  return(list(
    breakpoint_data = estimate_res$bp_dat,
    fitted_vals = estimate_res$pred,
    lm_left = estimate_res$lm_left,
    lm_right = estimate_res$lm_right,
    lm_simple = estimate_res$lm_simple,
    bp_plot = bp_plot
  ))
}

#' @keywords internal
#' Scores every division of the data into a left OLS segment and a right segment
#' constrained to join it at the breakpoint (the Jones-Molitoris model). The
#' O(n) search is delegated to piecewise_metrics_cumsum_jm() (cumulative
#' sufficient statistics), which is mathematically identical to fitting a fresh
#' left OLS and a through-join right regression at every split, but far faster.
loop_jm <- function(.data,
                    .x,
                    .y,
                    alpha_linearity = 0.05,
                    conf_level = 0.95) {
  loop_stats <- piecewise_metrics_cumsum_jm(
    .data = .data, .x = .x, .y = .y, conf_level = conf_level
  )
  if (is.null(loop_stats)) {
    return(NULL)
  }

  loop_stats %>%
    dplyr::select(
      idx, p, pos_change, pos_slope_after_bp, int_point_x, inside_ci
    )
}

#' @keywords internal
#' Similar to loop_jm(), but gets results at specific index
get_jm_res <- function(.data, bp_idx, .x, .y, bp,
                       alpha_linearity, pos_change,
                       pos_slope_after_bp,
                       est_ci = c("estimate", "lower_ci", "upper_ci"),
                       ...) {
  est_ci <- match.arg(est_ci, several.ok = FALSE)

  df_left <- .data[1:bp_idx, ] # x < x0
  df_right <- .data[bp_idx:nrow(.data), ] # x >= x0

  x0 <- .data[[.x]][bp_idx]

  lm_left <- paste0(.y, " ~ ", "1 + ", .x) %>%
    stats::lm(data = df_left)

  # get y value we will force lm_right to pass through
  b0_plus_b1x0 <- lm_left$coefficients[1] +
    lm_left$coefficients[2] * x0

  # use I() function to force model through (x0, b0_plus_b1x0)
  lm_right <- paste0(
    "I(", .y, "-", b0_plus_b1x0, ")",
    " ~ 0 + I(", .x, " - ", x0, ")"
  ) %>%
    stats::lm(data = df_right)

  lm_simple <- paste0(.y, " ~ ", "1 + ", .x) %>%
    stats::lm(data = .data)

  # check for a significant departure from linearity
  pw_stats <- piecewise_stats(lm_left, lm_right, lm_simple)
  # list2env(pw_stats, envir = environment())

  y_hat_left <- tibble::tibble("{.x}" := df_left[[.x]],
                               "{.y}" := lm_left$fitted.values,
                               algorithm = "jm"
  )
  y_hat_right <- tibble::tibble("{.x}" := df_right[[.x]],
                                "{.y}" := lm_right$fitted.values + b0_plus_b1x0,
                                algorithm = "jm"
  )
  pred <- dplyr::bind_rows(y_hat_left, y_hat_right)
  # need to use lm_right$coefficients[1] instead of [2] b/c there is only
  # a slope given we forced the line through the origin to satisfy the algorithm
  pct_slope_change <- 100 * (lm_right$coefficients[1] - lm_left$coefficients[2]) /
    abs(lm_left$coefficients[2])

  determinant_bp <- check_if_determinant_bp(
    p = pw_stats$pf_two,
    pct_slope_change = pct_slope_change,
    pos_change = pos_change,
    pos_slope_after_bp =
      pos_slope_after_bp,
    slope_after_bp = stats::coef(lm_right)[1],
    alpha = alpha_linearity
  )

  bp_dat <- find_threshold_vals(
    .data = .data, thr_x = x0,
    thr_y = b0_plus_b1x0, .x = .x,
    .y = .y, ...
  )
  bp_dat <- bp_dat %>%
    dplyr::mutate(
      algorithm = "jm",
      est_ci = est_ci,
      x_var = .x,
      y_var = .y,
      determinant_bp = determinant_bp,
      bp = bp,
      pct_slope_change = pct_slope_change,
      f_stat = pw_stats$ f_stat,
      p_val_f = pw_stats$pf_two
    ) %>%
    dplyr::relocate(
      bp, algorithm, determinant_bp, est_ci,
      pct_slope_change, f_stat, p_val_f
    )

  list(
    bp_dat = bp_dat,
    fitted_vals = pred,
    lm_left = lm_left,
    lm_right = lm_right,
    lm_simple = lm_simple,
    b0_plus_b1x0 = b0_plus_b1x0
  )
}
