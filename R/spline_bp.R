#' Finding a breakpoint by iteratively moving an indicator variable
#'
#' This function is similar to other piecewise regression formulas in that it iteratively moves the breakpoint along all n-2 inner data points. We recommend using this function for relationships such as VCO2 vs. VO2 (V-slope), VE vs. VCO2 (respiratory compensation point), and excess CO2, provided that the first ~1 minute of data is removed from excess CO2 prior to fitting the curve.
#'
#' @param .data Gas exchange data.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param degree The degree of polynomial spline to use. Default is 1 to mimic other algorithms.
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.
#' @param bp Is this algorithm being used to find vt1 or vt2?
#' @param vo2 The name of the vo2 column in \code{.data}
#' @param vco2 The name of the vco2 column in \code{.data}
#' @param ve The name of the ve column in \code{.data}
#' @param time The name of the time column in \code{.data}
#' @param pos_change Do you expect the change in slope to be positive (default) or negative? If a two-line regression explains significantly reduces the sum square error but the change in slope does not match the expected underlying physiology, the breakpoint will be classified as indeterminate.
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#' @param front_trim_vt1 How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param front_trim_vt2 How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param pos_slope_after_bp Should the slope after the breakpoint be positive? Default is `TRUE`. This catches cases when the percent change in slope is positive, but the second slope is still negative. Change to `FALSE` when PetCO2 is the y-axis variable.
#'
#' @return A list including slice of the original data frame at the threshold index with new columns `algorithm`, `determinant_bp`, `pct_slope_change`, `f_stat`, and `p_val_f.` The list also includes the fitted values, the left and right sides of the piecewise regression, and a simiple linear regression.
#' @importFrom rlang :=
#' @export
#'
#' @examples
#' # TODO write an example
spline_bp <- function(.data,
                      .x,
                      .y,
                      bp,
                      degree = 1,
                      vo2 = "vo2",
                      vco2 = "vco2",
                      ve = "ve",
                      time = "time",
                      front_trim_vt1 = 60,
                      front_trim_vt2 = 60,
                      alpha_linearity = 0.05,
                      pos_change = TRUE,
                      pos_slope_after_bp = TRUE,
                      ordering = c("by_x", "time"),
                      ...) {
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
    bp <- match.arg(bp, choices = c("vt1", "vt2"), several.ok = FALSE)
    ordering <- match.arg(ordering, several.ok = FALSE)

    .data <- order_cpet_df(.data, .x = .x , time = time,
                           ordering = ordering)
    plot_df <- .data

    front_trim <- set_front_trim(bp = bp,
                                 front_trim_vt1 = front_trim_vt1,
                                 front_trim_vt2 = front_trim_vt2)
    .data <- .data %>%
        dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))

    ss <- loop_spline_bp(.data = .data, .x, .y, degree = degree)
    min_ss_idx <- which.min(ss) # re
    .data <- .data %>%
        dplyr::mutate(s1 = dplyr::if_else(.data[[.x]] <= .data[[.x]][min_ss_idx], 0,
                            (.data[[.x]] - .data[[.x]][min_ss_idx])^degree))
    # create linear model
    lm_spline <- paste0(.y, " ~ ", "1 + poly(", .x, ", degree = ", degree,
                        ", raw = TRUE) + s1") %>%
        stats::lm(data = .data)
    lm_simple <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::lm(data = .data)

    RSS_simple <- sum(stats::resid(lm_simple)^2)
    RSS_two <- sum(stats::resid(lm_spline)^2)
    MSE_two <- RSS_two / (nrow(.data) - 4) # -4 b/c estimating 4 parameters
    f_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)
    pf_two <- stats::pf(f_stat, df1 = 2, df2 = nrow(.data) - 4, lower.tail = FALSE)

    # slope BEFORE breakpoint = β1
    # slope AFTER breakpoint = β1 + β2
    # this is wrong if degree > 1
    slope_left <- lm_spline$coefficients[2]
    slope_right <- lm_spline$coefficients[2]+lm_spline$coefficients[3]
    pct_slope_change <- (slope_right - slope_left)/abs(slope_left) * 100

    determinant_bp <- check_if_determinant_bp(
        p = pf_two,
        pct_slope_change = pct_slope_change,
        pos_change = pos_change,
        pos_slope_after_bp =
            pos_slope_after_bp,
        slope_after_bp = slope_right,
        alpha = alpha_linearity)

    # find last time the indicator variable == 0
    threshold_idx <- which(.data[["s1"]] != 0) %>% min() - 1
    threshold_x <- .data[[.x]][threshold_idx]
    threshold_y <- lm_spline$fitted.values[threshold_idx]

    bp_dat <- find_threshold_vals(.data = .data, thr_x = threshold_x,
                                  thr_y = threshold_y, .x = .x,
                                  .y = .y, ...)
    bp_dat <- bp_dat %>%
        dplyr::mutate(algorithm = "spline_bp",
                      x_var = .x,
                      y_var = .y,
                      determinant_bp = determinant_bp,
                      bp = bp,
                      pct_slope_change = pct_slope_change,
                      f_stat = f_stat,
                      p_val_f = pf_two) %>%
        dplyr::relocate(bp, algorithm, determinant_bp,
                        pct_slope_change, f_stat, p_val_f)

    # add indicator variable to plotting data frame
    plot_df <- plot_df %>%
        dplyr::mutate(s1 = dplyr::if_else(get(.x) <= threshold_x, 0,
                                   (get(.x) - threshold_x)^degree))

    pred <- tibble::tibble("{.x}" := plot_df[[.x]],
                           "{.y}" := stats::predict(lm_spline, newdata = plot_df),
                           algorithm = "spline_bp")

    bp_plot <- ggplot2::ggplot(data = plot_df, aes(x = .data[[.x]], y = .data[[.y]])) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::geom_line(aes(y = pred[[.y]])) +
        ggplot2::geom_vline(xintercept = threshold_x) +
        ggplot2::theme_minimal()

    return(list(breakpoint_data = bp_dat,
                fitted_vals = pred,
                lm_spline = lm_spline,
                lm_simple = lm_simple,
                bp_plot = bp_plot))
}

#' @keywords internal
loop_spline_bp <- function(.data, .x, .y, degree = 1) {
    ss_models <- numeric(length = nrow(.data))

    for(i in 1:nrow(.data)) {
        if(i == 1 | i == nrow(.data)) {
            ss_models[i] <- NA
            next
        }
        temp <- .data %>% # Generate data with indicator variable (knot)
            dplyr::mutate(s1 = dplyr::if_else(.data[[.x]] <= .data[[.x]][i], 0,
                                (.data[[.x]] - .data[[.x]][i])^degree))
        # create linear model
        lm_spline <- paste0(.y, " ~ ", "1 + poly(", .x, ", degree = ", degree,
                            ", raw = TRUE) + s1") %>%
            stats::lm(data = temp)
        # record residual sums of squares for later comparison
        ss_models[i] <- sum((lm_spline$residuals)^2)
    }

    if(any(ss_models == 0, is.na(ss_models))) {
        ss_models[which(ss_models == 0)] <- NA
    }

    return(ss_models)
}
