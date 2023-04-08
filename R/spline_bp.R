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
#' @param front_trim_vt1 How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param front_trim_vt2 How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param pos_slope_after_bp Should the slope after the breakpoint be positive? Default is `TRUE`. This catches cases when the percent change in slope is positive, but the second slope is still negative. Change to `FALSE` when PetCO2 is the y-axis variable.
#' @param ci Should the output include confidence interval data? Default is `FALSE`.
#' @param conf_level Confidence level to use if calculating confidence intervals.
#' @param plots Should this function generate plots? Set to `FALSE` to save time.
#' @param ... Dot dot dot mostly allows this function to work properly if this function gets passed arguments that are not strictly necessary.
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
                      ...,
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
                      ci = FALSE,
                      conf_level = 0.95,
                      plots = TRUE
                      ) {
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

    loop_res <- loop_spline_bp(.data = .data, .x, .y,
                         alpha_linearity = alpha_linearity,
                         conf_level = conf_level)

    best_idx <- get_best_piecewise_idx(loop_res,
                                       range(.data[[.x]]),
                                       alpha_linearity = alpha_linearity,
                                       pos_change = pos_change,
                                       pos_slope_after_bp = pos_slope_after_bp)

    estimate_res <- get_spline_bp_res(.data = .data,
                                      bp_idx = best_idx,
                                      .x = .x,
                                      .y = .y,
                                      bp = bp,
                                      alpha_linearity = alpha_linearity,
                                      pos_change = pos_change,
                                      pos_slope_after_bp = pos_slope_after_bp,
                                      est_ci = "estimate")
    # create indicator variable value
    threshold_x <- estimate_res$bp_dat[[.x]]
    # add indicator variable to plotting data frame
    plot_df <- plot_df %>%
        dplyr::mutate(s1 = dplyr::if_else(get(.x) <= threshold_x,
                                          0,
                                          get(.x) - threshold_x))

    pred <- tibble::tibble("{.x}" := plot_df[[.x]],
                           "{.y}" := stats::predict(estimate_res$lm_spline,
                                                    newdata = plot_df),
                           algorithm = "spline_bp")

    if(ci) {
        ci_lower_idx <- loop_res %>%
            dplyr::filter(inside_ci) %>%
            dplyr::filter(int_point_x == min(int_point_x)) %>%
            dplyr::select(idx) %>%
            dplyr::pull()

        lower_ci_res <- get_spline_bp_res(.data = .data,
                                          bp_idx = ci_lower_idx,
                                          .x = .x,
                                          .y = .y,
                                          bp = bp,
                                          alpha_linearity = alpha_linearity,
                                          pos_change = pos_change,
                                          pos_slope_after_bp = pos_slope_after_bp,
                                          est_ci = "lower")

        ci_upper_idx <- loop_res %>%
            dplyr::filter(inside_ci) %>%
            dplyr::filter(int_point_x == max(int_point_x)) %>%
            dplyr::select(idx) %>%
            dplyr::pull()

        upper_ci_res <- get_spline_bp_res(.data = .data,
                                          bp_idx = ci_upper_idx,
                                          .x = .x,
                                          .y = .y,
                                          bp = bp,
                                          alpha_linearity = alpha_linearity,
                                          pos_change = pos_change,
                                          pos_slope_after_bp = pos_slope_after_bp,
                                          est_ci = "upper")

        # combine estimate and both CI breakpoint res into one tibble
        estimate_res$bp_dat <- rbind(lower_ci_res$bp_dat,
                                     estimate_res$bp_dat,
                                     upper_ci_res$bp_dat)
    }

    if(plots){

        threshold_x <- estimate_res$bp_dat[[.x]]

        bp_plot <- ggplot2::ggplot(data = plot_df,
                                   ggplot2::aes(x = .data[[.x]], y = .data[[.y]])) +
            ggplot2::geom_point(alpha = 0.5) +
            ggplot2::geom_line(aes(y = pred[[.y]])) +
            ggplot2::geom_vline(xintercept = threshold_x) +
            ggplot2::theme_minimal()
    } else {
        bp_plot <- NULL
    }

    return(list(breakpoint_data = estimate_res$bp_dat,
                fitted_vals = pred,
                lm_spline = estimate_res$lm_spline,
                lm_simple = estimate_res$lm_simple,
                bp_plot = bp_plot))
}

#' @keywords internal
loop_spline_bp <- function(.data, .x, .y, alpha_linearity = 0.05,
                           conf_level = 0.95) {
    # perhaps later let people specify a degree, but for now, it'll be hard
    # breakpoint
    n_rows <- nrow(.data)

    # initialize empty vectors
    ss_models <- MSE_two <- f_stat <- pf_two <-
        pct_slope_change <- int_point_x <- numeric(n_rows)
    pos_change <- pos_slope_after_bp <- logical(n_rows)

    # create simple linear model and calculate its RSS
    lm_simple <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::lm(data = .data)
    RSS_simple <- sum(stats::resid(lm_simple)^2)

    for(i in 1:nrow(.data)) {
        if(i == 1 | i == nrow(.data)) {
            ss_models[i] <- MSE_two[i] <- f_stat[i] <-
                pf_two[i] <- pos_change[i] <- pos_slope_after_bp[i] <-
                int_point_x[i] <- NA
            next
        }

        int_point_x[i] <- .data[[.x]][i]
        temp <- .data %>% # Generate data with indicator variable (knot)
            dplyr::mutate(s1 = dplyr::if_else(.data[[.x]] <= .data[[.x]][i], 0,
                                .data[[.x]] - .data[[.x]][i]))
        # (.data[[.x]] - .data[[.x]][i])^degree
        # create linear model
        lm_spline <- paste0(.y, " ~ ", "1 + poly(", .x, # add back degree here
                            ", raw = TRUE) + s1") %>%
            stats::lm(data = temp)

        slope_left <- lm_spline$coefficients[2]
        slope_right <- lm_spline$coefficients[2]+lm_spline$coefficients[3]
        pct_slope_change[i] <- (slope_right - slope_left)/abs(slope_left) * 100
        pos_change[i] <- dplyr::if_else(pct_slope_change[i] > 0, TRUE, FALSE)
        pos_slope_after_bp[i] <- dplyr::if_else(slope_right > 0, TRUE, FALSE)

        # record residual sums of squares for later comparison
        ss_models[i] <- sum((lm_spline$residuals)^2)
        MSE_two[i] <- ss_models[i] / (nrow(lm_simple$model) - 4) # -4 b/c estimating 4 parameters
        f_stat[i] <- (RSS_simple - ss_models[i]) / (2 * MSE_two[i])
        pf_two[i] <- stats::pf(f_stat[i], df1 = 2,
                               df2 = nrow(lm_simple$model) - 4,
                               lower.tail = FALSE)
    }

    crit_F <- stats::qf(conf_level, 1, n_rows - 4, lower.tail = TRUE)
    inside_ci <- dplyr::if_else(
        ((ss_models - min(ss_models, na.rm = TRUE)) / MSE_two) < crit_F,
        TRUE, FALSE)
    # for debugging purposes, plot breakpoints inside 95% CI
    # plot((ss_models - min(ss_models, na.rm = TRUE)) / MSE_two)
    # abline(h = crit_F)

    loop_stats <- tibble::tibble(p = pf_two,
                                 pos_change = pos_change,
                                 pos_slope_after_bp = pos_slope_after_bp,
                                 int_point_x = int_point_x,
                                 inside_ci = inside_ci
    ) %>%
        dplyr::mutate(idx = dplyr::row_number())

    loop_stats
}

get_spline_bp_res <- function(.data, bp_idx, .x, .y, bp,
                              alpha_linearity, pos_change,
                              pos_slope_after_bp,
                              est_ci = c("estimate", "lower_ci", "upper_ci"),
                              ...) {
    est_ci <- match.arg(est_ci, several.ok = FALSE)

    .data <- .data %>% # Generate data with indicator variable (knot)
        dplyr::mutate(s1 = dplyr::if_else(.data[[.x]] <= .data[[.x]][bp_idx],
                                          0,
                                          .data[[.x]] - .data[[.x]][bp_idx]))
    # (.data[[.x]] - .data[[.x]][i])^degree
    # create linear model

    lm_spline <- paste0(.y, " ~ ", "1 + ", .x, # add back degree here
                        " + s1") %>%
        stats::lm(data = .data)
    lm_simple <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::lm(data = .data)

    n_rows <- nrow(.data)

    RSS_simple <- sum(stats::resid(lm_simple)^2)
    RSS_two <- sum(stats::resid(lm_spline)^2)
    MSE_two <- RSS_two / (n_rows - 4) # -4 b/c estimating 4 parameters
    f_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)
    pf_two <- stats::pf(f_stat, df1 = 2, df2 = n_rows - 4, lower.tail = FALSE)

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
                      est_ci = est_ci,
                      x_var = .x,
                      y_var = .y,
                      determinant_bp = determinant_bp,
                      bp = bp,
                      pct_slope_change = pct_slope_change,
                      f_stat = f_stat,
                      p_val_f = pf_two) %>%
        dplyr::relocate(bp, algorithm, determinant_bp, est_ci,
                        pct_slope_change, f_stat, p_val_f)

    list(bp_dat = bp_dat,
         lm_spline = lm_spline,
         lm_simple = lm_simple)
}
