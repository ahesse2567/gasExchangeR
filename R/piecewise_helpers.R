#' Statistics of Piecewise Regressions Compared to a Simple Linear Regression
#'
#' @keywords internal
piecewise_stats <- function(lm_left, lm_right, lm_simple) {
    RSS_simple <- sum(stats::resid(lm_simple)^2)
    RSS_two <- sum(stats::resid(lm_left)^2) + sum(stats::resid(lm_right)^2)
    MSE_two <- RSS_two / (nrow(lm_simple$model) - 4) # -4 b/c estimating 4 parameters
    f_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)
    pf_two <- stats::pf(f_stat, df1 = 2, df2 = nrow(lm_simple$model) - 4,
                        lower.tail = FALSE)

    list(RSS_simple = RSS_simple, RSS_two = RSS_two, MSE_two = MSE_two,
         f_stat = f_stat, pf_two = pf_two)
}

#' Intersection point between two lines
#'
#' This function takes two \code{lm} objects in the form of \code{y = mx + b} and returns a vector of their \code{(x, y)} intersection point.
#'
#' @keywords internal
#' @noRd
intersection_point <- function(lm1, lm2) {
    stopifnot(class(lm1) == "lm",
              class(lm2) == "lm",
              length(lm1$coefficients) == 2,
              length(lm2$coefficients) == 2)
    # y = mx + b
    # mx1 + b1 = mx2 + b2 at the intersection point
    # x = (b2 - b1) / (m1 - m2)
    x <- (lm2$coefficients[1] - lm1$coefficients[1]) /
        (lm1$coefficients[2] - lm2$coefficients[2])
    y <- lm1$coefficients[1] + lm1$coefficients[2] * x

    point <- c(x, y)
    names(point) <- c('x', 'y')
    point
}

#' @keywords internal
#' @noMd
make_piecewise_bp_plot <- function(.data, .x, .y, lm_left, lm_right,
                                   bp_dat, ...) {
    # right now this function adds upper and lower CI values if supplied
    # I should probably make that optional later
    plot_x <- .data[[.x]]
    plot_df <- tibble::tibble("{.x}" := plot_x) %>%
        dplyr::mutate(y_hat_left = stats::predict(lm_left,
                                           tibble::tibble("{.x}" := plot_x)),
               y_hat_right = stats::predict(lm_right,
                                            tibble::tibble("{.x}" := plot_x)))

    bp_plot <- ggplot2::ggplot(data = .data,
                               ggplot2::aes(x = .data[[.x]], y = .data[[.y]])) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::geom_line(data = plot_df, ggplot2::aes(x = get(.x), y = y_hat_left)) +
        ggplot2::geom_line(data = plot_df, ggplot2::aes(x = get(.x), y = y_hat_right)) +
        ggplot2::geom_vline(xintercept = bp_dat[[.x]]) +
        ggplot2::theme_minimal()
    bp_plot
}

#' Check if the breakpoint is a "true" breakpoint
#'
#' There must be a better fit than a single regression line p < alpha the slope must change in the direction indicated (usually positive) the direction of the second slope must match the expected direction  (we usually expect the second slope to be positive). This is necessary because you can have a positive percent change but still have a negative second slope.
#' @keywords internal
#' @noRd
#'
check_if_determinant_bp <- function(.data,
                                    .x,
                                    p,
                                    pct_slope_change,
                                    pos_change,
                                    pos_slope_after_bp,
                                    slope_after_bp,
                                    range_x,
                                    alpha = 0.05) {

    determinant_bp <- dplyr::if_else(
        all(p < alpha,
            pos_change == (pct_slope_change > 0),
            pos_slope_after_bp == (slope_after_bp > 0)),
        TRUE,
        FALSE)
    determinant_bp
}

#' @keywords internal
#' @noMd
get_best_piecewise_idx <- function(loop_res_df,
                                   data_range,
                                   alpha_linearity,
                                   pos_change,
                                   pos_slope_after_bp) {
    # filter loop results for the best-fit, determinant solution
    best_idx <- loop_res_df %>%
        dplyr::filter(p < alpha_linearity &
                          pos_change == pos_change &
                          pos_slope_after_bp == pos_slope_after_bp &
                          dplyr::between(int_point_x,
                                         data_range[1],
                                         data_range[2]),
                      inside_ci) %>%
        dplyr::filter(p == min(p)) %>% # filter by smallest p-value (lowest RSS)
        dplyr::select(idx) %>%
        dplyr::pull()

    if(length(best_idx) > 1) {
        best_idx <- sample(best_idx, 1)
    } else if(length(best_idx) == 0) {
        # no valid solution, but produce next-best result anyway
        best_idx <- which.min(loop_res_df$p)
    }

    best_idx
}
