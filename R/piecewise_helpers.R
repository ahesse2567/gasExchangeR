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

    plot_x <- .data[[.x]]
    plot_df <- tibble("{.x}" := plot_x) %>%
        mutate(y_hat_left = predict(lm_left, tibble("{.x}" := plot_x)),
               y_hat_right = predict(lm_right, tibble("{.x}" := plot_x)))

    bp_plot <- ggplot2::ggplot(data = .data,
                    aes(x = .data[[.x]], y = .data[[.y]])) +
        geom_point(alpha = 0.5) +
        geom_line(data = plot_df, aes(x = get(.x), y = y_hat_left)) +
        geom_line(data = plot_df, aes(x = get(.x), y = y_hat_right)) +
        geom_vline(xintercept = bp_dat[[.x]]) +
        theme_minimal()
    bp_plot
}
