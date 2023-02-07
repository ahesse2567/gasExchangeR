#' Finding a breakpoint using Orr's 'bruteforce' algorithm.
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
#' @param front_trim How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#'
#' @returns A list including slice of the original data frame at the threshold index with new columns `algorithm`, `determinant_bp`, `pct_slope_change`, `f_stat`, and `p_val_f.` The list also includes the fitted values, the left and right sides of the piecewise regression, and a simple linear regression.
#'
#' @importFrom rlang :=
#' @export
#'
#' @references
#' Orr, G. W., Green, H. J., Hughson, R. L., & Bennett, G. W. (1982). A computer linear regression model to determine ventilatory anaerobic threshold. Journal of Applied Physiology Respiratory Environmental and Exercise Physiology, 52(5), 1349–1352. https://doi.org/10.1152/jappl.1982.52.5.1349
#'
#' @examples
#' # TODO write examples
orr <- function(.data,
                .x,
                .y,
                bp,
                vo2 = "vo2",
                vco2 = "vco2",
                ve = "ve",
                time = "time",
                alpha_linearity = 0.05,
                front_trim = 60,
                pos_change = TRUE,
                ...) {
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]]) %>%
        filter(.data[[time]] >= min(.data[[time]] + front_trim))

    ss <- loop_orr(.data = .data, .x = .x, .y = .y)
    min_ss_idx <- which.min(ss)

    df_left <- .data[1:min_ss_idx,]
    df_right <- .data[min_ss_idx:nrow(.data),]

    lm_left <- stats::lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)
    lm_right <- stats::lm(df_right[[.y]] ~ 1 + df_right[[.x]], data = df_right)
    lm_simple <- stats::lm(.data[[.y]] ~ 1 + .data[[.x]], data = .data)

    RSS_simple <- sum(stats::resid(lm_simple)^2)
    RSS_two <- sum(stats::resid(lm_left)^2) + sum(stats::resid(lm_right)^2)
    MSE_two <- RSS_two / (nrow(.data) - 4) # -4 b/c estimating 4 parameters
    f_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)
    pf_two <- stats::pf(f_stat, df1 = 2, df2 = nrow(.data) - 4, lower.tail = FALSE)

    pct_slope_change <- 100*(lm_right$coefficients[2] - lm_left$coefficients[2]) /
        abs(lm_left$coefficients[2])

    y_hat_left <- tibble::tibble("{.x}" := df_left[[.x]],
                         "{.y}" := lm_left$fitted.values,
                         algorithm = "orr")
    y_hat_right <- tibble::tibble("{.x}" := df_right[[.x]],
                          "{.y}" := lm_right$fitted.values,
                          algorithm = "orr")
    pred <- dplyr::bind_rows(y_hat_left, y_hat_right)

    determinant_bp <- dplyr::if_else(pf_two < alpha_linearity &
                                         (pos_change == (pct_slope_change > 0)),
                                     TRUE, FALSE)

    int_point <- intersection_point(lm_left, lm_right)

    bp_dat <- .data %>%
        dplyr::mutate(dist_x_sq = (.data[[.x]] - int_point["x"])^2,
               dist_y_sq = (.data[[.y]] - int_point["y"])^2,
               sum_dist_sq = dist_x_sq + dist_y_sq) %>%
        dplyr::arrange(sum_dist_sq) %>%
        dplyr::slice(1) %>%
        dplyr::mutate(algorithm = "orr",
                      x_var = .x,
                      y_var = .y,
                      determinant_bp = determinant_bp,
                      bp = bp) %>%
        dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp) %>%
        dplyr::select(-c(dist_x_sq, dist_y_sq, sum_dist_sq))

    bp_dat <- bp_dat %>%
        dplyr::mutate(pct_slope_change = pct_slope_change,
               f_stat = f_stat,
               p_val_f = pf_two) %>%
        dplyr::relocate(bp, algorithm, determinant_bp, pct_slope_change, f_stat, p_val_f)

    return(list(breakpoint_data = bp_dat,
                fitted_vals = pred,
                lm_left = lm_left,
                lm_right = lm_right,
                lm_simple = lm_simple))

    # Should we add the three line regression code later?

}
#'
#' @keywords internal
loop_orr <- function(.data, .x, .y) {
    # browser()
    ss_both <- vector(length = nrow(.data))

    for(i in 1:nrow(.data)) {
        if(i == 1 | i == nrow(.data)) {
            ss_both[i] <- NA
            next
        }
        # split data into left and right halves
        # should these share the same point? Or should, they be different by 1 point?
        # Ekkekakis said WB makes them differ by 1 pt
        # I wonder exactly what the word 'contiguous' means in this context
        df_left <- .data[1:i,]
        df_right <- .data[i:nrow(.data),]

        # make linear models of the two halves
        lm_left <- stats::lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)
        lm_right <- stats::lm(df_right[[.y]] ~ 1 + df_right[[.x]], data = df_right)

        ss_left <- sum((lm_left$residuals)^2)
        ss_right <- sum((lm_right$residuals)^2)
        #ss_both <- c(ss_both, ss_left + ss_right)
        ss_both[i] <- (ss_left + ss_right)
    }
    return(ss_both)
}
