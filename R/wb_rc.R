#' Find the Respiratory Compensation Point as Implemented by WinBreak 3.7
#'
#' This is an adaptation of WinBreak 3.7's (Epistemic Mindworks, 2003) method to find the respiratory compensation point (RCP), also known as the second ventilatory threshold (VT2). Their method is an adaptation of the original V-slope method by Beaver et al. (1986).
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
#' @param min_pct_change The original V-slope method (Beaver et al., 1986) considers the respiratory compensation point as determinant "if the change in slope"...between the left and right regressions..."is greater than a preselected amount (15\% of the initial slope)." See details for the specific WinBreak implementation.
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#'
#' @returns A list including slice of the original data frame at the threshold index with new columns `algorithm`, `determinant_bp`, `pct_slope_change`, `f_stat`, and `p_val_f.` The list also includes the fitted values, the left and right sides of the piecewise regression, and a simple linear regression.
#'
#' @details This uses a method similar to Ekkekakis et al. (2004) that adjusted the Beaver et al. (1986) method. In the original V-slope method, the RCP is the *first* time is the considered the first time the percent change from the left to the right regression exceeds 15\% (Beaver et al., 1986). However, Ekkekakis et al. (2004) also required a significant departure from linearity *and* used the largest percent increase between slopes observed by all divisions of the data. The reason, they state, is "because using a fixed slope difference (such as the 15\% mentioned by Beaver et al.) occasionally resulted in untenable results (e.g. points below 50\% VO2max)."
#'
#' We note that this improvement can still yield innacurate results based on sudden changes in the tails of the data. For example the right regression may only contain two data points (the minimum to form a line) with a rather steep and positive slope because the very last data point is higher than the penultimate data point. The percent change from the left to the right regression at the second-to-last data point may therefore be the largest, and the function will choose this data point, despite the fact that the best-fit solution that still had a percent increase in slope of > 15\% was at a much earlier division of the data. Therefore, our implementation of WinBreak 3.7's (Epistemic Mindworks, 2003) RC method finds the division of the data that yields the smallest residual sums of squares, is a significant departure from linearity (significant F-test against a single regression line), and has a percent increase in slope by at least the stated percent (default = 15\%).
#'
#' @references
#' Beaver, W. L., Wasserman, K., & Whipp, B. J. (1986). A new method for detecting anaerobic threshold by gas exchange. Journal of Applied Physiology, 121(6), 2020–2027.
#' Ekkekakis, P., Lind, E., Hall, E. E., & Petruzzello, S. J. (2008). Do regression-based computer algorithms for determining the ventilatory threshold agree?. Journal of Sports Sciences, 26(9), 967-976.
#' Epistemic Mindworks. (2003). WinBreak User Guide. Epistemic Mindworks.
#'
#' @export
#'
#' @examples
#'
#' # TODO write an example
#'
wb_rc <- function(.data,
                  .x,
                  .y,
                  bp,
                  vo2 = "vo2",
                  vco2 = "vco2",
                  ve = "ve",
                  time = "time",
                  alpha_linearity = 0.05,
                  front_trim_vt1 = 60,
                  front_trim_vt2 = 60,
                  pos_change = TRUE,
                  min_pct_change = 0.15,
                  ...) {

    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
    bp = match.arg(bp, choices = c("vt1", "vt2"), several.ok = FALSE)
    if(.x != vco2 | .y != ve) {
        warning("The `x` variable may NOT be VCO2 or y variable may NOT be VE. This method is designed to work with x = VO2 and y = VE")
    }

    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])
    plot_df <- .data # save full data frame for plotting later

    front_trim <- set_front_trim(bp = bp,
                                 front_trim_vt1 = front_trim_vt1,
                                 front_trim_vt2 = front_trim_vt2)

    .data <- .data %>%
        dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))

    loop_idx <- loop_wb_rc(.data = .data, .x = .x, .y = .y,
                           min_pct_change = min_pct_change)

    lm_simple <- stats::lm(paste0(.y, " ~ ", .x), data = .data)

    if(is.na(loop_idx)) {
        bp_dat <- .data[0,] %>%
            dplyr::mutate(algorithm = "wb_rc",
                          x_var = .x,
                          y_var = .y,
                          determinant_bp = FALSE,
                          bp = bp,
                          pct_slope_change = NA,
                          f_stat = NA,
                          p_val_f = NA) %>%
            dplyr::relocate(bp, algorithm, determinant_bp,
                            pct_slope_change, f_stat, p_val_f)
        return(list(breakpoint_data = bp_dat,
                    fitted_vals = pred,
                    lm_left = NA,
                    lm_right = NA,
                    lm_simple = lm_simple,
                    bp_plot = NA))
    } else {
        df_left <- .data[1:loop_idx,]
        df_right <- .data[loop_idx:nrow(.data),]

        lm_left <- stats::lm(paste0(.y, " ~ ", .x), data = df_left)
        lm_right <- stats::lm(paste0(.y, " ~ ", .x), data = df_right)

        pw_stats <- piecewise_stats(lm_left, lm_right, lm_simple)
        list2env(pw_stats, envir = environment())

        pct_slope_change <- 100*(lm_right$coefficients[2] - lm_left$coefficients[2]) /
            abs(lm_left$coefficients[2])

        y_hat_left <- tibble::tibble("{.x}" := df_left[[.x]],
                                     "{.y}" := lm_left$fitted.values,
                                     algorithm = "wb_rc")
        y_hat_right <- tibble::tibble("{.x}" := df_right[[.x]],
                                      "{.y}" := lm_right$fitted.values,
                                      algorithm = "wb_rc")
        pred <- dplyr::bind_rows(y_hat_left, y_hat_right)

        determinant_bp <- dplyr::if_else(pf_two < alpha_linearity &
                                             (pos_change == (pct_slope_change > 0)),
                                         TRUE, FALSE)

        int_point <- intersection_point(lm_left, lm_right)

        bp_dat <- find_threshold_vals(.data = .data, thr_x = int_point["x"],
                                      int_point["y"], .x = .x,
                                      .y = .y, ...)
        bp_dat <- bp_dat %>%
            dplyr::mutate(algorithm = "wb_rc",
                          x_var = .x,
                          y_var = .y,
                          determinant_bp = determinant_bp,
                          bp = bp,
                          pct_slope_change = pct_slope_change,
                          f_stat = f_stat,
                          p_val_f = pf_two) %>%
            dplyr::relocate(bp, algorithm, determinant_bp,
                            pct_slope_change, f_stat, p_val_f)

        bp_plot <- make_piecewise_bp_plot(plot_df, .x, .y, lm_left, lm_right, bp_dat)

        return(list(breakpoint_data = bp_dat,
                    fitted_vals = pred,
                    lm_left = lm_left,
                    lm_right = lm_right,
                    lm_simple = lm_simple,
                    bp_plot = bp_plot))
    }
}

#' @keywords internal
loop_wb_rc <- function(.data, .x, .y, min_pct_change = 0.15) {

    # fit a simple linear model with all the data for a comparison
    lm_simple <- lm(paste0(.y, " ~ ", .x), data = .data)
    RSS_simple <- sum(stats::resid(lm_simple)^2)
    df2 <- nrow(lm_simple$model) - 4 # -4 b/c estimating 4 parameters

    # initialize results storage
    f_stats <- numeric(length = nrow(.data))
    p_values <- numeric(length = nrow(.data))
    pct_changes <- numeric(length = nrow(.data))

    # if at the first or last data points, set values to NA. Doing this so
    # the length of the vectors match the number of rows in the data frame
    # this way the index matches up
    f_stats[c(1, length(f_stats))] <- NA
    p_values[c(1, length(p_values))] <- NA
    pct_changes[c(1, length(pct_changes))] <- NA

    # loop through divisions in data and calculate metrics
    for (i in 2:(nrow(.data) - 1)) {
        # split data into left and right portions
        df_left <- .data[1:i, ]
        df_right <- .data[i:nrow(.data), ]

        # fit piecewise linear model to left and right portions using the .x and .y variables
        lm_left <- stats::lm(paste0(.y, " ~ ", .x), data = df_left)
        lm_right <- stats::lm(paste0(.y, " ~ ", .x), data = df_right)

        # calculate and store the pooled sums of squares in the left and right portions
        ss_left <- sum(stats::resid(lm_left)^2)
        ss_right <- sum(stats::resid(lm_right)^2)
        RSS_two <- sum(c(ss_left, ss_right))
        MSE_two <- RSS_two / df2
        # compare the sums of squares in the simple linear model against
        # the pooled sums of squares in the piecewise model using an F test
        f_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)
        p_value <- stats::pf(f_stat, df1 = 2, df2 = df2,
                             lower.tail = FALSE)
        # calculate and store the percent change from the left to the right regression
        pct_change <- 100 * (stats::coefficients(lm_right)[2] - stats::coefficients(lm_left)[2]) /
            abs(stats::coefficients(lm_left)[2])

        # store values
        f_stats[i] <- f_stat
        p_values[i] <- p_value
        pct_changes[i] <- pct_change
    }

    # return the index of the division of the data where there was a significant departure from linearity
    # and also a percent change equal to or greater than min_pct_change
    if (any(p_values < 0.05) & any(pct_changes >= min_pct_change)) {
        return(which(p_values == min(p_values, na.rm = TRUE) & pct_changes >= min_pct_change))
    } else {
        return(NA)
    }
}
