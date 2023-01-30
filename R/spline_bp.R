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
                      degree = 1,
                      vo2 = "vo2",
                      vco2 = "vco2",
                      ve = "ve",
                      time = "time",
                      alpha_linearity = 0.05,
                      bp,
                      pos_change = TRUE) {
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])

    ss <- loop_spline_bp(.data = .data, .x, .y, degree = degree)
    min_ss_idx <- which.min(ss) # re
    .data <- .data %>%
        dplyr::mutate(s1 = dplyr::if_else(.data[[.x]] <= .data[[.x]][min_ss_idx], 0,
                            (.data[[.x]] - .data[[.x]][min_ss_idx])^degree))
    lm_spline <- stats::lm(.data[[.y]] ~ 1 + poly(.data[[.x]], degree = degree) + s1,
                data = .data)
    lm_simple <- stats::lm(.data[[.y]] ~ 1 + .data[[.x]], data = .data)

    RSS_simple <- sum(stats::resid(lm_simple)^2)
    RSS_two <- sum(stats::resid(lm_spline)^2)
    MSE_two <- RSS_two / (nrow(.data) - 4) # -4 b/c estimating 4 parameters
    f_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)
    pf_two <- stats::pf(f_stat, df1 = 2, df2 = nrow(.data) - 4, lower.tail = FALSE)

    # slope BEFORE breakpoint = β1
    # slope AFTER breakpoint = β1 + β2
    # this is wrong if degree > 1
    ref <- lm_spline$coefficients[2]
    new <- lm_spline$coefficients[2]+lm_spline$coefficients[3]
    pct_slope_change <- (new - ref)/ref * 100

    determinant_bp <- dplyr::if_else(pf_two < alpha_linearity &
                                         (pos_change == (pct_slope_change > 0)),
                                     TRUE, FALSE)

    bp_dat <- .data[min_ss_idx,] %>%
        dplyr::select(-s1) %>%
        dplyr::mutate(bp = bp,
               algorithm = "spline_bp",
               x_var = .x,
               y_var = .y,
               determinant_bp = determinant_bp,
               pct_slope_change = pct_slope_change,
               f_stat = f_stat,
               p_val_f = pf_two) %>%
        dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp,
                 pct_slope_change, f_stat, p_val_f)

    pred <- tibble::tibble("{.x}" := .data[[.x]],
                   "{.y}" := lm_spline$fitted.values,
                   algorithm = "spline_bp")

    return(list(breakpoint_data = bp_dat,
                fitted_vals = pred,
                lm_spline = lm_spline,
                lm_simple = lm_simple))
}

#' @keywords internal
loop_spline_bp <- function(.data, .x, .y, degree = 1) {
    # browser()
    ss_models <- numeric(length = nrow(.data))

    for(i in 1:nrow(.data)) {
        if(i == 1 | i == nrow(.data)) {
            ss_models[i] <- NA
            next
        }
        temp <- .data %>%
            dplyr::mutate(s1 = dplyr::if_else(.data[[.x]] <= .data[[.x]][i], 0,
                                (.data[[.x]] - .data[[.x]][i])^degree))
        lm_spline <- stats::lm(temp[[.y]] ~ 1 + poly(temp[[.x]], degree = degree) + s1,
                    data = temp)
        ss_models[i] <- sum((lm_spline$residuals)^2)
    }

    if(any(ss_models == 0)) {
        ss_models[which(ss_models == 0)] <- NA
    }

    return(ss_models)
}
