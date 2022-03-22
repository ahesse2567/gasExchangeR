#' Finding a breakpoint by iteratively moving an indicator variable
#'
#' @param .data Gas exchange data.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param degree The degree of polynomial spline to use. Default is 1 to mimic other methods.
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.
#' @param bp Is this algorithm being used to find vt1 or vt2?
#' @param vo2 The name of the vo2 column in \code{.data}
#' @param vco2 The name of the vco2 column in \code{.data}
#' @param ve The name of the ve column in \code{.data}
#' @param time The name of the time column in \code{.data}
#'
#' @return
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
                      bp) {

    ss <- loop_spline_bp(.data = .data, .x, .y, degree = degree)
    min_ss_idx <- which.min(ss) # re
    .data <- .data %>%
        mutate(s1 = if_else(.data[[.x]] <= .data[[.x]][min_ss_idx], 0,
                            (.data[[.x]] - .data[[.x]][min_ss_idx])^degree))
    lm_spline <- lm(.data[[.y]] ~ 1 + poly(.data[[.x]], degree = degree) + s1,
                data = .data)
    lm_simple <- lm(.data[[.y]] ~ 1 + .data[[.x]], data = .data)

    RSS_simple <- sum(resid(lm_simple)^2)
    RSS_two <- sum(resid(lm_spline)^2)
    MSE_two <- RSS_two / (nrow(.data) - 4) # -4 b/c estimating 4 parameters
    f_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)
    pf_two <- pf(f_stat, df1 = 2, df2 = nrow(.data) - 4, lower.tail = FALSE)
    determinant_bp <- dplyr::if_else(pf_two > alpha_linearity, FALSE, TRUE)

    pred <- bind_rows(x = .data[[.x]],
                      y_hat = lm_spline$fitted.values,
                      method = "spline")

    # slope BEFORE breakpoint = β1
    # slope AFTER breakpoint = β1 + β2
    # this is wrong if degree > 1
    ref <- lm_spline$coefficients[2]
    new <- lm_spline$coefficients[2]+lm_spline$coefficients[3]
    pct_slope_change <- (new - ref)/ref * 100
    bp_dat <- .data[min_ss_idx,] %>%
        select(-s1) %>%
        mutate(bp = bp,
               method = "spline_bp",
               determinant_bp = determinant_bp,
               pct_slope_change = pct_slope_change,
               f_stat = f_stat,
               p_val_f = pf_two) %>%
        relocate(bp, method, determinant_bp, pct_slope_change, f_stat, p_val_f)

    return(list(breakpoint_data = bp_dat,
                fitted_vals = pred,
                lm_spline = lm_spline,
                lm_simple = lm_simple))
}

#' @keywords internal
loop_spline_bp <- function(.data, .x, .y, degree = 1){
    # browser()
    ss_models <- numeric(length = nrow(.data)-2)

    for(i in 2:(nrow(.data)-1)) {
        temp <- .data %>%
            mutate(s1 = if_else(.data[[.x]] <= .data[[.x]][i], 0,
                                (.data[[.x]] - .data[[.x]][i])^degree))
        lm_spline <- lm(temp[[.y]] ~ 1 + poly(temp[[.x]], degree = degree) + s1,
                    data = temp)
        ss_models[i] <- sum((lm_spline$residuals)^2)
    }

    if(any(ss_models == 0)) {
        ss_models[which(ss_models == 0)] <- NA
    }

    return(ss_models)
}
