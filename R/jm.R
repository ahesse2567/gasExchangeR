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
#'
#' @return
#' @export
#'
#' @references
#' Jones, R. H., & Molitoris, B. A. (1984). A statistical method for determining the breakpoint of two lines. Analytical Biochemistry, 141(1), 287–290. https://doi.org/10.1016/0003-2697(84)90458-5
#'
#' @examples
#' # TODO write examples
jm <- function(.data,
               .x,
               .y,
               vo2 = "vo2",
               vco2 = "vco2",
               ve = "ve",
               time = "time",
               alpha_linearity = 0.05,
               bp,
               pos_change = TRUE) {
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
    bp <- match.arg(bp, choices = c("vt1", "vt2"), several.ok = FALSE)
    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])

    ss <- loop_jm(.data = .data, .x = .x, .y = .y)
    min_ss_idx <- which.min(ss)

    # I don't think this actually intersects at x0 right now!!!!

    df_left <- .data[1:min_ss_idx,] # x < x0
    df_right <- .data[(min_ss_idx):nrow(.data),] # x >= x0

    x_knot <- .data[[.x]][min_ss_idx]

    df_right <- df_right %>%
        mutate(s1 = df_right[[.x]] - x_knot)

    lm_left <- lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)
    # according to the JM algorithm, the right regression line will have a constant equal to b0 + b1*x0
    b0_plus_b1x0 <- lm_left$coefficients[1] + lm_left$coefficients[2] * .data[[.x]][min_ss_idx]

    # for some reason, using I() function to try to force the regression to use b0_plus_b1x0 as an intercept doesn't work, while the offset argument does. The slope coefficient (b3) is the same, but the predicted values are not.
    lm_right <- lm(df_right[[.y]] ~ 0 + s1, data = df_right,
                   offset = rep(b0_plus_b1x0, nrow(df_right)))

    lm_simple <- lm(.data[[.y]] ~ 1 + .data[[.x]])

    RSS_simple <- sum(resid(lm_simple)^2)
    RSS_two <- sum(resid(lm_left)^2) + sum(resid(lm_right)^2)
    MSE_two <- RSS_two / (nrow(.data) - 4) # -4 b/c estimating 4 parameters
    f_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)
    pf_two <- pf(f_stat, df1 = 2, df2 = nrow(.data) - 4, lower.tail = FALSE)

    y_hat_left <- tibble(x = df_left[[.x]],
                         y_hat = lm_left$fitted.values,
                         algorithm = "jm")
    y_hat_right <- tibble(x = df_right[[.x]],
                          y_hat = lm_right$fitted.values,
                          algorithm = "jm")
    pred <- bind_rows(y_hat_left, y_hat_right)
    pct_slope_change <- 100*(lm_right$coefficients[1] - lm_left$coefficients[2]) /
        lm_left$coefficients[2]

    determinant_bp <- dplyr::if_else(pf_two < alpha_linearity &
                                         (pos_change == (pct_slope_change > 0)),
                                     TRUE, FALSE)

    bp_dat <- .data[min_ss_idx,] %>%
        mutate(algorithm = "jm",
               bp = bp,
               x_var = .x,
               y_var = .y,
               determinant_bp = determinant_bp,
               pct_slope_change = pct_slope_change,
               f_stat = f_stat,
               p_val_f = pf_two) %>%
        relocate(bp, algorithm, x_var, y_var, determinant_bp,
                 pct_slope_change, f_stat, p_val_f)

    return(list(breakpoint_data = bp_dat,
                fitted_vals = pred,
                lm_left = lm_left,
                lm_right = lm_right,
                lm_simple = lm_simple))

}

#' @keywords internal
loop_jm <- function(.data, .x, .y) {

    ss_both <- numeric(length = nrow(.data))
    for(i in 1:nrow(.data)) {
        if(i == 1 | i == nrow(.data)) {
            ss_both[i] <- NA
            next
        }
        # split data into left and right halves. x0 = i
        df_left <- .data[1:i,] # x <= x0
        df_right <- .data[i:nrow(.data),] # x > x0

        # the JM formula is y = b0 + b1*x0 + b3(x-x0). Therefore, we will create a new column in df_right that is equal to x - x0.
        df_right <- df_right %>%
            mutate(s1 = df_right[[.x]] - .data[[.x]][i])

        df_right$s1 = df_right[[.x]] - .data[[.x]][i]

        lm_left <- lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)

        # according to the JM algorithm, the right regression line will have a constant equal to b0 + b1*x0
        b0_plus_b1x0 <- lm_left$coefficients[1] + lm_left$coefficients[2] * .data[[.x]][i]

        # for some reason, using I() function to try to force the regression to use b0_plus_b1x0 as an intercept doesn't work, while the offset argument does. The slope coefficient (b3) is the same, but the predicted values are not.
        lm_right <- lm(df_right[[.y]] ~ 0 + s1, data = df_right,
                       offset = rep(b0_plus_b1x0, nrow(df_right)))

        # This should be the same as using the offset argument, but it's not
        # lm_right <- lm(I(.data[[.y]] - b0_plus_b1x0) ~ 0 + .data[[.x]], data = df_right)

        ss_left <- sum((lm_left$residuals)^2)
        ss_right <- sum((lm_right$residuals)^2)
        ss_both[i] <- (ss_left + ss_right) # is this off by 1?
    }
    ss_both[which(ss_both == 0)] <- NA
    ss_both
}
