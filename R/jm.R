#' Finding a breakpoint using the Jones-Molitoris method.
#'
#' @param .data Gas exchange data.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param vo2 The name of the \code{vo2} variable.
#' @param vco2 The name of the \code{vco2} variable.
#' @param ve The name of the \code{ve} variable.
#' @param time The name of the \code{time} variable.
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.
#' @param bp Is this algorithm being used to find vt1 or vt2?
#'
#' @return
#' @export
#'
#' @references
#' Jones, R. H., & Molitoris, B. A. (1984). A statistical method for determining the breakpoint of two lines. Analytical Biochemistry, 141(1), 287–290. https://doi.org/10.1016/0003-2697(84)90458-5
#'
#' @examples
#' #' # TODO write examples
jm <- function(.data,
               .x,
               .y,
               vo2 = "vo2",
               vco2 = "vco2",
               ve = "ve",
               time = "time",
               alpha_linearity = 0.05,
               bp) {
    # TODO add which bp argument and determinate/indeterminate output
    # browser()
    ss <- loop_jm(.data = .data, .x = .x, .y = .y)
    min_ss_idx <- which.min(ss)

    # I don't think this actually intersects at x0 right now!!!!

    df_left <- .data[1:min_ss_idx,] # x < x0
    df_right <- .data[(min_ss_idx+1):nrow(.data),] # x >= x0

    x_knot <- .data[[.x]][min_ss_idx+1]

    df_right <- df_right %>%
        mutate(s1 = df_right[[.x]] - x_knot)

    lm_left <- lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)
    # according to the JM method, the right regression line will have a constant equal to b0 + b1*x0
    b0_plus_b1x0 <- lm_left$coefficients[1] + lm_left$coefficients[2] * .data[[.x]][min_ss_idx+1]

    # for some reason, using I() function to try to force the regression to use b0_plus_b1x0 as an intercept doesn't work, while the offset argument does. The slope coefficient (b3) is the same, but the predicted values are not.
    lm_right <- lm(df_right[[.y]] ~ 0 + s1, data = df_right,
                   offset = rep(b0_plus_b1x0, nrow(df_right)))

    lm_simple <- lm(.data[[.y]] ~ 1 + .data[[.x]])

    RSS_simple <- sum(resid(lm_simple)^2)
    RSS_two <- sum(resid(lm_left)^2) + sum(resid(lm_right)^2)
    MSE_two <- RSS_two / (nrow(df_avg) - 4) # -4 b/c estimating 4 parameters
    f_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)
    pf_two <- pf(f_stat, df1 = 2, df2 = nrow(df_avg) - 4, lower.tail = FALSE)

    y_hat_left <- tibble(x = df_left[[.x]],
                         y_hat = lm_left$fitted.values,
                         method = "jm")
    y_hat_right <- tibble(x = df_right[[.x]],
                          y_hat = lm_right$fitted.values,
                          method = "jm")
    pred <- bind_rows(y_hat_left, y_hat_right)
    pct_slope_change <- 100*(lm_right$coefficients[1] - lm_left$coefficients[2]) /
        lm_left$coefficients[2]

    bp_dat <- .data[min_ss_idx+1,] %>%
        select(time, vo2, vco2, ve) %>%
        mutate(method = "jm",
               bp = bp,
               pct_slope_change = pct_slope_change,
               f_stat = f_stat,
               p_val_f = pf_two) %>%
        relocate(bp, method)

    return(list(breakpoint_data = bp_dat,
                fitted_vals = pred,
                lm_left = lm_left,
                lm_right = lm_right,
                lm_simple = lm_simple))

}

#' @keywords internal
loop_jm <- function(.data,
                    .x,
                    .y) {
    # browser()
    # rearranging data by the x variables was maybe a bad idea?
    # .data <- .data %>%
    #     arrange(.data[[.x]], .data[[time]])

    ss_both <- numeric(length = nrow(.data)-2)
    for(i in 2:(nrow(.data)-1)) {
        # algorithm data into left and right halves. x0 = i+1
        .data_left <- .data[1:i,] # x < x0
        .data_right <- .data[(i+1):nrow(.data),] # x >= x0

        # the JM formula is y = b0 + b1*x0 + b3(x-x0). Therefore, we will create a new column in .data_right that is equal to x - x0.
        .data_right <- .data_right %>%
            mutate(s1 = .data_right[[.x]] - .data[[.x]][i+1])

        .data_right$s1 = .data_right[[.x]] - .data[[.x]][i+1]

        lm_left <- lm(.data_left[[.y]] ~ 1 + .data_left[[.x]], data = .data_left)

        # according to the JM method, the right regression line will have a constant equal to b0 + b1*x0
        b0_plus_b1x0 <- lm_left$coefficients[1] + lm_left$coefficients[2] * .data[[.x]][i+1]

        # for some reason, using I() function to try to force the regression to use b0_plus_b1x0 as an intercept doesn't work, while the offset argument does. The slope coefficient (b3) is the same, but the predicted values are not.
        lm_right <- lm(.data_right[[.y]] ~ 0 + s1, data = .data_right,
                       offset = rep(b0_plus_b1x0, nrow(.data_right)))

        # This should be the same as using the offset argument, but it's not
        # lm_right <- lm(I(.data[[.y]] - b0_plus_b1x0) ~ 0 + .data[[.x]], data = .data_right)

        ss_left <- sum((lm_left$residuals)^2)
        ss_right <- sum((lm_right$residuals)^2)
        ss_both[i-1] <- (ss_left + ss_right)

        #ss_both <- c(ss_both, ss_left + ss_right)
    }

    ss_both
}
