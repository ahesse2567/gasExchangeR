#' Finding a breakpoint using Orr's 'bruteforce' method
#'
#' @param .data Gas exchange data.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param vo2 The name of the \code{vo2} variable.
#' @param vco2 The name of the \code{vco2} variable.
#' @param ve The name of the \code{ve} variable.
#' @param time The name of the \code{time} variable.
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simplier model.
#'
#' @return
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
                vo2 = "vo2",
                vco2 = "vco2",
                ve = "ve",
                time = "time",
                alpha_linearity = 0.05) {
    browser()
    lm_simple <- lm(.data[[.y]] ~ 1 + .data[[.x]])
    RSS_simple <- sum(resid(lm_simple)^2)

    ss <- loop_orr(.data = .data, .x = .x, .y = .y)
    min_ss_idx <- which.min(ss)

    df_left <- df_avg[1:ss_min_idx,]
    df_right <- df_avg[ss_min_idx:nrow(df_avg),]

    lm_left <- lm(vco2 ~ 1 + vo2_abs, data = df_left)
    lm_right <- lm(vco2 ~ 1 + vo2_abs, data = df_right)

    RSS_two <- sum(resid(lm_left)^2) + sum(resid(lm_right)^2)
    MSE_two <- RSS_two / (nrow(df_avg) - 4) # -4 b/c estimating 4 parameters
    F_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)

    pf_two <- pf(F_stat, df1 = 2, df2 = nrow(df_avg) - 4, lower.tail = FALSE)
    if(pf_two > alpha_linearity) {
        message("No breakpoint detected")
        return(NULL) # not sure if this should be stop() or something else
    }



    # (P < 0.01) reduction in the total sum of squares is achieved by the addition of         # the second and third line segments

    # now we need to find the closest data point to the breakpoint
}
#'
#' @keywords internal
loop_orr <- function(.data, .x, .y) {
    # browser()
    ss_both <- vector(length = nrow(.data)-2)
    for(i in 2:(nrow(.data)-1)) {
        # split data into left and right halves
        # should these share the same point? Or should, they be different by 1 point?
        df_left <- .data[1:i,]
        df_right <- .data[(i):nrow(.data),] # each side of the regression is contiguous

        # make linear models of the two halves
        lm_left <- lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)
        lm_right <- lm(df_right[[.y]] ~ 1 + df_right[[.x]], data = df_right)

        ss_left <- sum((lm_left$residuals)^2)
        ss_right <- sum((lm_right$residuals)^2)
        #ss_both <- c(ss_both, ss_left + ss_right)
        ss_both[i-1] <- (ss_left + ss_right)
    }
    return(ss_both)
}
