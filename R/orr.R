#' Finding a breakpoint using Orr's 'bruteforce' method
#'
#' @param .data Gas exchange data.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param vo2 The name of the \code{vo2} variable.
#' @param vco2 The name of the \code{vco2} variable.
#' @param ve The name of the \code{ve} variable.
#' @param time The name of the \code{time} variable.
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
                time = "time") {
    browser()
    lm_simple <- lm(.data[[.y]] ~ 1 + .data[[.x]])
    ss <- loop_orr(.data = .data, .x = .x, .y = .y)

    min_ss_idx <- which.min(ss)
    min_ss_idx

    # get sst and sse for simple and the two-line regression
    # compare those two regressions with an F test.
    # does that F test require the total sums of squares, or just the mse?
    # how can I do the likelihood ratio test by hand?
    # An analysis of variance determines whether a significant
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
