#' Finding a breakpoint using Beaver's V-slope algorithm
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
#' Beaver, W. L., Wasserman, K., & Whipp, B. J. (1986). A new method for detecting anaerobic threshold by gas exchange. Journal of Applied Physiology, 121(6), 2020–2027.
#'
#' @examples
v_slope <- function(.data,
                    .x,
                    .y,
                    vo2 = "vo2",
                    vco2 = "vco2",
                    ve = "ve",
                    time = "time") {
    browser()
    lm_simple <- lm(.data[[.y]] ~ 1 + .data[[.x]])


}

loop_v_slope <- function(.data, .x, .y) {
    lm_simple <- lm(.data[[.y]] ~ 1 + .data[[.x]], data = .data)
    dist_MSE_ratio <- vector(length = nrow(.data)-3)
    for(i in 2:(nrow(.data)-2)) { # nrow(df)-2) b/c these DON'T share a point
        # the curve is divided into two regions, each of which is fitted by linear regression

        df_left <- .data[1:i,] # split data into left half
        df_right <- .data[(i+1):nrow(.data),] # split data into right half

        # make linear models of the two regressions
        lm_left <- lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)
        lm_right <- lm(df_right[[.y]] ~ 1 + df_right[[.x]], data = df_right)

        # find intersection point of left and right regressions
        # x = (b2 - b1) / (m1 - m2)
        x_left_right <- (lm_right$coefficients[1] - lm_left$coefficients[1]) /
            (lm_left$coefficients[2] - lm_right$coefficients[2])
        y_left_right <- lm_left$coefficients[1] + lm_left$coefficients[2]*x_left_right

        # find slope of line perpandicular to slope of lm_simple
        recip_slope <- (-1 / lm_simple$coefficients[2])
        b_recip <- recip_slope * (-1) * x_left_right + y_left_right

        x_simple_recip <- (b_recip - lm_simple$coefficients[1]) /
            (lm_simple$coefficients[2] - (-1 / lm_simple$coefficients[2]))
        y_simple_recip <- lm_simple$coefficients[1] +
            lm_simple$coefficients[2]*x_simple_recip

        d <- sqrt((x_simple_recip - x_left_right)^2 +
                      (y_simple_recip - y_left_right)^2)

        #######
        # x = (b2 - b1) / (m1 - m2)
        # treat b2 and m2 as from the simple regression line
        x_left_right <- (lm_right$coefficients[1] - lm_left$coefficients[1]) /
            (lm_left$coefficients[2] - lm_right$coefficients[2])
        y_left_right <- lm_left$coefficients[1] + lm_left$coefficients[2]*x_left_right

        x_left_simple <- (lm_simple$coefficients[1] - lm_left$coefficients[1]) /
            (lm_left$coefficients[2] - lm_simple$coefficients[2])
        y_left_simple <- lm_left$coefficients[1] + lm_left$coefficients[2]*x_left_simple

        x_right_simple <- (lm_simple$coefficients[1] - lm_right$coefficients[1]) /
            (lm_right$coefficients[2] - lm_simple$coefficients[2])
        y_right_simple <- lm_left$coefficients[1] + lm_left$coefficients[2]*x_right_simple

        dist_left_simple <- sqrt((x_left_simple - x_left_right)^2 +
                                     (y_left_simple - y_left_right)^2)

        dist_right_simple <- sqrt((x_right_simple - x_left_right)^2 +
                                      (y_right_simple - y_left_right)^2)

        if(dist_left_simple >= dist_right_simple) {
            dist <- dist_left_simple
        } else {
            dist <- dist_right_simple
        }

        # we need to calculate the mean square error of the regression, but the Beaver paper doesn't indicate which regression.
        #Is it the MSE of the single regression line? Or is it the MSE of the left and right regression lines?
        #Let's assume for now that it's the single regression line b/c that's easier

        deg_free <- nrow(df) - 2 # we're subtracting 2 b/c we're estimating 2 parameters: the intercept and the slope
        MSE_simple <- sum(lm_simple$residuals^2) / deg_free
        dist_MSE_ratio[i-1] <- (dist / MSE_simple)
        #dist_MSE_ratio <- c(dist_MSE_ratio, dist / MSE_simple)
    }

    return(dist_MSE_ratio)
}
