#' Finding a breakpoint using Beaver's V-slope algorithm
#'
#' @param .data Gas exchange data.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param slope_change_lim The absolute amount by which the slope should change from the left to the right regression line. Default per Beaver (1986) is 0.1.
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
#' # TODO write an example
v_slope <- function(.data,
                    .x,
                    .y,
                    slope_change_lim = 0.1,
                    vo2 = "vo2",
                    vco2 = "vco2",
                    ve = "ve",
                    time = "time") {
    # browser()
    # TODO exclude data at the beginning if the VCO2 vs. VO2 slope is < 0.6
    ss <- loop_v_slope(.data = .data, .x = .x, .y = .y)
    slope_change <- 0
    i <- 1
    while(slope_change < slope_change_lim) {
        bp_idx <- order(-ss)[i]
        df_left <- df_avg[1:bp_idx,]
        df_right <- df_avg[(bp_idx+1):nrow(df_avg),]
        lm_left <- lm(vco2 ~ 1 + vo2_abs, data = df_left)
        lm_right <- lm(vco2 ~ 1 + vo2_abs, data = df_right)

        slope_change <- lm_right$coefficients[2] - lm_left$coefficients[2]
        i <- i + 1
        if (i > length(ss)) {
            stop("No breakpoint found. Change between slopes was never >= 0.1.")
        }
    }

    df_left <- .data[1:bp_idx,] # split data into left half
    df_right <- .data[(bp_idx+1):nrow(.data),] # split data into right half

    # make linear models of the two regressions
    lm_left <- lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)
    lm_right <- lm(df_right[[.y]] ~ 1 + df_right[[.x]], data = df_right)

    # find intersection point of left and right regressions
    lr_intersect <- intersection_point(lm_left, lm_right)

    # find closest data point to intersection point
    bp <- .data %>%
        mutate(dist_x_sq = (.data[[.x]] - lr_intersect["x"])^2,
               dist_y_sq = (.data[[.y]] - lr_intersect["y"])^2,
               d = sqrt(dist_x_sq + dist_y_sq)) %>%
        as_tibble() %>%
        arrange(d) %>%
        slice(1)

    bp
}

#' @keywords internal
loop_v_slope <- function(.data, .x, .y) {
    # browser()
    lm_simple <- lm(.data[[.y]] ~ 1 + .data[[.x]], data = .data)
    # find slope of line perpendicular to slope of lm_simple
    recip_slope <- (-1 / lm_simple$coefficients[2]) # used in for loop
    dist_MSE_ratio <- vector(length = nrow(.data)-3)
    for(i in 2:(nrow(.data)-2)) { # nrow(df)-2) b/c these DON'T share a point
        # the curve is divided into two regions, each of which is fitted by linear regression

        df_left <- .data[1:i,] # split data into left half
        df_right <- .data[(i+1):nrow(.data),] # split data into right half

        # make linear models of the two regressions
        lm_left <- lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)
        lm_right <- lm(df_right[[.y]] ~ 1 + df_right[[.x]], data = df_right)

        # find intersection point of left and right regressions
        lr_intersect <- intersection_point(lm_left, lm_right)

        b_recip <- recip_slope * (-1) * lr_intersect["x"] + lr_intersect["y"]

        x_simple_recip <- (b_recip - lm_simple$coefficients[1]) /
            (lm_simple$coefficients[2] - (-1 / lm_simple$coefficients[2]))
        y_simple_recip <- lm_simple$coefficients[1] +
            lm_simple$coefficients[2]*x_simple_recip

        d <- sqrt((x_simple_recip - lr_intersect["x"])^2 +
                      (y_simple_recip - lr_intersect["y"])^2)

        dist_MSE_ratio[i-1] <- d / anova(lm_simple)['Residuals', 'Mean Sq']
    }

    dist_MSE_ratio
}
