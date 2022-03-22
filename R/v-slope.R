#' Finding a breakpoint using Beaver's V-slope algorithm
#'
#' @param .data Gas exchange data.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param slope_change_lim The absolute amount by which the slope should change from the left to the right regression line. Default per Beaver (1986) is 0.1.
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.
#' @param bp Is this algorithm being used to find vt1 or vt2?
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
                    alpha_linearity = 0.05,
                    bp) {
    # browser()
    # TODO exclude data at the beginning if the VCO2 vs. VO2 slope is < 0.6
    dist_MSE_ratio <- loop_v_slope(.data = .data, .x = .x, .y = .y)
    slope_change <- 0
    i <- 1
    while(slope_change < slope_change_lim) {
        bp_idx <- order(-dist_MSE_ratio)[i]
        df_left <- .data[1:bp_idx,]
        df_right <- .data[(bp_idx+1):nrow(.data),]
        lm_left <- lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)
        lm_right <- lm(df_right[[.y]] ~ 1 + df_right[[.x]], data = df_right)

        slope_change <- lm_right$coefficients[2] - lm_left$coefficients[2]
        i <- i + 1
        if (i > length(dist_MSE_ratio)) {
            message("No breakpoint found. Change between slopes was never >= 0.1.")
            bp_dat <- .data %>%
                dplyr::slice(1) %>%
                dplyr::select(time, vo2, vco2, ve) %>%
                dplyr::mutate(bp = bp,
                       method = "v-slope",
                       determinant_bp = FALSE,
                       pct_slope_change = NA,
                       f_stat = NA,
                       p_val_f = NA) %>%
                dplyr::relocate(bp, method, determinant_bp) %>%
                map_df(num_to_na)

            return(list(breakpoint_data = bp_dat,
                        # fitted_vals = pred, # dTODO how to return fitted values?
                        lm_left = NULL,
                        lm_right = NULL,
                        lm_simple = NULL))
        }
    }

    df_left <- .data[1:bp_idx,] # split data into left half
    df_right <- .data[(bp_idx+1):nrow(.data),] # split data into right half

    # make linear models of the two regressions
    lm_left <- lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)
    lm_right <- lm(df_right[[.y]] ~ 1 + df_right[[.x]], data = df_right)
    # sinle line regression
    lm_simple <- lm(.data[[.y]] ~ 1 + .data[[.x]], data = .data)

    # check for a significant departure from linearity
    RSS_simple <- sum(resid(lm_simple)^2)
    RSS_two <- sum(resid(lm_left)^2) + sum(resid(lm_right)^2)
    MSE_two <- RSS_two / (nrow(.data) - 4) # -4 b/c estimating 4 parameters
    f_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)
    pf_two <- pf(f_stat, df1 = 2, df2 = nrow(.data) - 4, lower.tail = FALSE)
    determinant_bp <- dplyr::if_else(pf_two > alpha_linearity, FALSE, TRUE)

    pct_slope_change <- 100*(lm_right$coefficients[2] - lm_left$coefficients[2]) /
        lm_left$coefficients[2]

    # find intersection point of left and right regressions
    lr_intersect <- intersection_point(lm_left, lm_right)
    # find closest data point to intersection point and prepare output
    bp_dat <- .data %>%
        dplyr::mutate(dist_x_sq = (.data[[.x]] - lr_intersect["x"])^2,
               dist_y_sq = (.data[[.y]] - lr_intersect["y"])^2,
               d = sqrt(dist_x_sq + dist_y_sq)) %>%
        tibble::as_tibble() %>%
        dplyr::arrange(d) %>%
        dplyr::slice(1) %>%
        dplyr::mutate(method = "v_slope",
                      determinant_bp = determinant_bp,
                      bp = bp) %>%
        dplyr::select(bp, method, determinant_bp, time, vo2, vco2, ve) %>%
        mutate(pct_slope_change = pct_slope_change,
           f_stat = f_stat,
           p_val_f = pf_two)

    return(list(breakpoint_data = bp_dat,
                # fitted_vals = pred, # TODO how to return fitted values?
                lm_left = lm_left,
                lm_right = lm_right,
                lm_simple = lm_simple))
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

#' @keywords internal
num_to_na <- function(x) {
    if(is.numeric(x)) {
        x <- NA
    }
    x
}
