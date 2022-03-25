#' Finding a breakpoint using Beaver's V-slope algorithm
#'
#' V-slope uses the v-slope \emph{algorithm} to find a breakpoint. Note, the v-slope algorithm is different from the v-slope \emph{method} as a whole because the latter also includes a specific set of data processing steps (see reference). In order to use the v-slope algorithm as originally described by Beaver et al. (1986), pass \code{"v-slope"} to the \code{method} argument in the \code{breakpoint} function. Without indicating the \code{method} as \code{"v-slope"}, the function returns the breakpoint associated with ratio of the highest distance between the single regression line to the intersection point of the two regression lines, to the mean square error of the single regression line. This can sometimes lead to odd behavior and it is therefore generally recommended to specify \code{"v-slope"} in the \code{method} argument of \code{breakpoint} when using this function. See the Warnings and Details sections for details.
#'
#' @section Warning:
#' Using the v-slope algorithm as originally described by Beaver et al. (1986) is only intended to be used with the VCO2 vs. VO2 relationship when both are absolute measures of the same units (e.g. mL/min). The quality control criteria for this algorithm are an the increase in slope from the left to right regression line of > 0.1 and the slope of the left regression > 0.6. Using other x-y relationship leads to odd behavior that may not satisfy those criteria.
#'
#' @details
#' If the \code{method} argument is set to \code{"v-slope"}, the function enters a while loop to satisfy the quality control criteria (see Warnings section for details). Otherwise, the function returns the breakpoint as described the initial description.
#'
#' @param .data Gas exchange data.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param slope_change_lim The absolute amount by which the slope should change from the left to the right regression line. Default per Beaver (1986) is 0.1.
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.
#' @param bp Is this algorithm being used to find vt1 or vt2?
#' @param vo2 The name of the vo2 column in \code{.data}
#' @param vco2 The name of the vco2 column in \code{.data}
#' @param ve The name of the ve column in \code{.data}
#' @param time The name of the time column in \code{.data}
#' @param method Pass \code{"v-slope"} to the method argument to exactly reproduce the procedure from the original paper. See the Warnings section for details.
#' @param left_slope_lim The original paper requires that the left regression line have a slope of > \code{0.6}.
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
                    vo2 = "vo2",
                    vco2 = "vco2",
                    ve = "ve",
                    time = "time",
                    method = NULL,
                    slope_change_lim = 0.1,
                    left_slope_lim = 0.6,
                    alpha_linearity = 0.05,
                    bp) {
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])

    dist_MSE_ratio <- loop_v_slope(.data = .data, .x = .x, .y = .y)
    slope_change <- 0
    i <- 1
    bp_idx <- order(-dist_MSE_ratio)[i]
    df_left <- .data[1:bp_idx,]
    lm_left <- lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)

    # TODO the slope_change < slope_change_lim and left slope > 0.6 only apply
    # specifically when x = VO2 and y = VCO2. How do we check for that?
    # I think actually specifically apply when method == "v-slope"
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
                dplyr::mutate(bp = bp,
                              x_var = .x,
                              y_var = .y,
                              algorithm = "v-slope",
                              determinant_bp = FALSE,
                              pct_slope_change = NA,
                              f_stat = NA,
                              p_val_f = NA) %>%
                dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp) %>%
                purrr::map_df(num_to_na)

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

    pct_slope_change <- 100*(lm_right$coefficients[2] - lm_left$coefficients[2]) /
        lm_left$coefficients[2]

    determinant_bp <- dplyr::if_else(pf_two > alpha_linearity, FALSE, TRUE)

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
        dplyr::mutate(algorithm = "v_slope",
                      bp = bp,
                      x_var = .x,
                      y_var = .y,
                      determinant_bp = determinant_bp,
                      pct_slope_change = pct_slope_change,
                      f_stat = f_stat,
                      p_val_f = pf_two) %>%
        dplyr::select(-c(dist_x_sq, dist_y_sq, d)) %>%
        relocate(bp, algorithm, x_var, y_var, determinant_bp,
                 pct_slope_change, f_stat, p_val_f)

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
    dist_MSE_ratio <- vector(length = nrow(.data))

    for(i in 1:nrow(.data)) {
        if(i == 1 | i == nrow(.data)) {
            dist_MSE_ratio[i] <- NA
            next
        }

        df_left <- .data[1:i,] # split data into left half
        df_right <- .data[i:nrow(.data),] # split data into right half

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

        dist_MSE_ratio[i] <- d / anova(lm_simple)['Residuals', 'Mean Sq']
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
