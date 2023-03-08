#' Finding a breakpoint using Beaver's V-slope algorithm
#'
#' V-slope uses the v-slope \emph{algorithm} to find a breakpoint. Note, the v-slope algorithm is different from the v-slope \emph{method} as a whole because the latter also includes a specific set of data processing steps (Beaver et al., 1986). In order to use the v-slope algorithm as originally described by Beaver et al. (1986), pass \code{"v-slope"} to the \code{method} argument in the \code{breakpoint} function. Without indicating the \code{method} as \code{"v-slope"}, the function returns the breakpoint associated with ratio of the highest distance between the single regression line to the intersection point of the two regression lines, to the mean square error of the piecewise regression lines. This can sometimes lead to odd behavior and it is therefore generally recommended to specify \code{"v-slope"} in the \code{method} argument of \code{breakpoint} when using this function. See the Warnings and Details sections for details.
#'
#' @section Warning:
#' Using the v-slope algorithm as originally described by Beaver et al. (1986) is only intended to be used with the VCO2 vs. VO2 relationship when both are absolute measures of the same units (e.g. mL/min). The quality control criteria for this algorithm are an the increase in slope from the left to right regression line of > 0.1 and the slope of the left regression > 0.6. Using other x-y relationship leads to odd behavior that may not satisfy those criteria.
#'
#' @details
#' If the \code{method} argument is set to \code{"v-slope"}, the function enters a while loop to satisfy the quality control criteria (see Warnings section for details). Otherwise, the function returns the breakpoint as described the initial description.
#'
#' According to the Beaver et al. (1986, p. 2023) paper, the "intersection between the two regression lines is the tentative AT point." This usually provides a reasonable solution, but in some cases, the intersection of the best-fit model lies beyond the range of the x-axis. Instead of returning that value, we report the solution with the highest ratio of distance between the intersection point of the left and right regressions and the single regression line to the mean square error of the piecewise regression that also satisfies the following criteria: 1) the change in slope from the left to the right regression matches the anticipated change, usually positive; 2) the slope of the right regression matches the anticipated slope, usually positive; 3) the p-value of the F-test based on the extra sums of squares principle is statistically significant; and 4) the x-coordinate of the intersection point is within the range of the observed x-values. This differs slightly from the original method, but provides potentially plausible answers when the alternative is a negative VO2 or a VO2 well beyond VO2max.
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
#' @param pos_change Do you expect the change in slope to be positive (default) or negative? If a two-line regression explains significantly reduces the sum square error but the change in slope does not match the expected underlying physiology, the breakpoint will be classified as indeterminate.
#' @param front_trim_vt1 How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param front_trim_vt2 See front_trim_vt1. This is here so this this function plays nice with `set_front_trim()`. The v-slope method should be used for finding VT1 only, so leave this alone.
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param pos_slope_after_bp Should the slope after the breakpoint be positive? Default is `TRUE`. This catches cases when the percent change in slope is positive, but the second slope is still negative. Change to `FALSE` when PetCO2 is the y-axis variable.
#'
#' @return A list including slice of the original data frame at the threshold index with new columns `algorithm`, `determinant_bp`, `pct_slope_change`, `f_stat`, and `p_val_f.` The list also includes the fitted values, the left and right sides of the piecewise regression, and a simple linear regression.
#' @importFrom rlang :=
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
                    bp,
                    ...,
                    vo2 = "vo2",
                    vco2 = "vco2",
                    ve = "ve",
                    time = "time",
                    method = NULL,
                    slope_change_lim = 0.1,
                    left_slope_lim = 0.6,
                    front_trim_vt1 = 60,
                    front_trim_vt2 = 60,
                    alpha_linearity = 0.05,
                    pos_change = TRUE,
                    pos_slope_after_bp = TRUE,
                    ordering = c("by_x", "time")
                    ) {
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
    if(.x != vo2 | .y != vco2) {
        warning("x variable may NOT be (absolute) VO2 or y variable may NOT be VCO2. This method is designed to work with x = absolute VO2 and y = VCO2")
    }

    ordering <- match.arg(ordering, several.ok = FALSE)
    .data <- order_cpet_df(.data, .x = .x , time = time,
                           ordering = ordering)
    plot_df <- .data

    front_trim <- set_front_trim(bp = bp,
                                 front_trim_vt1 = front_trim_vt1,
                                 front_trim_vt2 = front_trim_vt2)
    .data <- .data %>%
        dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))

    dist_MSE_ratio_idx <- loop_v_slope(.data = .data, .x = .x, .y = .y,
                                       pos_change, pos_slope_after_bp)
    slope_change <- 0
    i <- 1
    bp_idx <- dist_MSE_ratio_idx[i]
    df_left <- .data[1:bp_idx,]
    lm_left <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::lm(data = df_left)

    n_rows <- nrow(.data)

    # browser()
    while((slope_change < slope_change_lim) |
          (lm_left$coefficients[2] < left_slope_lim)) {
        # enter this loop if slopes did not change be enough or if the left slope
        # was too low
        bp_idx <- dist_MSE_ratio_idx[i] # find the next-best fit
        df_left <- .data[1:bp_idx,]
        df_right <- .data[bp_idx:n_rows,]
        lm_left <- stats::lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)
        lm_right <- stats::lm(df_right[[.y]] ~ 1 + df_right[[.x]], data = df_right)

        slope_change <- lm_right$coefficients[2] - lm_left$coefficients[2]

        if (i > length(dist_MSE_ratio_idx[!is.na(dist_MSE_ratio_idx)])) {
            message(paste0("No VT1 breakpoint found with v-slope method because the change between slopes was never >= ",
                           slope_change_lim, "."))
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
                        # fitted_vals = pred, # TODO how to return fitted values?
                        lm_left = NULL,
                        lm_right = NULL,
                        lm_simple = NULL))
        }
        i <- i + 1
    }

    df_left <- .data[1:bp_idx,] # split data into left portion
    df_right <- .data[(bp_idx+1):n_rows,] # split data into right portion

    # make linear models of the two regressions
    lm_left <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::lm(data = df_left)
    lm_right <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::lm(data = df_right)
    # simple linear regression
    lm_simple <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::lm(data = .data)

    # check for a significant departure from linearity
    pw_stats <- piecewise_stats(lm_left, lm_right, lm_simple)
    list2env(pw_stats, envir = environment())

    pct_slope_change <- 100*(lm_right$coefficients[2] - lm_left$coefficients[2]) /
        abs(lm_left$coefficients[2])

    determinant_bp <- check_if_determinant_bp(p = pf_two,
                                              pct_slope_change = pct_slope_change,
                                              pos_change = pos_change,
                                              pos_slope_after_bp =
                                                  pos_slope_after_bp,
                                              slope_after_bp = stats::coef(lm_right)[2],
                                              alpha = alpha_linearity)

    y_hat_left <- tibble::tibble("{.x}" := df_left[[.x]],
                         "{.y}" := lm_left$fitted.values,
                         algorithm = "v-slope")
    y_hat_right <- tibble::tibble("{.x}" := df_right[[.x]],
                          "{.y}" := lm_right$fitted.values,
                          algorithm = "v-slope")
    pred <- dplyr::bind_rows(y_hat_left, y_hat_right)

    # find intersection point of left and right regressions
    int_point <- intersection_point(lm_left, lm_right)
    # find closest data point to intersection point and prepare output
    bp_dat <- find_threshold_vals(.data = .data, thr_x = int_point["x"],
                                  thr_y = int_point["y"], .x = .x,
                                  .y = .y, ...)
    bp_dat <- bp_dat %>%
        dplyr::mutate(algorithm = "v-slope",
                      x_var = .x,
                      y_var = .y,
                      determinant_bp = determinant_bp,
                      bp = bp,
                      pct_slope_change = pct_slope_change,
                      f_stat = f_stat,
                      p_val_f = pf_two) %>%
        dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp,
                        pct_slope_change, f_stat, p_val_f)

    bp_plot <- make_piecewise_bp_plot(plot_df, .x, .y, lm_left, lm_right, bp_dat)

    return(list(breakpoint_data = bp_dat,
                fitted_vals = pred,
                lm_left = lm_left,
                lm_right = lm_right,
                lm_simple = lm_simple,
                bp_plot = bp_plot))
}

#' @keywords internal
loop_v_slope <- function(.data, .x, .y, pos_change, pos_slope_after_bp) {

    n_rows <- nrow(.data)

    # initialize empty vectors
    ss_left <- ss_right <- ss_both <- RSS_two <- MSE_two <- f_stat <- pf_two <-
        pct_slope_change <- int_point_x <- dist_MSE_ratio <-  numeric(n_rows)
    pos_change_vec <- pos_slope_after_bp_vec <- logical(n_rows)

    # browser()
    lm_simple <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::lm(data = .data)
    RSS_simple <- sum(stats::resid(lm_simple)^2)
    recip_slope <- (-1 / lm_simple$coefficients[2]) # used in for loop
    # find slope of line perpendicular to slope of lm_simple
    # There's an easier way to find the distance between a point and a line
    # by putting the equation of the simple regression into the standard form
    # of Ax + By + C = 0
    # https://www.mathportal.org/calculators/analytic-geometry/line-point-distance.php
    # https://brilliant.org/wiki/dot-product-distance-between-point-and-a-line/

    for(i in 1:n_rows) {
        if(i == 1 | i == n_rows) {
            ss_left[i] <- NA
            ss_right[i] <- NA
            ss_both[i] <- NA
            RSS_two[i] <- NA
            MSE_two[i] <- NA
            f_stat[i] <- NA
            pf_two[i] <- NA
            pos_change_vec[i] <- NA
            pos_slope_after_bp_vec[i] <- NA
            int_point_x[i] <- NA
            dist_MSE_ratio[i] <- NA
            next
        }

        df_left <- .data[1:i,] # split data into left half
        df_right <- .data[i:n_rows,] # split data into right half

        # make linear models of the two regressions
        lm_left <- paste0(.y, " ~ ", "1 + ", .x) %>%
            stats::lm(data = df_left)
        lm_right <- paste0(.y, " ~ ", "1 + ", .x) %>%
            stats::lm(data = df_right)

        ss_left[i] <- sum((lm_left$residuals)^2)
        ss_right[i] <- sum((lm_right$residuals)^2)
        ss_both[i] <- ss_left[i] + ss_right[i]

        RSS_two[i] <- sum(stats::resid(lm_left)^2) + sum(stats::resid(lm_right)^2)
        MSE_two[i] <- RSS_two[i] / (nrow(lm_simple$model) - 4) # -4 b/c estimating 4 parameters
        f_stat[i] <- (RSS_simple - RSS_two[i]) / (2 * MSE_two[i])
        pf_two[i] <- stats::pf(f_stat[i], df1 = 2, df2 = nrow(lm_simple$model) - 4,
                               lower.tail = FALSE)

        pct_slope_change[i] <- 100*(lm_right$coefficients[2] - lm_left$coefficients[2]) /
            abs(lm_left$coefficients[2])
        pos_change_vec[i] <- if_else(pct_slope_change[i] > 0, TRUE, FALSE)
        pos_slope_after_bp_vec[i] <- if_else(lm_right$coefficients[2] > 0, TRUE, FALSE)

        # Calculate MSE according to the JM algorithm: MSE = ss_both / (n - 4)
        # we're subtracting for because we're estimating 4 parameters: two slopes
        # and two intercepts
        MSE <- ss_both[i] / (n_rows - 4)

        # find intersection point of left and right regressions
        lr_intersect <- intersection_point(lm_left, lm_right)
        int_point_x[i] <- lr_intersect["x"]

        b_recip <- recip_slope * (-1) * lr_intersect["x"] + lr_intersect["y"]

        x_simple_recip <- (b_recip - lm_simple$coefficients[1]) /
            (lm_simple$coefficients[2] - (-1 / lm_simple$coefficients[2]))
        y_simple_recip <- lm_simple$coefficients[1] +
            lm_simple$coefficients[2]*x_simple_recip

        d <- sqrt((x_simple_recip - lr_intersect["x"])^2 +
                      (y_simple_recip - lr_intersect["y"])^2)

        dist_MSE_ratio[i] <- d / MSE
    }

    v_slope_stats <- tibble::tibble(p = pf_two,
                                pos_change = pos_change,
                                pos_slope_after_bp = pos_slope_after_bp,
                                dist_MSE_ratio = dist_MSE_ratio,
                                int_point_x = int_point_x,
                                # inside_F95 = inside_F95
    ) %>%
        dplyr::mutate(idx = dplyr::row_number())

    range_x <- range(.data[[.x]])

    dist_MSE_ratio_idx <- v_slope_stats %>%
        dplyr::filter(p < 0.05 &
                          pos_change_vec == TRUE &
                          pos_slope_after_bp_vec == TRUE &
                          dplyr::between(int_point_x,
                                         range_x[1],
                                         range_x[2])) %>%
        arrange(desc(dist_MSE_ratio)) %>%
        select(idx) %>%
        pull()

    if(length(dist_MSE_ratio_idx) > 0) {
        return(dist_MSE_ratio_idx)
    } else {
        return(2)
    }
}

#' @keywords internal
num_to_na <- function(x) {
    if(is.numeric(x)) {
        x <- NA
    }
    x
}
