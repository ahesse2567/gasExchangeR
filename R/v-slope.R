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
#' @param ci Should the output include confidence interval data? Default is `FALSE`.
#' @param conf_level Confidence level to use if calculating confidence intervals.
#' @param plots Should this function generate plots? Set to `FALSE` to save time.
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
                    vo2 = "vo2",
                    vco2 = "vco2",
                    ve = "ve",
                    time = "time",
                    ordering = c("by_x", "time"),
                    method = NULL, # can I get rid of this?
                    slope_change_lim = 0.1,
                    left_slope_lim = 0.6,
                    front_trim_vt1 = 60,
                    front_trim_vt2 = 60,
                    alpha_linearity = 0.05,
                    pos_change = TRUE,
                    pos_slope_after_bp = TRUE,
                    ci = FALSE,
                    conf_level = 0.95,
                    plots = TRUE,
                    ...
                    ) {
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
    if(.x != vo2 | .y != vco2) {
        warning("x variable may NOT be (absolute) VO2 or y variable may NOT be VCO2. This method is designed to work with x = absolute VO2 and y = VCO2")
    }

    ordering <- match.arg(ordering, several.ok = FALSE)
    .data <- order_cpet_df(.data, .x = .x , time = time,
                           ordering = ordering)
    plot_df <- .data

    # for now, we're ignoring the initial segment < 0.6 b/c it's unclear
    # what counted as the "initial segment". Did they create linear regressions
    # iteratively from the first data point through all the remaining ones and
    # filter the slope that way?
    # calculate all possible slopes to find initial segments of < 0.6
    # slopes <- intercepts <- numeric(length = nrow(plot_df))
    # slopes[1] <- intercepts[1] <-  NA
    # for(i in 2:length(slopes)) {
    #     lm_initial_segment <- paste0(.y, " ~ ", "1 + ", .x) %>%
    #         stats::lm(data = plot_df[1:i,])
    #     slopes[i] <- coef(lm_initial_segment)[2]
    #     intercepts[i] <- coef(lm_initial_segment)[1]
    # }
    #
    # left_slope_idx <- 1
    #
    # for(i in 2:length(slopes)) {
    #     color = dplyr::if_else(slopes[i] > 0.6, "black", "red")
    #     abline(a = intercepts[i], b = slopes[i], col = color)
    # }
    #
    # for (i in seq_along(slopes)) {
    #     if (all(slopes[(i+1):length(slopes)] > 0.6)) {
    #         left_slope_idx <- i + 1
    #         break
    #     }
    # }
    #
    # remove data from the left side if "any initial segment of the curve above
    # [one minute has] a slope < 0.6" (Beaver et al., 1986, p. 2023)
    # .data <- .data[left_slope_idx:nrow(.data),]


    front_trim <- set_front_trim(bp = bp,
                                 front_trim_vt1 = front_trim_vt1,
                                 front_trim_vt2 = front_trim_vt2)
    .data <- .data %>%
        dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))

    loop_res <- loop_v_slope(.data = .data,
                             .x = .x,
                             .y = .y,
                             alpha_linearity = alpha_linearity,
                             conf_level = conf_level)

    # return quick summary if generating models fails
    if(is.null(loop_res)) {
        bp_dat <- return_indeterminant_findings(
            bp = bp,
            algorithm = as.character(match.call()[[1]]),
            .x = .x,
            .y = .y,
            est_ci = "estimate")

        return(list(breakpoint_data = bp_dat))
    } else {
        best_idx <- numeric()
    }

    # need to go by dist_MSE_ratio? Basically how does this interact with while
    # loop()
    range_x <- range(.data[[.x]])

    dist_MSE_ratio_idx <- loop_res %>%
        dplyr::filter(p < alpha_linearity &
                          pos_change == pos_change &
                          pos_slope_after_bp == pos_slope_after_bp &
                          dplyr::between(int_point_x,
                                         range_x[1],
                                         range_x[2]) &
                          inside_ci) %>%
        dplyr::arrange(dplyr::desc(dist_MSE_ratio)) %>%
        dplyr::select(idx) %>%
        dplyr::pull()

    slope_change <- 0
    i <- 1
    # bp_idx <- dist_MSE_ratio_idx[i]
    # df_left <- .data[1:bp_idx,]
    # lm_left <- paste0(.y, " ~ ", "1 + ", .x) %>%
    #     stats::lm(data = df_left)

    n_rows <- nrow(.data)

    # browser()
    while(slope_change < slope_change_lim) {
        # enter this loop if slopes did not change be enough or if the left slope
        # was too low
        # check that we haven't run out of possible options
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
        bp_idx <- dist_MSE_ratio_idx[i] # find the next-best fit
        df_left <- .data[1:bp_idx,]
        df_right <- .data[bp_idx:n_rows,]
        lm_left <- stats::lm(df_left[[.y]] ~ 1 + df_left[[.x]], data = df_left)
        lm_right <- stats::lm(df_right[[.y]] ~ 1 + df_right[[.x]], data = df_right)

        slope_change <- lm_right$coefficients[2] - lm_left$coefficients[2]

        i <- i + 1
    }

    # THIS ISN'T GIVING THE SAME INT POINT X AS THE LOOP
    estimate_res <- get_v_slope_res(.data = .data,
                                    bp_idx = bp_idx,
                                    .x = .x,
                                    .y = .y,
                                    bp = bp,
                                    alpha_linearity = alpha_linearity,
                                    pos_change = pos_change,
                                    pos_slope_after_bp = pos_slope_after_bp,
                                    est_ci = "estimate")

    # browser()

    if(ci) {

        # broom::tidy(estimate_res$lm_left)
        # broom::tidy(estimate_res$lm_right)
        # conf_stats_left <- stats::confint(estimate_res$lm_left)
        # conf_stats_left
        # conf_stats_right <- stats::confint(estimate_res$lm_right)
        # conf_stats_right
        #
        # plot(.data[[.x]], .data[[.y]])
        # plot_lims <- 3000
        # plot(.data[[.x]], .data[[.y]],
        #      xlim = c(1600, 2000),
        #      ylim = c(1400, 1600))
        # abline(lm_left, col = "blue") # add estimated left regression
        # abline(lm_right, col = "blue") # add estimated right regression
        # x_est <- intersection_point(lm_left, lm_right)["x"]
        # abline(v = x_est, col = "blue")
        # # from left regression
        # # add lower intercept and lower slope
        # abline(a = conf_stats_left[1,1], b = conf_stats_left[2,1])
        # # add lower intercept and higher slope
        # abline(a = conf_stats_left[1,1], b = conf_stats_left[2,2])
        # # add higher intercept and lower slope
        # abline(a = conf_stats_left[1,2], b = conf_stats_left[2,1])
        # # add higher intercept and higher slope
        # abline(a = conf_stats_left[1,2], b = conf_stats_left[2,2])
        #
        # # from right regression
        # # add lower intercept and lower slope
        # abline(a = conf_stats_right[1,1],
        #        b = conf_stats_right[2,1], col = "red")
        # # add lower intercept and higher slope
        # abline(a = conf_stats_right[1,1],
        #        b = conf_stats_right[2,2], col = "red")
        # # add higher intercept and lower slope
        # abline(a = conf_stats_right[1,2],
        #        b = conf_stats_right[2,1], col = "red")
        # # add higher intercept and higher slope
        # abline(a = conf_stats_right[1,2],
        #        b = conf_stats_right[2,2], col = "red")
        #
        # x_lower <- (conf_stats_left[1,1] - conf_stats_right[1,2]) /
        #     (conf_stats_right[2,2] -
        #          conf_stats_left[2,1])
        # x_lower
        # x_upper <- (conf_stats_left[1,2] - conf_stats_right[1,1]) /
        #     (conf_stats_right[2,1] -
        #          conf_stats_left[2,2])
        # x_upper
        # x_est <- intersection_point(lm_left, lm_right)["x"]
        #
        # abline(v = x_lower)
        # abline(v = x_upper)
        # abline(v = x_est, col = "blue")
        #
        # # it's still a little anecdotal, but it seems like the
        # # confidence interval using the combinations of lowest and highest
        # # confidence intervals for the slopes and intercepts of the
        # # two regression lines gives estimates that extend beyond the
        # # range of the x-axis
        #
        # # what about just using the different intercepts and not the slopes
        # # to create the confidence interval?
        #
        # plot(.data[[.x]], .data[[.y]])
        # # plot estimate lines
        # abline(lm_left, col = "blue")
        # abline(lm_right, col = "blue")
        # # left regression lines
        # # lower intercept
        # abline(a = conf_stats_left[1,1],
        #        b = broom::tidy(estimate_res$lm_left)$estimate[2])
        # # upper intercept
        # abline(a = conf_stats_left[1,2],
        #        b = broom::tidy(estimate_res$lm_left)$estimate[2])
        # # right regression lines
        # # lower intercept
        # abline(a = conf_stats_right[1,1],
        #        b = broom::tidy(estimate_res$lm_right)$estimate[2],
        #        col = "red")
        # # upper intercept
        # abline(a = conf_stats_right[1,2],
        #        b = broom::tidy(estimate_res$lm_right)$estimate[2],
        #        col = "red")
        #
        # x_lower <- (conf_stats_left[1,1] - conf_stats_right[1,2]) /
        #     (broom::tidy(estimate_res$lm_right)$estimate[2] -
        #          broom::tidy(estimate_res$lm_left)$estimate[2])
        # x_lower
        # x_upper <- (conf_stats_left[1,2] - conf_stats_right[1,1]) /
        #     (broom::tidy(estimate_res$lm_right)$estimate[2] -
        #          broom::tidy(estimate_res$lm_left)$estimate[2])
        # x_upper
        # x_est <- intersection_point(lm_left, lm_right)["x"]
        #
        # abline(v = x_lower)
        # abline(v = x_upper)
        # abline(v = x_est, col = "blue")
        #
        # # I suspect the following
        # # farthest left: left lm lower lower, right lm higher higher
        # # farthest right: lm left higher higher, lm right lower lower
        # # but maybe we should just calculate all combinations and search
        # # for the lowest and highest values, respectively
        #
        # ?abline
        # crit_t_left <- qt(1 - 0.05 / 2, estimate_res$lm_left$df.residual)
        # stats::confint(estimate_res$lm_left)
        #
        # std_err_int_left <- broom::tidy(estimate_res$lm_left)$std.error[1]
        # broom::tidy(estimate_res$lm_left)$estimate[1] +
        #     c(-1, 1) * crit_t_left * std_err_int_left
        #
        #
        # broom::tidy(estimate_res$lm_left)$std.error[1]
        # broom::glance(estimate_res$lm_left) %>% View
        # broom::augment(estimate_res$lm_left) %>% View

        ci_lower_idx <- loop_res %>%
            dplyr::filter(inside_ci) %>%
            dplyr::filter(int_point_x == min(int_point_x)) %>%
            dplyr::select(idx) %>%
            dplyr::pull()

        lower_ci_res <- get_v_slope_res(
            .data = .data,
            bp_idx = ci_lower_idx,
            .x = .x,
            .y = .y,
            bp = bp,
            alpha_linearity = alpha_linearity,
            pos_change = pos_change,
            pos_slope_after_bp = pos_slope_after_bp,
            est_ci = "lower_ci")

        ci_upper_idx <- loop_res %>%
            dplyr::filter(inside_ci) %>%
            dplyr::filter(int_point_x == max(int_point_x)) %>%
            dplyr::select(idx) %>%
            dplyr::pull()

        upper_ci_res <- get_v_slope_res(.data = .data,
                                    bp_idx = ci_upper_idx,
                                    .x = .x,
                                    .y = .y,
                                    bp = bp,
                                    alpha_linearity = alpha_linearity,
                                    pos_change = pos_change,
                                    pos_slope_after_bp = pos_slope_after_bp,
                                    est_ci = "upper_ci")

        # combine estimate and both CI breakpoint res into one tibble
        estimate_res$bp_dat <- dplyr::bind_rows(lower_ci_res$bp_dat,
                                     estimate_res$bp_dat,
                                     upper_ci_res$bp_dat)
    }

    if(plots) {
        bp_plot <- make_piecewise_bp_plot(plot_df,
                                          .x,
                                          .y,
                                          estimate_res$lm_left,
                                          estimate_res$lm_right,
                                          estimate_res$bp_dat)
    } else {
        bp_plot <- NULL
    }

    return(list(breakpoint_data = estimate_res$bp_dat,
                fitted_vals = estimate_res$pred,
                lm_left = estimate_res$lm_left,
                lm_right = estimate_res$lm_right,
                lm_simple = estimate_res$lm_simple,
                bp_plot = bp_plot))
}

#' @keywords internal
loop_v_slope <- function(.data,
                         .x,
                         .y,
                         pos_change,
                         pos_slope_after_bp,
                         alpha_linearity = 0.05,
                         conf_level = 0.95) {
    # you can't fit a model if you don't have any data
    if(nrow(.data) == 0) return(NULL)

    n_rows <- nrow(.data)

    # initialize empty vectors
    RSS_two <- MSE_two <- f_stat <- pf_two <- pct_slope_change <-
        int_point_x <- dist_MSE_ratio <- numeric(n_rows)
    pos_change_vec <- pos_slope_after_bp_vec <- logical(n_rows)

    lm_simple <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::lm(data = .data)

    RSS_simple <- sum(stats::resid(lm_simple)^2)
    # recip_slope <- (-1 / lm_simple$coefficients[2]) # used in for loop
    # find slope of line perpendicular to slope of lm_simple
    # There's an easier way to find the distance between a point and a line
    # by putting the equation of the simple regression into the standard form
    # of Ax + By + C = 0
    # https://www.mathportal.org/calculators/analytic-geometry/line-point-distance.php
    # https://brilliant.org/wiki/dot-product-distance-between-point-and-a-line/

    for(i in 1:n_rows) {
        # setting first and last iteration to NA keeps the index lengths equal
        if(i == 1 | i == n_rows) {
            RSS_two[i] <- MSE_two[i] <- f_stat[i] <- pf_two[i] <-
                pos_change_vec[i] <- pos_slope_after_bp_vec[i] <-
                int_point_x[i] <- dist_MSE_ratio[i] <- NA
            next
        }

        # split data into left and right portions, sharing division point
        df_left <- .data[1:i,]
        df_right <- .data[i:n_rows,]

        # make linear models of the two regressions
        lm_left <- paste0(.y, " ~ ", "1 + ", .x) %>%
            stats::lm(data = df_left)
        lm_right <- paste0(.y, " ~ ", "1 + ", .x) %>%
            stats::lm(data = df_right)

        RSS_two[i] <- sum(stats::resid(lm_left)^2) +
            sum(stats::resid(lm_right)^2)
        # divide by -4 b/c estimating 4 parameters
        MSE_two[i] <- RSS_two[i] / (nrow(lm_simple$model) - 4)
        f_stat[i] <- (RSS_simple - RSS_two[i]) / (2 * MSE_two[i])
        pf_two[i] <- stats::pf(
            f_stat[i], df1 = 2, df2 = nrow(lm_simple$model) - 4,
                               lower.tail = FALSE)

        pct_slope_change[i] <- 100*(lm_right$coefficients[2] -
                                        lm_left$coefficients[2]) /
            abs(lm_left$coefficients[2])
        pos_change_vec[i] <- dplyr::if_else(pct_slope_change[i] > 0, TRUE, FALSE)
        pos_slope_after_bp_vec[i] <- dplyr::if_else(
            lm_right$coefficients[2] > 0, TRUE, FALSE)

        # find intersection point of left and right regressions
        lr_intersect <- intersection_point(lm_left, lm_right)
        int_point_x[i] <- lr_intersect["x"]

        # solve for distance using standard form of line Ax + By + C = 0
        d <- abs(lm_simple$coefficients[2] * lr_intersect["x"] -
                lr_intersect["y"] + lm_simple$coefficients[1]) /
            sqrt(lm_simple$coefficients[2]^2 + (-1)^2)

        dist_MSE_ratio[i] <- d / MSE_two[i]
    }
    # calculate cutoff for finding approximate confidence interval
    crit_F <- stats::qf(conf_level, 1, n_rows - 4, lower.tail = TRUE)
    # using min(MSE_two) and min(RSS_two) based on JM paper
    inside_ci <- dplyr::if_else(
        (RSS_two - min(RSS_two, na.rm = TRUE)) /
            min(MSE_two, na.rm =TRUE) < crit_F,
        TRUE, FALSE)

    # for debugging purposes, plot breakpoints inside 95% CI
    # plot((RSS_two - min(RSS_two, na.rm = TRUE)) / MSE_two)
    # abline(h = crit_F)
    # add x_range constraint
    loop_stats <- tibble::tibble(p = pf_two,
                                pos_change = pos_change_vec,
                                pos_slope_after_bp = pos_slope_after_bp_vec,
                                dist_MSE_ratio = dist_MSE_ratio,
                                int_point_x = int_point_x,
                                inside_ci = inside_ci
    ) %>%
        dplyr::mutate(idx = dplyr::row_number())

    loop_stats
}

#' @keywords internal
num_to_na <- function(x) {
    if(is.numeric(x)) {
        x <- NA
    }
    x
}

#' @keywords internal
get_v_slope_res <- function(.data, bp_idx, .x, .y, bp,
                            alpha_linearity, pos_change,
                            pos_slope_after_bp,
                            est_ci = c("estimate", "lower_ci", "upper_ci"),
                            ...) {

    df_left <- .data[1:bp_idx,] # split data into left portion
    df_right <- .data[bp_idx:nrow(.data),] # split data into right portion

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

    int_point <- intersection_point(lm_left, lm_right)

    determinant_bp <- and(check_if_determinant_bp(
        p = pf_two,
        pct_slope_change = pct_slope_change,
        pos_change = pos_change,
        pos_slope_after_bp = pos_slope_after_bp,
        slope_after_bp = stats::coef(lm_right)[2],
        alpha = alpha_linearity),
        # the intersection point needs to be in the range of x
        dplyr::between(int_point["x"],
                       range(.data[[.x]])[1],
                       range(.data[[.x]])[2]))

    y_hat_left <- tibble::tibble("{.x}" := df_left[[.x]],
                                 "{.y}" := lm_left$fitted.values,
                                 algorithm = "v-slope")
    y_hat_right <- tibble::tibble("{.x}" := df_right[[.x]],
                                  "{.y}" := lm_right$fitted.values,
                                  algorithm = "v-slope")
    pred <- dplyr::bind_rows(y_hat_left, y_hat_right)

    # find closest data point to intersection point and prepare output
    bp_dat <- find_threshold_vals(.data = .data, thr_x = int_point["x"],
                                  thr_y = int_point["y"], .x = .x,
                                  .y = .y, ...)
    bp_dat <- bp_dat %>%
        dplyr::mutate(algorithm = "v-slope",
                      est_ci = est_ci,
                      x_var = .x,
                      y_var = .y,
                      determinant_bp = determinant_bp,
                      bp = bp,
                      pct_slope_change = pct_slope_change,
                      f_stat = f_stat,
                      p_val_f = pf_two) %>%
        dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp, est_ci,
                        pct_slope_change, f_stat, p_val_f)
    list(bp_dat = bp_dat,
         fitted_vals = pred,
         lm_left = lm_left,
         lm_right = lm_right,
         lm_simple = lm_simple)
}
