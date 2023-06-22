#' Finding a breakpoint using Orr's 'bruteforce' algorithm.
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
#' @param front_trim_vt1 How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param front_trim_vt2 How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param pos_slope_after_bp Should the slope after the breakpoint be positive? Default is `TRUE`. This catches cases when the percent change in slope is positive, but the second slope is still negative. Change to `FALSE` when PetCO2 is the y-axis variable.
#' @param ci Should the output include confidence interval data? Default is `FALSE`.
#' @param conf_level Confidence level to use if calculating confidence intervals.
#' @param plots Should this function generate plots? Set to `FALSE` to save time.
#'
#' @returns A list including slice of the original data frame at the threshold index with new columns `algorithm`, `determinant_bp`, `pct_slope_change`, `f_stat`, and `p_val_f.` The list also includes the fitted values, the left and right sides of the piecewise regression, and a simple linear regression.
#'
#' @details
#' According to the Orr et al. (1982, p. 1350) paper, the "anaerobic threshold is then reported as the first intersection point [of the lines] of the more appropriate model." In an email to Dr. Hughson, a weakness of this method is that the intersection point hardley ever lands on an existing data point. In some cases, the intersection of the best-fit model lies beyond the range of the x-axis. Instead of returning that value, we report the solution with the lowest pooled sums of squares that satisfies the following criteria: 1) the change in slope from the left to the right regression matches the anticipated change, usually positive; 2) the slope of the right regression matches the anticipated slope, usually positive; 3) the p-value of the F-test based on the extra sums of squares principle is statistically significant; and 4) the x-coordinate of the intersection point is within the range of the observed x-values. This differs slightly from the original method, but provides potentially plausible answers when the alternative is a negative VO2 or a VO2 well beyond VO2max.
#'
#' The text of the original paper states that "Regression lines are calculated for all possible divisions of the data into two \emph{contiguous} groups." We interpret contiguous to mean "shares a boundary". That is, the left and right portions of the data each contain one identical data point at the division.
#'
#' @importFrom rlang :=
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
                bp,
                ...,
                vo2 = "vo2",
                vco2 = "vco2",
                ve = "ve",
                time = "time",
                ordering = c("by_x", "time"),
                alpha_linearity = 0.05,
                front_trim_vt1 = 60,
                front_trim_vt2 = 60,
                pos_change = TRUE,
                pos_slope_after_bp = TRUE,
                ci = FALSE,
                conf_level = 0.95,
                plots = TRUE
                ) {

    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
    bp <- match.arg(bp, choices = c("vt1", "vt2"), several.ok = FALSE)
    ordering <- match.arg(ordering, several.ok = FALSE)

    .data <- order_cpet_df(.data, .x = .x , time = time,
                           ordering = ordering)
    plot_df <- .data

    front_trim <- set_front_trim(bp = bp,
                                 front_trim_vt1 = front_trim_vt1,
                                 front_trim_vt2 = front_trim_vt2)
    .data <- .data %>%
        dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))

    loop_res <- loop_orr(.data = .data,
                         .x = .x,
                         .y = .y,
                         alpha_linearity = alpha_linearity,
                         conf_level = conf_level)

    if(!is.null(loop_res)) {
        best_idx <- get_best_piecewise_idx(
            loop_res,
            range(.data[[.x]]),
            alpha_linearity = alpha_linearity,
            pos_change = pos_change,
            pos_slope_after_bp = pos_slope_after_bp)
    } else {
        best_idx <- numeric()
    }

    # return quick summary if generating models fails
    if(is.null(loop_res) | length(best_idx) == 0) {
        bp_dat <- return_indeterminant_findings(
            .data = plot_df,
            bp = bp,
            algorithm = as.character(match.call()[[1]]),
            .x = .x,
            .y = .y,
            est_ci = "estimate")

        return(list(breakpoint_data = bp_dat))
    }

    estimate_res <- get_orr_res(.data = .data,
                           bp_idx = best_idx,
                           .x = .x,
                           .y = .y,
                           bp = bp,
                           alpha_linearity = alpha_linearity,
                           pos_change = pos_change,
                           pos_slope_after_bp = pos_slope_after_bp,
                           est_ci = "estimate")

    if(ci) {
        ci_lower_idx <- loop_res %>%
            dplyr::filter(inside_ci) %>%
            dplyr::filter(int_point_x == min(int_point_x)) %>%
            dplyr::select(idx) %>%
            dplyr::pull()

        lower_ci_res <- get_orr_res(.data = .data,
                                    bp_idx = ci_lower_idx,
                                    .x = .x,
                                    .y = .y,
                                    bp = bp,
                                    alpha_linearity = alpha_linearity,
                                    pos_change = pos_change,
                                    pos_slope_after_bp = pos_slope_after_bp,
                                    est_ci = "lower")

        ci_upper_idx <- loop_res %>%
            dplyr::filter(inside_ci) %>%
            dplyr::filter(int_point_x == max(int_point_x)) %>%
            dplyr::select(idx) %>%
            dplyr::pull()

        upper_ci_res <- get_orr_res(.data = .data,
                                        bp_idx = ci_upper_idx,
                                        .x = .x,
                                        .y = .y,
                                        bp = bp,
                                        alpha_linearity = alpha_linearity,
                                        pos_change = pos_change,
                                        pos_slope_after_bp = pos_slope_after_bp,
                                        est_ci = "upper")

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
                lm_right = estimate_res$m_right,
                lm_simple = estimate_res$lm_simple,
                bp_plot = bp_plot))

    # Should we add the three line regression code later?
}

#' @keywords internal
#' This function loops through all possible divisions of the data and finds the
#' best solution
loop_orr <- function(.data,
                     .x,
                     .y,
                     alpha_linearity = 0.05,
                     conf_level = 0.95) {
    # you can't fit a model if you don't have any data
    if(nrow(.data) == 0) return(NULL)

    n_rows <- nrow(.data)

    # initialize empty vectors
    RSS_two <- MSE_two <- f_stat <- pf_two <-
        pct_slope_change <- int_point_x <- numeric(n_rows)
    pos_change <- pos_slope_after_bp <- logical(n_rows)

    # create simple linear model and calculate its RSS
    lm_simple <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::lm(data = .data)
    RSS_simple <- sum(stats::resid(lm_simple)^2)

    for(i in 1:n_rows) {
        # setting first and last iteration to NA keeps the index lengths equal
        if(i == 1 | i == n_rows) {
            RSS_two[i] <- MSE_two[i] <- f_stat[i] <- pf_two[i] <-
                pos_change[i] <- pos_slope_after_bp[i] <-
                int_point_x[i] <- NA
            next
        }
        # split data into left and right halves
        # should these share the same point? Or should, they be different by 1 point?
        # Ekkekakis said WB makes them differ by 1 pt
        # I wonder exactly what the word 'contiguous' means in this context
        df_left <- .data[1:i,]
        df_right <- .data[i:n_rows,]
        # NOTE: WB'S implementation of Orr does NOT share a point at i

        # make linear models of the two halves
        lm_left <- paste0(.y, " ~ ", "1 + ", .x) %>%
            stats::lm(data = df_left)
        lm_right <- paste0(.y, " ~ ", "1 + ", .x) %>%
            stats::lm(data = df_right)

        # ss_left[i] <- sum((lm_left$residuals)^2)
        # ss_right[i] <- sum((lm_right$residuals)^2)

        RSS_two[i] <- sum(stats::resid(lm_left)^2) +
            sum(stats::resid(lm_right)^2)
        MSE_two[i] <- RSS_two[i] / (nrow(lm_simple$model) - 4) # -4 b/c estimating 4 parameters
        f_stat[i] <- (RSS_simple - RSS_two[i]) / (2 * MSE_two[i])
        pf_two[i] <- stats::pf(f_stat[i], df1 = 2,
                               df2 = nrow(lm_simple$model) - 4,
                            lower.tail = FALSE)

        pct_slope_change[i] <- 100*(lm_right$coefficients[2] -
                                        lm_left$coefficients[2]) /
            abs(lm_left$coefficients[2])
        pos_change[i] <- dplyr::if_else(pct_slope_change[i] > 0, TRUE, FALSE)
        pos_slope_after_bp[i] <- dplyr::if_else(lm_right$coefficients[2] > 0,
                                                TRUE, FALSE)

        int_point_x[i] <- intersection_point(lm_left, lm_right)["x"]
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

    loop_stats <- tibble::tibble(p = pf_two,
                        pos_change = pos_change,
                        pos_slope_after_bp = pos_slope_after_bp,
                        int_point_x = int_point_x,
                        inside_ci = inside_ci
                        ) %>%
        dplyr::mutate(idx = dplyr::row_number())

    loop_stats
}

#' @keywords internal
#' Similar to loop_orr(), but gets results at specific index
get_orr_res <- function(.data, bp_idx, .x, .y, bp,
                        alpha_linearity, pos_change,
                        pos_slope_after_bp,
                        est_ci = c("estimate", "lower_ci", "upper_ci"),
                        ...) {
    est_ci <- match.arg(est_ci, several.ok = FALSE)

    df_left <- .data[1:bp_idx,]
    df_right <- .data[bp_idx:nrow(.data),]

    lm_left <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::lm(data = df_left)
    lm_right <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::lm(data = df_right)
    lm_simple <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::lm(data = .data)

    pw_stats <- piecewise_stats(lm_left, lm_right, lm_simple)
    list2env(pw_stats, envir = environment())

    pct_slope_change <- 100*(lm_right$coefficients[2] -
                                 lm_left$coefficients[2]) /
        abs(lm_left$coefficients[2])

    y_hat_left <- tibble::tibble("{.x}" := df_left[[.x]],
                                 "{.y}" := lm_left$fitted.values,
                                 algorithm = "orr")
    y_hat_right <- tibble::tibble("{.x}" := df_right[[.x]],
                                  "{.y}" := lm_right$fitted.values,
                                  algorithm = "orr")
    pred <- dplyr::bind_rows(y_hat_left, y_hat_right)

    int_point <- intersection_point(lm_left, lm_right)

    determinant_bp <- and(check_if_determinant_bp(
        p = pf_two,
        pct_slope_change = pct_slope_change,
        pos_change = pos_change,
        pos_slope_after_bp = pos_slope_after_bp,
        slope_after_bp = stats::coef(lm_right)[2],
        alpha = alpha_linearity),
        dplyr::between(int_point["x"],
                       range(.data[[.x]])[1],
                       range(.data[[.x]])[2]))

    bp_dat <- find_threshold_vals(.data = .data,
                                  thr_x = int_point["x"],
                                  thr_y = int_point["y"],
                                  .x = .x,
                                  .y = .y,
                                  ...)

    bp_dat <- bp_dat %>%
        dplyr::mutate(algorithm = "orr",
                      est_ci = est_ci,
                      x_var = .x,
                      y_var = .y,
                      determinant_bp = determinant_bp,
                      bp = bp,
                      pct_slope_change = pct_slope_change,
                      f_stat = f_stat,
                      p_val_f = pf_two) %>%
        dplyr::relocate(bp, algorithm, determinant_bp, est_ci,
                        pct_slope_change, f_stat, p_val_f)

    list(bp_dat = bp_dat,
         fitted_vals = pred,
         lm_left = lm_left,
         lm_right = lm_right,
         lm_simple = lm_simple)
}
