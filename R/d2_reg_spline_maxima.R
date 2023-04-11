#' Find a breakpoint using the maxima in the 2nd derivative of a regression spline
#'
#' #' Use regression splines of with an increasing number of knots to find the best-fit line for the data, take the second derivative, and find the highest maxima of the second derivative. This method is traditionally used with ventilatory equivalents (VE/VO2 or VE/VCO2), PetO2, or PetCO2 vs. time or vs. VO2.
#'
#' @details
#' #' Several works use an algorithm similar to this one. However, it is somewhat unclear if those previous works were truly using "smoothing splines" or "regression splines" (they are *not* the same). Previous works may refer to a third or fifth order "polynomial smoothing splines", which likely means "regression splines", but we cannot be certain at this time.
#'
#' Unlike several published works, this method iteratively finds the best-fit regression by increasing the number of knots and using the likelihood ratio test with the \code{anova} function. When the likelihood ratio test does *not* find a statistically significant difference between the previous and the newest regression spline, this returns the previous regression spline. We found that this is a more conservative approach and yields fewer possible second derivative maxima to choose between. However, users can specify their desired number of knots.
#'
#' @param .data Gas exchange data frame or tibble.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param bp Is this the first (\code{vt1}) or the second (\code{vt2}) breakpoint?
#' @param degree The polynomial degree the function should fit to the curve. Most previous studies using this method use 3 or 5.
#' @param df If the user wants to use a specific number of knots, your `df` should equal the number of desired knots +`degree` to account for the df "consumed" by the polynomial regressions.
#' @param vo2 Name of the \code{vo2} variable
#' @param vco2 Name of the \code{vco2} variable
#' @param ve Name of the \code{ve} variable
#' @param time Name of the \code{time} variable
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.
#' @param pos_change Do you expect the slope to be increasing at the breakpoint? This helps with filtering maxima.
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param front_trim_vt1 How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param front_trim_vt2 How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param ci Should the output include confidence interval data? Default is `FALSE`.
#' @param conf_level Confidence level to use if calculating confidence intervals.
#' @param plots Should this function generate plots? Set to `FALSE` to save time.
#'
#' #' @references
#' Leo, J. A., Sabapathy, S., Simmonds, M. J., & Cross, T. J. (2017). The Respiratory Compensation Point is Not a Valid Surrogate for Critical Power. Medicine and science in sports and exercise, 49(7), 1452-1460.
#' Sherrill, D. L., Anderson, S. J., & Swanson, G. E. O. R. G. E. (1990). Using smoothing splines for detecting ventilatory thresholds. Medicine and Science in Sports and Exercise, 22(5), 684-689.
#' Wade, T. D., Anderson, S. J., Bondy, J., Ramadevi, V. A., Jones, R. H., & Swanson, G. D. (1988). Using smoothing splines to make inferences about the shape of gas-exchange curves. Computers and biomedical research, 21(1), 16-26.
#'
#' @returns A slice of the original data frame at the threshold index with a new `algorithm` column.
#' @export
#'
#' @examples
#' # TODO write an example
#'
d2_reg_spline_maxima <- function(.data,
                                 .x,
                                 .y,
                                 bp,
                                 degree = 5,
                                 df = NULL,
                                 vo2 = "vo2",
                                 vco2 = "vco2",
                                 ve = "ve",
                                 time = "time",
                                 alpha_linearity = 0.05,
                                 pos_change = TRUE,
                                 ordering = c("by_x", "time"),
                                 front_trim_vt1 = 60,
                                 front_trim_vt2 = 60,
                                 ci = FALSE,
                                 conf_level = 0.95,
                                 plots = TRUE,
                                 ...
                                 ) {

    # check if there is crucial missing data
    stopifnot(!any(missing(.data),
                   missing(.x),
                   missing(.y),
                   missing(bp)))

    ordering <- match.arg(ordering, several.ok = FALSE)
    .data <- order_cpet_df(.data, .x = .x , time = time,
                           ordering = ordering)

    plot_df <- .data

    front_trim <- set_front_trim(bp = bp,
                                 front_trim_vt1 = front_trim_vt1,
                                 front_trim_vt2 = front_trim_vt2)
    .data <- .data %>%
        dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))

    # find best number of knots
    lm_spline <- loop_d2_reg_spline(.data = .data, .x = .x, .y = .y,
                                    df = df, degree = degree)

    # It feels like there's a better way than splinefun to find a more exact
    # 2nd derivative, but I don't know it

    # find maxima in 2nd derivative. This is accomplished by making
    # a spline interpolation function. Use this with optimize() function.

    # get new values at equal spacing for a smoother splinefun result
    equi_spaced_x <- seq(from = min(.data[[.x]]), to = max(.data[[.x]]),
                         length.out = range(.data[[vo2]]) %>%
                             diff() %>%
                             round())
    pred <- stats::predict(
        lm_spline,
        newdata = tibble::tibble("{.x}" := equi_spaced_x))

    spline_func <- stats::splinefun(x = equi_spaced_x, y = pred)

    # find number of maxima
    if(pos_change) {
        sign_change_idx <- slope_sign_changes(y = spline_func(x = equi_spaced_x,
                                                              deriv = 2),
                                              change = "pos_to_neg")
        # filter by expected slope (usually positive)
        # there can be a spike in accel during the initial drop in vent eqs
        # however, only remove these thresholds with a negative derivative
        # if there is more than one value with a negative slope
        if(length(sign_change_idx) > 1) {
            # see which derivatives at the sign changes are negative
            neg_slope_idx <- which(spline_func(
                x = equi_spaced_x[sign_change_idx],
                deriv = 1) < 0)
            while(length(neg_slope_idx) > 0 & length(sign_change_idx) > 1) {
                # remove one value with a negative slope at a time in case
                # all are negative. If the derivative is negative but it's
                # right befor a clear upturn, that's probably the threshold
                sign_change_idx <- sign_change_idx[-neg_slope_idx[1]]
                # find which values have negative signs to exit loop
                neg_slope_idx <- which(spline_func(
                    x = equi_spaced_x[sign_change_idx],
                    deriv = 1) < 0)
            }
        }
    } else { # this only applies when y var = petco2
        sign_change_idx <- slope_sign_changes(y = spline_func(x = equi_spaced_x,
                                                              deriv = 2),
                                              change = "neg_to_pos")
        # filter by expected slope (in this case negative)
        if(length(sign_change_idx) > 1) {
            # see which derivatives at the sign changes are positive
            pos_slope_idx <- which(spline_func(
                x = equi_spaced_x[sign_change_idx],
                deriv = 1) > 0)
            while(length(pos_slope_idx) > 0 & length(sign_change_idx) > 1) {
                # remove values with a positive slope by logical indexing
                # should I still be removing these incrementally? in case the
                # nadir isocapnic buffering period is VERY short-lived?
                # I could see how the deriv could be negative, but the
                # accel could be very positive. this would be okay if
                # it were the only maxima (I think)
                sign_change_idx <- sign_change_idx[-pos_slope_idx[1]]
                # find values with positive slope again to exit loop
                pos_slope_idx <- which(spline_func(
                    x = equi_spaced_x[sign_change_idx],
                    deriv = 1) > 0)
            }
        }
    }

    # there can be multiple local extrema. Choose the extrema closer to the
    # nadir or peak because this often represents the beginning of a systematic
    # rise

    # choosing closest to the nadir/peak probably isn't a good idea when
    # the graph is NOT a ventilatory equivalents graph b/c the nadir will
    # be super close to the beginning of the test
    # At least for now, for better or for worse, I'd say choose the extrema
    # with the largest y-magnitude

    y_val_sign_change <- spline_func(x = equi_spaced_x[sign_change_idx],
                                     deriv = 1)

    if(length(sign_change_idx) > 1) {
        if(pos_change) {

            sign_change_idx <- sign_change_idx[which.max(
                spline_func(x = equi_spaced_x[sign_change_idx], deriv = 2)
            )]
            # OLD code
            # find nadir of y values
            # nadir <- min(pred)
            # sign_change_idx <- sign_change_idx[which.min(abs(
            #     nadir - equi_spaced_x[sign_change_idx]))]
        } else {
            sign_change_idx <- sign_change_idx[which.min(
                spline_func(x = equi_spaced_x[sign_change_idx], deriv = 2)
            )]
            # OLD code
            # find peak (when using petco2 basically)
            # peak <- max(pred)
            # sign_change_idx <- sign_change_idx[which.min(abs(
            #     peak - equi_spaced_x[sign_change_idx]))]
        }
    }

    # OLD CODE, KEEP FOR A BIT
    # there are often two local maxima. Choose the higher x-axis value
    # for VT2/RC, and the lower for VT1. This is similar to
    # Sherril et al. (1990) and Cross et al. (2012)
    # if(length(sign_change_idx) > 1) {
    #     if(bp == "vt1") {
    #         sign_change_idx <- sign_change_idx[which.min(sign_change_idx)]
    #     }
    #     if(bp == "vt2") {
    #         sign_change_idx <- sign_change_idx[which.max(sign_change_idx)]
    #     }
    # }

    # calculate y-value at the threshold
    # get value on either side of index
    # interval <- c(equi_spaced_x[sign_change_idx[i] - 1],
    #               equi_spaced_x[sign_change_idx[i] + 1])
    #
    # extrema_deriv2 <- stats::optimize(f = spline_func, interval = interval,
    #                                   deriv = 2, maximum = TRUE)
    #
    # y_val_sign_changes <- extrema_deriv2$objective

    if(length(sign_change_idx) > 0) {
        x_threshold <- equi_spaced_x[sign_change_idx]
        y_hat_threshold <- stats::predict(lm_spline,
                                   tibble::tibble("{.x}" := x_threshold))
        # get values at threshold
        bp_dat <- find_threshold_vals(.data = .data, thr_x = x_threshold,
                                      thr_y = y_hat_threshold, .x = .x,
                                      .y = .y, ...)

    } else {
        bp_dat <- tibble::tibble()
    }

    if(nrow(bp_dat) == 0) { # no breakpoint found
        bp_dat <- bp_dat %>%
            dplyr::add_row() %>%
            dplyr::mutate(determinant_bp = FALSE)
        # add character or factor columns that are all the same value (e.g. ids)
        non_numeric_df <- .data %>%
            dplyr::select(tidyselect::where(
                function(x) is.character(x) | is.factor(x) &
                             all(x == x[1]))) %>%
            dplyr::slice(1)
        bp_dat <- dplyr::bind_cols(bp_dat, non_numeric_df)
    } else { # breakpoint found
        bp_dat <- bp_dat %>%
            dplyr::mutate(determinant_bp = TRUE)
    }

    bp_dat <- bp_dat %>%
        dplyr::mutate(algorithm = "d2_reg_spline_maxima",
                      x_var = .x,
                      y_var = .y,
                      bp = bp,
                      est_ci = "estimate")

    if(ci) {
        # nonparametric bootstrapping
        boot_res <- boot::boot(data = .data,
                               statistic = bootstrap_ci_d2_reg_spline,
                               R = 100,
                               sim = "ordinary",
                               .x = .x,
                               .y = .y,
                               df = df,
                               degree = degree,
                               pos_change = pos_change,
                               parallel = "multicore")
        # calculate percentile CI's
        boot_ci <- broom::tidy(boot_res,
                               conf.int = TRUE,
                               conf.method = "perc")
        # calculate new values at threshold
        ci_lower_x <- boot_ci$conf.low
        ci_lower_y <- stats::predict(lm_spline,
                       tibble::tibble("{.x}" := ci_lower_x))

        lower_ci_res <- find_threshold_vals(.data = .data,
                                            thr_x = ci_lower_x,
                                            thr_y = ci_lower_y,
                                            .x = .x,
                                            .y = .y,
                                            ...) %>%
            dplyr::mutate(bp = bp,
                          algorithm = "d2_reg_spline_maxima",
                          x_var = .x,
                          y_var = .y,
                          est_ci = "lower_ci"
            )

        ci_upper_x <- boot_ci$conf.high
        ci_upper_y <- stats::predict(lm_spline,
                                     tibble::tibble("{.x}" := ci_upper_x))

        upper_ci_res <- find_threshold_vals(.data = .data,
                                            thr_x = ci_upper_x,
                                            thr_y = ci_upper_y,
                                            .x = .x,
                                            .y = .y,
                                            ...) %>%
            dplyr::mutate(bp = bp,
                          algorithm = "d2_reg_spline_maxima",
                          x_var = .x,
                          y_var = .y,
                          est_ci = "upper_ci"
            )

        bp_dat <- dplyr::bind_rows(lower_ci_res,
                                   bp_dat,
                                   upper_ci_res)
        # if the estimate was true, set the others to TRUE
        if(any(bp_dat[["determinant_bp"]], na.rm = TRUE)) {
            bp_dat[["determinant_bp"]] <- TRUE
        }
    }

    bp_dat <- bp_dat %>%
        dplyr::relocate(bp, algorithm, x_var, y_var, est_ci, determinant_bp)

    if(plots) {
        bp_plot <- ggplot2::ggplot(data = plot_df,
                                   ggplot2::aes(x = .data[[.x]],
                                                y = .data[[.y]])) +
            ggplot2::geom_point(alpha = 0.5) +
            ggplot2::geom_line(
                data = tibble::tibble(x = equi_spaced_x,
                                      y = pred),
                ggplot2::aes(x = equi_spaced_x, y = pred)) +
            ggplot2::geom_vline(xintercept = bp_dat[[.x]]) +
            ggplot2::theme_minimal()
    } else {
        bp_plot <- NULL
    }

    return(list(breakpoint_data = bp_dat,
                lm_reg_spline = lm_spline,
                spline_deriv_func = spline_func,
                bp_plot = bp_plot))
}

#' @keywords internal
loop_d2_reg_spline <- function(.data, .x, .y, df = NULL,
                               degree = 5, alpha_linearity = 0.05) {
    # TODO allow users to specify b-spline or natural-spline basis
    # would that use do.call()?

    # if statement for if users specify the knots or df
    if(!is.null(df)) {
        lm_spline <- paste0(.y, " ~ ", "1 + ",
                            "splines::bs(", .x,
                            ", df = ", df,
                            ", degree = ", degree, ")") %>%
            stats::lm(data = .data)
    } else {
        spline_mod_list = vector(mode = "list", length = 0)
        cont <- TRUE
        i <- 1
        # reference model with one interior knot
        lm_spline <- paste0(.y, " ~ ", "1 + ",
                            "splines::bs(", .x,
                            ", df = ", i + degree,
                            ", degree = ", degree, ")") %>%
            stats::lm(data = .data)
        spline_mod_list <- append(spline_mod_list, list(lm_spline))
        # while loop beings with 1 knot (3 df already used assuming a 3rd spline)
        while(cont == TRUE) {
            # browser()
            i <- i + 1
            # TODO add options for user-defined knots and df
            lm_spline <- paste0(.y, " ~ ", "1 + ",
                                "splines::bs(", .x,
                                ", df = ", i + degree,
                                ", degree = ", degree, ")") %>%
                stats::lm(data = .data)
            spline_mod_list <- append(spline_mod_list, list(lm_spline))
            lrt <- stats::anova(spline_mod_list[[i-1]], spline_mod_list[[i]])
            if (is.na(lrt$`Pr(>F)`[2]) | lrt$`Pr(>F)`[2] >= alpha_linearity) {
                cont = FALSE
                lm_spline <- spline_mod_list[[i-1]] # take the previous model
            }
        }
    }
    lm_spline
}

#' @keywords internal
bootstrap_ci_d2_reg_spline <- function(.data,
                                       i,
                                       .x,
                                       .y,
                                       df,
                                       degree,
                                       pos_change) {

    # find best number of knots
    lm_spline <- loop_d2_reg_spline(.data = .data[i,], .x = .x, .y = .y,
                                    df = df, degree = degree)

    # It feels like there's a better way than splinefun to find a more exact
    # 2nd derivative, but I don't know it

    # find maxima in 2nd derivative. This is accomplished by making
    # a spline interpolation function. Use this with optimize() function.

    # get new values at equal spacing for a smoother splinefun result
    equi_spaced_x <- seq(from = min(.data[i,][[.x]]),
                         to = max(.data[i,][[.x]]),
                         length.out = range(.data[i,][[.x]]) %>%
                             diff() %>%
                             round())
    pred <- stats::predict(lm_spline,
                           newdata = tibble::tibble("{.x}" := equi_spaced_x))

    spline_func <- stats::splinefun(x = equi_spaced_x, y = pred)

    # find number of maxima
    if(pos_change) {
        sign_change_idx <- slope_sign_changes(y = spline_func(x = equi_spaced_x,
                                                              deriv = 2),
                                              change = "pos_to_neg")
        # filter by expected slope (usually positive)
        # there can be a spike in accel during the initial drop in vent eqs
        # however, only remove these thresholds with a negative derivative
        # if there is more than one value with a negative slope
        if(length(sign_change_idx) > 1) {
            # see which derivatives at the sign changes are negative
            neg_slope_idx <- which(spline_func(
                x = equi_spaced_x[sign_change_idx],
                deriv = 1) < 0)
            while(length(neg_slope_idx) > 0 & length(sign_change_idx) > 1) {
                # remove one value with a negative slope at a time in case
                # all are negative. If the derivative is negative but it's
                # right befor a clear upturn, that's probably the threshold
                sign_change_idx <- sign_change_idx[-neg_slope_idx[1]]
                # find which values have negative signs to exit loop
                neg_slope_idx <- which(spline_func(
                    x = equi_spaced_x[sign_change_idx],
                    deriv = 1) < 0)
            }
        }
    } else { # this only applies when y var = petco2
        sign_change_idx <- slope_sign_changes(y = spline_func(x = equi_spaced_x,
                                                              deriv = 2),
                                              change = "neg_to_pos")
        # filter by expected slope (in this case negative)
        if(length(sign_change_idx) > 1) {
            # see which derivatives at the sign changes are positive
            pos_slope_idx <- which(spline_func(
                x = equi_spaced_x[sign_change_idx],
                deriv = 1) > 0)
            while(length(pos_slope_idx) > 0 & length(sign_change_idx) > 1) {
                # remove values with a positive slope by logical indexing
                # should I still be removing these incrementally? in case the
                # nadir isocapnic buffering period is VERY short-lived?
                # I could see how the deriv could be negative, but the
                # accel could be very positive. this would be okay if
                # it were the only maxima (I think)
                sign_change_idx <- sign_change_idx[-pos_slope_idx[1]]
                # find values with positive slope again to exit loop
                pos_slope_idx <- which(spline_func(
                    x = equi_spaced_x[sign_change_idx],
                    deriv = 1) > 0)
            }
        }
    }

    # there can be multiple local extrema. Choose the extrema closer to the
    # nadir or peak because this often represents the beginning of a systematic
    # rise

    # choosing closest to the nadir/peak probably isn't a good idea when
    # the graph is NOT a ventilatory equivalents graph b/c the nadir will
    # be super close to the beginning of the test
    # At least for now, for better or for worse, I'd say choose the extrema
    # with the largest y-magnitude

    y_val_sign_change <- spline_func(x = equi_spaced_x[sign_change_idx],
                                     deriv = 1)

    if(length(sign_change_idx) > 1) {
        if(pos_change) {

            sign_change_idx <- sign_change_idx[which.max(
                spline_func(x = equi_spaced_x[sign_change_idx], deriv = 2)
            )]
            # OLD code
            # find nadir of y values
            # nadir <- min(pred)
            # sign_change_idx <- sign_change_idx[which.min(abs(
            #     nadir - equi_spaced_x[sign_change_idx]))]
        } else {
            sign_change_idx <- sign_change_idx[which.min(
                spline_func(x = equi_spaced_x[sign_change_idx], deriv = 2)
            )]
            # OLD code
            # find peak (when using petco2 basically)
            # peak <- max(pred)
            # sign_change_idx <- sign_change_idx[which.min(abs(
            #     peak - equi_spaced_x[sign_change_idx]))]
        }
    }

    if(length(sign_change_idx) > 0) {
        return(equi_spaced_x[sign_change_idx])
    } else {
        return(NA)
    }

}
