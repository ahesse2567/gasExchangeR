#' Find a breakpoint using the maxima in the 2nd derivative of a polynomial regression
#'
#' Use polynomial regression of increasing order to find the best-fit line for the data, take the second derivative, and find the highest maxima of the second derivative. This method is traditionally used with ventilatory equivalents (VE/VO2 or VE/VCO2), PetO2, or PetCO2 vs. time or vs. VO2.
#'
#' Note, polynomial regressions are usually inferior to regression splines or to smoothing splines because polynomial regressions because polynomial regressions can overfit the data more easily. However, we offer this method in case users want to reproduce previous works and while we develop this package.
#'
#' Unlike several published works, this method iteratively finds the best-fit polynomial regression by increasing the polynomial order and using the likelihood ratio test with the \code{anova} function. When the likelihood ratio test does *not* find a statistically significant difference between the previous and the newest polynomial, the function uses the previous polynomial. However, users can specify the polynomial \code{degree} if desired.
#'
#' #' @details
#' Our implementation is similar to two studies by Wisén and Wohlfart et al. in that this uses a polynomial regression. However, they only found the \emph{first} derivative of the VO2 and VCO2 vs. time curves and denote the threshold as their crossing point. Using ventilatory equivalents vs. time, they considered the threshold as when the \emph{first} derivative of each ventilatory equivalent increased above zero (for the last time). In contrast, the majority of other literature we are aware of that uses regression splines or smoothing splines finds the \emph{second} derivative. We found that this is a more conservative approach and yields fewer possible second derivative maxima to choose between. This function finds the second derivative by default to match the majority of other literature.
#'
#' Taking either the first or the second derivative are both valid approaches, but we prefer the reasoning behind taking the second derivative. After crossing a threshold, ventilation increases out of proportion to VO2 (VT1) or to VCO2 (VT2). After the threshold, the rate of increase is faster. Put another way, the slope after the breakpoint is higher than the slope before the breakpoint. However, since we are more interested in when that \emph{slope changes}, we want to know the slope of the slope, i.e., the acceleration, or the second derivative.
#'
#' Following the logic of Leo et al. (2017), this method then finds the local maxima in the acceleration because this shows when the slope is changing most rapidly. The idea is that above and below each threshold, our body is in two different states, and by finding the where slope changes most rapidly, we are locating that tipping point. Although other research shows that "thresholds" are \emph{not} tipping points but rather \emph{transitions} with gray zones between them (Ozkaya et al., 2022), it is nevertheless convenient to denote the threshold with a single point.
#'
#' To implement a method similar to Leo et al. (2017), this function actually takes four derivatives of the best-fit polynomial. We find the roots of the the third derivative to find the real parts of the local maxima and minima. For each of those real maxima and minima, we use the fourth derivative to determine if they represent maxima (negative sign), or minima (positive sign). The function then finds the highest of any real local maxima. The values associated with the closest data point to highest real local maxima are taken as the threshold.
#'
#' @param .data Gas exchange data frame or tibble.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param bp Is this the first (\code{vt1}) or the second (\code{vt2}) breakpoint?
#' @param degree The polynomial degree the function should fit to the curve. If left \code{NULL} this function finds the best fit.
#' @param vo2 Name of the \code{vo2} variable
#' @param vco2 Name of the \code{vco2} variable
#' @param ve Name of the \code{ve} variable
#' @param time Name of the \code{time} variable
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param pos_change Do you expect the slope to be increasing at the breakpoint? This helps with filtering maxima.
#' @param ci Should the output include confidence interval data? Default is `FALSE`.
#' @param conf_level Confidence level to use if calculating confidence intervals.
#' @param plots Should this function generate plots? Set to `FALSE` to save time.
#'
#' @returns A slice of the original data frame at the threshold index with a new `algorithm` column.
#'
#' @export
#'
#' @references
#' Leo, J. A., Sabapathy, S., Simmonds, M. J., & Cross, T. J. (2017). The Respiratory Compensation Point is Not a Valid Surrogate for Critical Power. Medicine and science in sports and exercise, 49(7), 1452-1460.
#' Ozkaya, O., Balci, G. A., As, H., Cabuk, R., & Norouzi, M. (2022). Grey zone: a gap between heavy and severe exercise domain. Journal of Strength and Conditioning Research, 36(1), 113-120.
#' Wisén, A. G., & Wohlfart, B. (2004). A refined technique for determining the respiratory gas exchange responses to anaerobic metabolism during progressive exercise–repeatability in a group of healthy men. Clinical physiology and functional imaging, 24(1), 1-9.
#' Wisén, A. G., & Wohlfart, B. (2004). Aerobic and functional capacity in a group of healthy women: reference values and repeatability. Clinical physiology and functional imaging, 24(6), 341-351.
#'
#'
#' @examples
#' # TODO write an example
#'
d2_poly_reg_maxima <- function(.data,
                            .x,
                            .y,
                            bp,
                            ...,
                            degree = NULL,
                            vo2 = "vo2",
                            vco2 = "vco2",
                            ve = "ve",
                            time = "time",
                            alpha_linearity = 0.05,
                            pos_change = TRUE,
                            ordering = c("by_x", "time"),
                            # TODO ADD FRONT TRIM ARGUMENTS,
                            ci = FALSE,
                            conf_level = 0.95,
                            plots = TRUE
                            ) {

    # check if there is crucial missing data
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))

    ordering <- match.arg(ordering, several.ok = FALSE)
    .data <- order_cpet_df(.data, .x = .x , time = time,
                           ordering = ordering)

    lm_poly <- loop_d2_poly_reg_maxima(.data = .data, .x = .x, .y = .y,
                             degree = degree,
                             alpha_linearity = alpha_linearity)

    equi_spaced_x <- seq(from = min(.data[[.x]]), to = max(.data[[.x]]),
                         length.out = range(.data[[vo2]]) %>%
                             diff() %>%
                             round())

    pred <- stats::predict(lm_poly,
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

    # get values at threshold
    if(length(sign_change_idx) > 0) {
        x_threshold <- equi_spaced_x[sign_change_idx]
        y_hat_threshold <- stats::predict(lm_poly,
                                          tibble::tibble("{.x}" := x_threshold))

        bp_dat <- find_threshold_vals(.data = .data, thr_x = x_threshold,
                                      thr_y = y_hat_threshold, .x = .x,
                                      .y = .y, ...)

        # get values at threshold
        bp_dat <- bp_dat %>%
            dplyr::mutate(bp = bp,
                          algorithm = "d2_reg_spline_maxima",
                          x_var = .x,
                          y_var = .y,
            )
    } else {
        bp_dat <- tibble::tibble()
    }

    if(nrow(bp_dat) == 0) { # no breakpoint found
        bp_dat <- bp_dat %>%
            dplyr::add_row() %>%
            dplyr::mutate(determinant_bp = FALSE)
        # add character or factor columns that are all the same value (e.g. ids)
        non_numeric_df <- .data %>%
            dplyr::select(tidyselect::where(function(x) is.character(x) | is.factor(x) &
                                                all(x == x[1]))) %>%
            dplyr::slice(1)
        bp_dat <- dplyr::bind_cols(bp_dat, non_numeric_df)
    } else { # breakpoint found
        bp_dat <- bp_dat %>%
            dplyr::mutate(determinant_bp = TRUE)
    }

    bp_dat <- bp_dat %>%
        dplyr::mutate(algorithm = "d2_poly_reg_maxima",
                      est_ci = "estimate",
                      x_var = .x,
                      y_var = .y,
                      bp = bp)

    if(ci) {
        # nonparametric bootstrapping
        boot_res <- boot::boot(data = .data,
                               statistic = boot_ci_d2_poly_reg,
                               R = 100,
                               sim = "ordinary",
                               .x = .x,
                               .y = .y,
                               degree = degree,
                               pos_change = pos_change,
                               alpha_linearity = alpha_linearity,
                               parallel = "multicore")
        # calculate percentile CI's
        boot_ci <- broom::tidy(boot_res,
                               conf.int = TRUE,
                               conf.method = "perc")
        # calculate new values at threshold
        ci_lower_x <- boot_ci$conf.low
        ci_lower_y <- stats::predict(lm_poly,
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
        ci_upper_y <- stats::predict(lm_poly,
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
        if(any(bp_dat[["determinant_bp"]])) {
            bp_dat[["determinant_bp"]] <- TRUE
        }
    }

    bp_dat <- bp_dat %>%
        dplyr::relocate(bp, algorithm, x_var, y_var, est_ci, determinant_bp)

    if(plots) {
        bp_plot <- ggplot2::ggplot(data = .data,
                                   ggplot2::aes(x = .data[[.x]], y = .data[[.y]])) +
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
                lm_poly_reg = lm_poly,
                deriv_func = spline_func,
                # deriv2_expr = deriv2,
                bp_plot = bp_plot))
}

#' @keywords internal
loop_d2_poly_reg_maxima <- function(.data, .x, .y,
                          degree = NULL, alpha_linearity = 0.05) {

    if (!is.null(degree)) {
        # if the user specifies a degree, find that and be done with it
        lm_poly <- paste0(.y, " ~ ", "1 + ", "poly(", .x, ", degree = ",
                          degree, ", raw = TRUE)") %>%
            stats::lm(data = .data)
    } else {
        # if the user does NOT specify a degree (default),
        # find the best degree using likelihood ratio test
        degree = 5 # start at degree = 5 so you can take 4 derivatives
        lm_list = vector(mode = "list", length = 0) # hold lm data
        # keep raw = TRUE for now b/c it's way easier for me to test
        # if the derivative expressions are working
        # starting with degree = 5 b/c that's the minimum I've seen in
        # other papers that seemed to use polynomial fits
        lm_poly <- paste0(.y, " ~ ", "1 + ", "poly(", .x, ", degree = ",
                          degree, ", raw = TRUE)") %>%
            stats::lm(data = .data)
        lm_list <- append(lm_list, list(lm_poly))

        cont <- TRUE
        i <- 1 # start at 2 b/c we already made linear (degree = 1) model
        while(cont == TRUE) {
            lm_poly <- paste0(.y, " ~ ", "1 + ", "poly(", .x, ", degree = ",
                              degree + i, ", raw = TRUE)") %>%
                stats::lm(data = .data)
            lm_list <- append(lm_list, list(lm_poly))
            lrt <- stats::anova(lm_list[[i]], lm_list[[i+1]])
            if (is.na(lrt$`Pr(>F)`[2]) | lrt$`Pr(>F)`[2] >= alpha_linearity) {
                cont = FALSE
                lm_poly <- lm_list[[i]] # take the previous model
            }
            i <- i + 1
        }
    }

    lm_poly
}


boot_ci_d2_poly_reg <- #' @keywords internal
    bootstrap_ci_d2_reg_spline <- function(.data,
                                           i,
                                           .x,
                                           .y,
                                           degree,
                                           pos_change,
                                           alpha_linearity) {
        # find best number of knots
        lm_poly <- loop_d2_poly_reg_maxima(.data = .data[i,],
                                             .x = .x,
                                             .y = .y,
                                             degree = degree,
                                             alpha_linearity = alpha_linearity)

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
        pred <- stats::predict(lm_poly,
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




# OLD CODE USING SYMBOLIC DIFFERENTIATION

# TODO what about if the best model is linear? Then this method probably
# won't work well and you should get a warning or an error about that
# also, I think you need a minimum of a 4th order equation if you want to
# find the maxima in the acceleration (2nd derivative) of the relationship
# b/c you essentially need to take a third derivative and still have an x term

# poly_expr <- expr_from_coefs(lm_poly$coefficients)
# deriv1 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 1) # slope
# deriv2 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 2) # acceleration
# deriv3 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 3) # jerk. Use to find maxima in 2nd deriv

# find x-values of roots
# roots_deriv3 <- deriv3 %>%
#     Ryacas::y_fn("Simplify") %>% # simplify derivative for ease of extraction coefficients
#     Ryacas::yac_str() %>%
#     stringr::str_extract_all("(?<!\\^)-?\\d+\\.?\\d*e?-?\\d*") %>% # extract coefficients
#     purrr::map(as.numeric) %>%
#     unlist() %>%
#     rev() %>% # reverse order for polyroot()
#     polyroot() %>%
#     find_real()

# filter by roots within range of x values
# roots_deriv3 <- roots_deriv3[roots_deriv3 >= min(.data[[.x]]) &
#                                  roots_deriv3 <= max(.data[[.x]])]

# use pos_change to search for values above or below 0
# if(pos_change) {
#     # use 4th derivative to find which are maxima (y @ roots < 0)
#     local_maxima_x <- roots_deriv3[grad(expr_to_func(deriv3),
#                                         roots_deriv3) < 0]
#     # uncomment if we need local_maxima_y ever later
#     # local_maxima_y <- eval(deriv2, envir = list(x = local_maxima_x))
# } else {
#     # use 4th derivative to find which are minima (y @ roots > 0)
#     local_maxima_x <- roots_deriv3[grad(expr_to_func(deriv3),
#                                         roots_deriv3) > 0]
    # local_maxima_y <- eval(deriv2, envir = list(x = local_maxima_x))
# }


# you would generally expect the threshold to be the highest of any possible
# local maxima (at least for VT2). In this case, highest means highest x-val
# (and usually highest y-val, too)

# select higher of the two local maxima b/c we expect just one large change
# I feel like this goes against my visual detection
# according to Cross et al. (2012), one should select the lower of any two
# maxima for VT1 and the higher of any two maxima for VT2
# if(length(local_maxima_x) > 1) {
#     if(bp == "vt1") {
#         local_maxima_x <- local_maxima_x[which.min(local_maxima_x)]
#     }
#     if(bp == "vt2") {
#         # techincally, this is a local minima
#         local_maxima_x <- local_maxima_x[which.max(local_maxima_x)]
#     }
# }
#
# y_hat_threshold <- eval(poly_expr, envir = list(x = local_maxima_x))

# if(length(local_maxima_x) > 0) {
#     x_threshold <- equi_spaced_x[sign_change_idx]
#     y_hat_threshold <- stats::predict(lm_poly,
#                                       tibble::tibble("{.x}" := x_threshold))
#     bp_dat <- find_threshold_vals(.data = .data, thr_x = local_maxima_x,
#                                   thr_y = y_hat_threshold, .x = .x,
#                                   .y = .y, ...)
# } else {
#     bp_dat <- tibble::tibble()
# }

# make better plotting data by making more, evenly spaced points
# equi_spaced_x <- seq(from = min(.data[[.x]]), to = max(.data[[.x]]),
#                      length.out = range(.data[[vo2]]) %>%
#                          diff() %>%
#                          round())
#
# pred <- stats::predict(lm_poly,
#                        newdata = tibble::tibble("{.x}" := equi_spaced_x))
