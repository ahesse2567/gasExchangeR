x_threshold_boot <- function(.data, i, .x, .y, df, degree) {

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
copy_df <- .data

boot_res <- boot::boot(data = copy_df,
           statistic = x_threshold_boot,
           R = 500,
           sim = "ordinary",
           .x = .x,
           .y = .y,
           df = df,
           degree = degree,
           parallel = "multicore")
plot(boot_res)
boot::boot.ci(boot_res, type = "perc")
boot_ci <- boot::boot.ci(boot_res, type = "perc",
                         h = log,
                         hinv = exp)
x_threshold

se_boot <- sd(boot_res$t) / sqrt(length(boot_res$t))
ci_95 <- x_threshold + c(-1, 1) * 1.96 * se_boot
ci_95

x <- c(0.2, 0.528, 0.11, 0.260, 0.091,
       1.314, 1.52, 0.244, 1.981, 0.273,
       0.461, 0.366, 1.407, 0.79, 2.266)


x <- boot_res$t


boxcox(model, data = x)

quantile(boot_res$t, probs = c(2.5, 97.5)/100)

library(MASS)
b <- boxcox(lm(x ~ 1), data = data.frame(x = x))
b

boxcox(lm(boot_res$t ~ 1, data = data.frame(x = boot_res$t)))

b <- boxcox(lm(x ~ 1, data = data.frame(x = boot_res$t)))
b
lambda <- b$x[which.max(b$y)]
lambda

hist((boot_res$t ^ lambda - 1) / lambda)
