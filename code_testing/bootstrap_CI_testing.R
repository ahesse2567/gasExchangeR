x_threshold_boot <- function(data, i, .x, .y, df, degree) {

    # find best number of knots
    lm_spline <- loop_d2_reg_spline(.data = data[i,], .x = .x, .y = .y,
                                    df = df, degree = degree)

    # It feels like there's a better way than splinefun to find a more exact
    # 2nd derivative, but I don't know it

    # find maxima in 2nd derivative. This is accomplished by making
    # a spline interpolation function. Use this with optimize() function.

    # get new values at equal spacing for a smoother splinefun result
    equi_spaced_x <- seq(from = min(data[i,][[.x]]),
                         to = max(data[i,][[.x]]),
                         length.out = range(data[i,][[vo2]]) %>%
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
                # remove the left-most value
                neg_slope_idx <- neg_slope_idx[-1]
                sign_change_idx <- sign_change_idx[neg_slope_idx]
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
                # remove the left-most value
                pos_slope_idx <- pos_slope_idx[-1]
                sign_change_idx <- sign_change_idx[pos_slope_idx]
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

    dplyr::if_else(length(sign_change_idx) > 0,
            equi_spaced_x[sign_change_idx],
            NA)
}
copy_df <- .data

boot::boot(data = copy_df,
           statistic = x_threshold_boot,
           R = 10,
           .x = .x,
           .y = .y,
           df = df,
           degree = degree)
