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
                                 alpha_linearity = 0.05, # change to just alpha?
                                 pos_change = TRUE,
                                 ...) {

    # check if there is crucial missing data
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))

    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])

    # find best number of knots
    lm_spline <- loop_d2_reg_spline(.data = .data, .x = .x, .y = .y,
                                    df = df, degree = degree)

    # It feels like there's a better way than splinefun to find a more exact
    # 2nd derivative, but I don't know it

    # find maxima in 2nd derivative. This is accomplished by making
    # a spline interpolation function. Use this with optimize() function.

    # get new values at equal spacing for a smoother splinefun result
    equi_spaced_x <- seq(from = min(.data[[.x]]), to = max(.data[[.x]]),
                         length.out = nrow(.data))
    pred <- predict(lm_spline, newdata = tibble::tibble("{.x}" := equi_spaced_x))

    spline_func <- stats::splinefun(x = equi_spaced_x, y = pred)
    # find number of maxima
    sign_change_idx <- slope_sign_changes(y = spline_func(x = .data[[.x]],
                                                          deriv = 2),
                                          change = "pos_to_neg")
    # filter by expected slope (usually positive)
    # there can be a spike in accel. during the initial drop in vent eqs
    # pos_change &
    if(pos_change) {
        sign_change_idx <- sign_change_idx[spline_func(
            x = .data[[.x]][sign_change_idx],
            deriv = 1) > 0]
    } else {
        sign_change_idx <- sign_change_idx[spline_func(
            x = .data[[.x]][sign_change_idx],
            deriv = 1) < 0]
    }

    y_val_sign_changes <- numeric(length = length(sign_change_idx))
    for(i in seq_along(sign_change_idx)) {
        interval <- c(.data[[.x]][sign_change_idx[i]-1],
                      .data[[.x]][sign_change_idx[i]+1])

        extrema_deriv2 <- stats::optimize(f = spline_func,
                                   interval = interval,
                                   deriv = 2, maximum = TRUE)
        y_val_sign_changes[i] <- extrema_deriv2$objective
    }

    # use highest maxima as threshold
    threshold_idx <- sign_change_idx[which.max(y_val_sign_changes)]

    # x_threshold <- .data[[.x]][threshold_idx]
    # y_hat_threshold <- predict(lm_spline,
    #                                tibble("{.x}" := x_threshold))

    # get values at threshold
    bp_dat <- .data[threshold_idx,] %>%
        dplyr::mutate(bp = bp,
                      algorithm = "d2_reg_spline_maxima",
                      x_var = .x,
                      y_var = .y,
                      # determinant_bp = determinant_bp,
                      # pct_slope_change = pct_slope_change,
                      # f_stat = f_stat,
                      # p_val_f = pf_two,
        ) %>%
        dplyr::relocate(bp, algorithm, x_var, y_var,
                        # determinant_bp,
                        # pct_slope_change, f_stat, p_val_f
        )

    bp_plot <- ggplot2::ggplot(data = .data,
                               aes(x = .data[[.x]], y = .data[[.y]])) +
        geom_point(alpha = 0.5) +
        geom_line(aes(x = equi_spaced_x, y = pred)) +
        geom_vline(xintercept = bp_dat[[.x]]) +
        theme_minimal()

    return(list(bp_dat = bp_dat,
                lm_reg_spline = lm_spline,
                spline_deriv_func = spline_func,
                bp_plot = bp_plot))
}

#' @keywords internal
loop_d2_reg_spline <- function(.data, .x, .y, df = NULL,
                               degree = 3, alpha_linearity = 0.05) {
    # TODO allow users to specify b-spline or natural-spline basis
    # would that use do.call()?

    # if statement for if users specify the knots or df
    if(!is.null(df)) {
        lm_spline <- paste0(.y, " ~ ", "1 + ",
                            "bs(", .x,
                            ", df = ", df,
                            ", degree = ", degree, ")") %>%
            stats::as.formula() %>%
            stats::lm(data = .data)
    } else {
        spline_mod_list = vector(mode = "list", length = 0)
        cont <- TRUE
        i <- 1
        # reference model with one interior knot
        lm_spline <- paste0(.y, " ~ ", "1 + ",
                            "bs(", .x,
                            ", df = ", i + degree,
                            ", degree = ", degree, ")") %>%
            stats::as.formula() %>%
            stats::lm(data = .data)
        spline_mod_list <- append(spline_mod_list, list(lm_spline))
        # while loop beings with 1 knot (3 df already used assuming a 3rd spline)
        while(cont == TRUE) {
            # browser()
            i <- i + 1
            # TODO add options for user-defined knots and df
            lm_spline <- paste0(.y, " ~ ", "1 + ",
                                "bs(", .x,
                                ", df = ", i + degree,
                                ", degree = ", degree, ")") %>%
                stats::as.formula() %>%
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
