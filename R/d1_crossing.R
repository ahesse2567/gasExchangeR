#' Find the gas exchange threshold (GET) (~1st ventilatory threshold (VT1)) using the 1st derivative crossing point
#'
#' Use polynomial regression to find the gas exchange threshold (GET) (practically equivalent to the first ventilatory threshold (VT1)). This method by Wisén and Wohlfart (2004) is similar to the V-slope method by Beaver et al. (1986) in that it finds when the rate of CO2 production outpaces the rate of O2 production. In contrast to piecewise regression in the original V-slope method, this uses the first derivative of VO2 and VCO2 vs. time. The threshold is when dVCO2/dt surpasses dVO2/dt for the final time.
#'
#' @param .data Gas exchange data frame or tibble.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param degree The polynomial degree the function should fit to the curve. If left \code{NULL} this function finds the best fit.
#' @param vo2 Name of the \code{vo2} variable
#' @param vco2 Name of the \code{vco2} variable
#' @param ve Name of the \code{ve} variable
#' @param time Name of the \code{time} variable
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.
#' @param bp Is this the first (\code{vt1}) or the second (\code{vt2}) breakpoint? This method is specifically designed for finding VT1 (GET)
#' @param pos_change Do you expect the change in slope to be positive (default) or negative? If a two-line regression explains significantly reduces the sum square error but the change in slope does not match the expected underlying physiology, the breakpoint will be classified as indeterminate.
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param ci Should the output include confidence interval data? Default is `FALSE`.
#' @param conf_level Confidence level to use if calculating confidence intervals.
#' @param plots Should this function generate plots? Set to `FALSE` to save time.
#'
#' @returns A list include breakpoint data, best-fit and related functions, and a plot.
#'
#' @export
#'
#' @references
#' Beaver, W. L., Wasserman, K. A. R. L. M. A. N., & Whipp, B. J. (1986). A new method for detecting anaerobic threshold by gas exchange. Journal of applied physiology, 60(6), 2020-2027.
#' Wisén, A. G., & Wohlfart, B. (2004). A refined technique for determining the respiratory gas exchange responses to anaerobic metabolism during progressive exercise–repeatability in a group of healthy men. Clinical physiology and functional imaging, 24(1), 1-9.
#'
#' @examples
#'
#' # TODO write an example
#'
d1_crossing <- function(.data,
                        .x,
                        .y,
                        bp = "vt1",
                        ...,
                        degree = NULL,
                        vo2 = "vo2",
                        vco2 = "vco2",
                        ve = "ve",
                        time = "time",
                        alpha_linearity = 0.05, # change to just alpha?
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

    # best-fit polynomial for vo2
    lm_poly_vo2 <- loop_poly_d1_crossing(.data = .data, .x = time, .y = vo2,
                                 degree = degree,
                                 alpha_linearity = alpha_linearity)

    # best-fit polynomial for vco2
    lm_poly_vco2 <- loop_poly_d1_crossing(.data = .data, .x = time, .y = vco2,
                                  degree = degree,
                                  alpha_linearity = alpha_linearity)
    # 1st derivative for vo2
    poly_expr_vo2 <- expr_from_coefs(lm_poly_vo2$coefficients)
    deriv1_vo2 <- Deriv::Deriv(poly_expr_vo2, x = "x", nderiv = 1) # slope

    # 1st derivative for vco2
    poly_expr_vco2 <- expr_from_coefs(lm_poly_vco2$coefficients)
    deriv1_vco2 <- Deriv::Deriv(poly_expr_vco2, x = "x", nderiv = 1) # slope

    # turn derivative expressions into functions to use with uniroot.all()
    deriv1_vo2_func <- expr_to_func(deriv1_vo2)
    deriv1_vco2_func <- expr_to_func(deriv1_vco2)

    # calculate derivative of difference in 1st derivative functions
    # this lets us know if vco2 was surpassing vo2 or the other way around
    deriv_diff_deriv1 <- paste0(deriv1_vo2 %>%
                                    Ryacas::y_fn("Simplify") %>%
                                    Ryacas::yac_str(),
                                " - (",
                                deriv1_vco2 %>%
                                    Ryacas::y_fn("Simplify") %>%
                                    Ryacas::yac_str(),
                                ")") %>%
        parse(text = .) %>%
        Deriv::Deriv(x = "x", nderiv = 1)

    roots <- rootSolve::uniroot.all(
        function(x) deriv1_vo2_func(x) - deriv1_vco2_func(x),
        interval = c(min(.data[[.x]]), max(.data[[.x]])))

    # filter by roots within range of x values
    roots <- roots[roots >= min(.data[[.x]]) &
                       roots <= max(.data[[.x]])]

    # filter by which roots have a negative derivative. This indicates that
    # co2 is rising above o2. Finding max finds the last time this occurs.
    final_crossing <- roots[eval(deriv_diff_deriv1,
                                 envir = list(x = roots)) < 0] %>%
        max()

    # find the fitted VO2 value and use that to find threshold values.
    # You could use VCO2 or somehow integrate both VO2 & VCO2, but that might be a future update

    if(length(final_crossing) > 0) {
        y_hat_threshold <- eval(poly_expr_vo2,
                                envir = list(x = final_crossing))
        # get values at threshold
        bp_dat <- find_threshold_vals(.data = .data, thr_x = final_crossing,
                                      thr_y = y_hat_threshold, .x = .x,
                                      .y = .y, ...)

        bp_dat <- bp_dat %>%
            dplyr::mutate(bp = bp,
                          algorithm = "d1_crossing",
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
    } else { # breakpoint found
        bp_dat <- bp_dat %>%
            dplyr::mutate(determinant_bp = TRUE)
    }
    bp_dat <- bp_dat %>%
        dplyr::mutate(algorithm = "d1_crossing",
                      x_var = .x,
                      y_var = .y,
                      bp = bp) %>%
        dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp)

    if(plots) {
        bp_plot <- ggplot2::ggplot(data = .data,
                                   ggplot2::aes(x = .data[[.x]],
                                                y = .data[[.y]])) +
            ggplot2::geom_point(ggplot2::aes(y = vo2, color = "vo2"), alpha = 0.5) +
            ggplot2::geom_point(ggplot2::aes(y = vco2, color = "vco2"), alpha = 0.5) +
            ggplot2::geom_line(ggplot2::aes(y = eval(poly_expr_vo2,
                                                     envir = list(x = .data[[.x]])))) +
            ggplot2::geom_line(ggplot2::aes(y = eval(poly_expr_vco2,
                                                     envir = list(x = .data[[.x]])))) +
            ggplot2::geom_vline(xintercept = bp_dat[[.x]]) +
            ggplot2::scale_color_manual(values = c("vo2" = "red", "vco2" = "blue")) +
            ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) +
            ggplot2::theme_minimal()
    } else {
        bp_plot <- NULL
    }

    return(list(breakpoint_data = bp_dat,
                lm_poly_vo2 = lm_poly_vo2,
                deriv1_vo2 = deriv1_vo2,
                lm_poly_vco2 = lm_poly_vco2,
                deriv1_vco2 = deriv1_vco2,
                bp_plot = bp_plot))

}

#' @keywords internal
loop_poly_d1_crossing <- function(.data, .x, .y,
                          degree = NULL, alpha_linearity = 0.05) {
    # browser()
    # if the user specifies a degree, find that and be done with it
    if (!is.null(degree)) {
        lm_poly <- paste0(.y, " ~ ", "1 + ",
                          "poly(", .x, ", degree = ", degree, ", raw = TRUE)") %>%
            stats::lm(data = .data)
        # if the user does NOT specify a degree, find the best degree using
        # likelihood ratio test
    } else {
        degree = 5 # from testing and previous papers it seems like you need to
        # force a higher derivative if you want a local maxima within the
        # range of x values. That doesn't feel wonderful.
        lm_list = vector(mode = "list", length = 0) # hold lm data
        # keep raw = TRUE for now b/c it's way easier for me to test
        # if the derivative expressions are working
        # starting with degree = 5 b/c that's the minimum I've seen in
        # other papers that seemed to use polynomial fits

        lm_poly <- paste0(.y, " ~ ", "1 + ",
                          "poly(", .x, ", degree = ", degree, ", raw = TRUE)") %>%
            stats::lm(data = .data)

        lm_list <- append(lm_list, list(lm_poly))

        cont <- TRUE
        i <- 1 # start at 2 b/c we already made linear (degree = 1) model
        while(cont == TRUE) {
            lm_poly <- paste0(.y, " ~ ", "1 + ",
                              "poly(", .x, ", degree = ",
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




## plotting code to figure out later

# max_diff_deriv1 <- max(eval(deriv_diff_deriv1, envir = list(x = df_avg[["time"]])))
#
# sum_range_deriv1 <- eval(deriv_diff_deriv1,
#                          envir = list(x = df_avg[["time"]])) %>%
#     range() %>%
#     abs() %>%
#     sum()
#
# scale_factor <- max(c(.data[[vo2]], .data[[vco2]])) /
#     sum_range_deriv2
#
#
# bp_plot +
#     scale_y_continuous(name = paste(vo2, "&", vco2),
#                        sec.axis = sec_axis(
#                            ~ . / scale_factor*2 - sum_range_deriv2,
#                            name = paste("1st Derivative of",
#                                         vo2, "%", vco2))) +
#     geom_line(ggplot2::aes(y = eval(deriv_diff_deriv1,
#                            envir = list(x = df_avg[["time"]]))),
#               linetype = "dashed")

# this doesn't give the negative values we'll need though
