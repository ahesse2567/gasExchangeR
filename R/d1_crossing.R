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
#'
#' @returns
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
                        degree = NULL,
                        vo2 = "vo2",
                        vco2 = "vco2",
                        ve = "ve",
                        time = "time",
                        alpha_linearity = 0.05, # change to just alpha?
                        pos_change = TRUE,
                        bp = "vt1") {

    # check if there is crucial missing data
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))

    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])

    # best-fit polynomial for vo2
    lm_poly_vo2 <- loop_poly_reg(.data = .data, .x = time, .y = vo2,
                                 degree = degree,
                                 alpha_linearity = alpha_linearity)

    # best-fit polynomial for vco2
    lm_poly_vco2 <- loop_poly_reg(.data = .data, .x = time, .y = vco2,
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
        Deriv(x = "x", nderiv = 1)

    roots <- rootSolve::uniroot.all(function(x) deriv1_vo2_func(x) - deriv1_vco2_func(x),
                                    interval = c(min(.data[[.x]]), max(.data[[.x]])))

    # filter by roots within range of x values
    roots <- roots[roots >= min(.data[[.x]]) &
                       roots <= max(.data[[.x]])]

    # filter by which roots have a negative derivative. This indicates that
    # co2 is rising above o2. Finding max finds the last time this occurs.
    final_crossing <- roots[eval(deriv_diff_deriv1, envir = list(x = roots)) < 0] %>%
        max()

    threshold_idx <- which.min(abs(.data[[.x]] - final_crossing))

    # get values at threshold
    bp_dat <- .data[threshold_idx,] %>%
        dplyr::mutate(bp = bp,
               algorithm = "d1_crossing",
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
    bp_dat # return breakpoint data
}

#' @keywords internal
loop_poly_reg <- function(.data, .x, .y,
                          degree = NULL, alpha_linearity = 0.05) {
    # browser()
    # if the user specifies a degree, find that and be done with it
    if (!is.null(degree)) {
        lm_poly <- paste0(.y, " ~ ", "1 + ",
                          "poly(", .x, ", degree = ", degree, ", raw = TRUE)") %>%
            as.formula() %>%
            lm(data = .data)
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
            as.formula() %>%
            lm(data = .data)

        lm_list <- append(lm_list, list(lm_poly))

        cont <- TRUE
        i <- 1 # start at 2 b/c we already made linear (degree = 1) model
        while(cont == TRUE) {
            lm_poly <- paste0(.y, " ~ ", "1 + ",
                              "poly(", .x, ", degree = ", degree + i, ", raw = TRUE)") %>%
                as.formula() %>%
                lm(data = .data)
            lm_list <- append(lm_list, list(lm_poly))
            lrt <- anova(lm_list[[i]], lm_list[[i+1]])
            if (is.na(lrt$`Pr(>F)`[2]) | lrt$`Pr(>F)`[2] >= alpha_linearity) {
                cont = FALSE
                lm_poly <- lm_list[[i]] # take the previous model
            }
            i <- i + 1
        }
    }

    lm_poly
}
