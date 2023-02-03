#' Second derivative inflection
#'
#' Find a breakpoint in ventilatory data based on the 2nd derivative inflection point after Santos et al. (2004). In their paper, they use pulmonary ventilation (VE) vs. VO2 to find the first ventilatory (VT1) / gas exchange threshold (GET) and VE vs. VCO2 to find the respiratory compensation point (RCP) / second ventilatory threshold (VT2). It is generally suitable to also use VE vs. time to find VT2. One can truncate the test and rerun VE vs. time to find VT1.
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

#'
#' @returns A data frame at the threshold, with a few columns added to note the algorithm used.
#' @export
#'
#' @details To mimic the Santos et al. (2004) paper, this fits a 5th-order polynomial regression to the VE vs. VO2 or VCO2 data. Users can set `degree = NULL` to iteratively increase the polynomial if the likelihood ratio test is significant. After establishing the polynomial, we find the real roots of the 2nd derivative within the given range of the `.x` variable to find the inflection points. We then filter by concavity changes that go from up to down because this matches the expected shape of the pulmonary ventilation curves.
#'
#' We are unsure how to handle multiple inflection points at this time, but this is unlikely to be an issue if the polynomial degree is kept equal to 5.
#'
#' @references
#' Santos, E. L., & Giannella-Neto, A. (2004). Comparison of computerized methods for detecting the ventilatory thresholds. European journal of applied physiology, 93, 315-324.
#'
#' @examples
#' # TODO write an example
#'
d2_inflection <- function(.data,
                          .x,
                          .y,
                          bp,
                          degree = 5,
                          vo2 = "vo2",
                          vco2 = "vco2",
                          ve = "ve",
                          time = "time",
                          alpha_linearity = 0.05) # change to just alpha?
{

    # check if there is crucial missing data
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))

    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])

    lm_poly <- loop_d2_inflection(.data = .data, .x = .x, .y = .y,
                                  degree = degree,
                                  alpha_linearity = alpha_linearity)

    poly_expr <- expr_from_coefs(lm_poly$coefficients)
    deriv1 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 1) # slope
    deriv2 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 2) # acceleration

    # find x-values of roots of second derivative (inflection points)
    inflection_points <- deriv2 %>%
        Ryacas::y_fn("Simplify") %>% # simplify derivative for ease of extraction coefficients
        Ryacas::yac_str() %>%
        stringr::str_extract_all("(?<!\\^)-?\\d+\\.?\\d*e?-?\\d*") %>% # extract coefficients
        purrr::map(as.numeric) %>%
        unlist() %>%
        rev() %>% # reverse order for polyroot()
        polyroot() %>%
        find_real()

    # filter by roots within range of x values
    inflection_points <- inflection_points[inflection_points >= min(.data[[.x]]) &
                                               inflection_points <= max(.data[[.x]])]

    concavity_chgs <- concavity_changes(poly_expr,
                                        inflection_points = inflection_points)
    # filter by concave up to concave down
    inflection_points <- inflection_points[concavity_chgs == "up to down"]

    # what do we do if there's more than one inflection point?
    # do we take the one with the highest 1st derivative slope? But then,
    # how different is that from just finding the highest slope?
    threshold_idx <- which.min(abs(.data[[.x]] - inflection_points))

    # get values at threshold
    bp_dat <- .data[threshold_idx,] %>%
        dplyr::mutate(bp = bp,
                      algorithm = "d2_inflection",
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
    bp_dat

}

#' @keywords internal
loop_d2_inflection <- function(.data, .x, .y,
                               degree = NULL, alpha_linearity = 0.05) {
    # if the user specifies a degree, find that and be done with it
    if (!is.null(degree)) {
        lm_poly <- paste0(.y, " ~ ", "1 + ",
                          "poly(", .x, ", degree = ", degree, ", raw = TRUE)") %>%
            as.formula() %>%
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
            as.formula() %>%
            stats::lm(data = .data)

        lm_list <- append(lm_list, list(lm_poly))

        cont <- TRUE
        i <- 1 # start at 2 b/c we already made linear (degree = 1) model
        while(cont == TRUE) {
            lm_poly <- paste0(.y, " ~ ", "1 + ",
                              "poly(", .x, ", degree = ", degree + i, ", raw = TRUE)") %>%
                as.formula() %>%
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

#' Find when a function changes from concave up to concave down or vice versa
#' @keywords internal
#' @noRd
#'
concavity_changes <- function(f, inflection_points, tol = 0.1) {
    # find second derivative
    f_dd <- Deriv::Deriv(f, x = "x", nderiv = 2)
    if (is.expression(f_dd)) { # change expression to function if needed
        f_dd <- expr_to_func(f_dd)
    }
    left_sign <- sign(f_dd(inflection_points - tol))
    right_sign <- sign(f_dd(inflection_points + tol))
    conc_changes <- if_else(left_sign > 0 & right_sign < 0,
                            "up to down", "down to up")
    conc_changes
}
