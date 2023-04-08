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
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param ci Should the output include confidence interval data? Default is `FALSE`.
#' @param conf_level Confidence level to use if calculating confidence intervals.
#' @param plots Should this function generate plots? Set to `FALSE` to save time.
#'
#' @returns A list containing breakpoint data, best-fit and related functions, and a plot.
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
                          ...,
                          degree = NULL,
                          vo2 = "vo2",
                          vco2 = "vco2",
                          ve = "ve",
                          time = "time",
                          alpha_linearity = 0.05,
                          ordering = c("by_x", "time"),
                          # TODO ADD FROM TRIM ARGUMENTS,
                          ci = FALSE,
                          conf_level = 0.95,
                          plots = TRUE
                          ) {

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
    inflection_points_x <- deriv2 %>%
        Ryacas::y_fn("Simplify") %>% # simplify derivative for ease of extraction coefficients
        Ryacas::yac_str() %>%
        stringr::str_extract_all("(?<!\\^)-?\\d+\\.?\\d*e?-?\\d*") %>% # extract coefficients
        purrr::map(as.numeric) %>%
        unlist() %>%
        rev() %>% # reverse order for polyroot()
        polyroot() %>%
        find_real()

    # filter by roots within range of x values
    inflection_points_x <- inflection_points_x[inflection_points_x >=
                                                   min(.data[[.x]]) &
                                               inflection_points_x <=
                                                   max(.data[[.x]])]

    concavity_chgs <- concavity_changes(poly_expr,
                                        inflection_points = inflection_points_x)
    # filter by concave up to concave down
    inflection_points_x <- inflection_points_x[concavity_chgs == "up to down"]
    # if there's more than one up to down inflection point, choose whichever
    # point has the higher slope
    inflection_points_x <- inflection_points_x[which.max(
        eval(deriv1, envir = list(x = inflection_points_x)))]

    inflection_points_y <- eval(poly_expr, envir = list(x = inflection_points_x))

    if(length(inflection_points_x) > 0 & length(inflection_points_y) > 0) {
        bp_dat <- find_threshold_vals(.data = .data, thr_x = inflection_points_x,
                                      thr_y = inflection_points_y, .x = .x,
                                      .y = .y, ...)
    } else {
        bp_dat <- tibble::tibble()
    }

    # get values at threshold. Make this a deriv helpers function to prep bp_dat?

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
        dplyr::mutate(bp = bp,
                      algorithm = "d2_inflection",
                      x_var = .x,
                      y_var = .y) %>%
        dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp)

    # make plot for output
    if(plots) {
        bp_plot <- ggplot2::ggplot(data = .data,
                                   ggplot2::aes(x = .data[[.x]], y = .data[[.y]])) +
            ggplot2::geom_point(alpha = 0.5) +
            ggplot2::geom_line(ggplot2::aes(y = eval(poly_expr,
                                                     envir = list(x = .data[[.x]])))) +
            ggplot2::geom_vline(xintercept = bp_dat[[.x]]) +
            ggplot2::theme_minimal()
    } else {
        bp_plot <- NULL
    }

    return(list(breakpoint_data = bp_dat,
                lm_poly_reg = lm_poly,
                deriv1_expr = deriv1,
                deriv2_expr = deriv2,
                bp_plot = bp_plot))
}

#' @keywords internal
loop_d2_inflection <- function(.data, .x, .y,
                               degree = NULL, alpha_linearity = 0.05) {
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
    conc_changes <- dplyr::if_else(left_sign > 0 & right_sign < 0,
                            "up to down", "down to up")
    conc_changes
}
