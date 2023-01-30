#' Find a breakpoint using the maxima in the 2nd derivative
#'
#' Use polynomial regression of increasing order to find the best-fit line for the data, take the second derivative, and find the highest maxima of the second derivative. This method is traditionally used with ventilatory equivalents (VE/VO2 or VE/VCO2), PetO2, or PetCO2 vs. time or vs. VO2.
#'
#' Note, polynomial regressions are usually inferior to regression splines or to smoothing splines because polynomial regressions because polynomial regressions can overfit the data more easily. However, we offer this method in case users want to reproduce previous works and while we develop this package.
#'
#' Unlike several published works, this method iteratively finds the best-fit polynomial regression by increasing the polynomial order and using the likelihood ratio test with the \code{anova} function. When the likelihood ratio test does not find a statistically significant different difference between the previous and the newest polynomial, the function uses the previous polynomial. However, uses can specify the polynomial \code{degree} if desired.
#'
#' #' @details
#' Our implementation is similar to two studies by Wisén and Wohlfart et al. in that this uses a polynomial regression. However, they only found the \emph{first} derivative of the VO2 and VCO2 vs. time curves and denote the threshold as their crossing point. Using ventilatory equivalents vs. time, they considered the threshold as when the \emph{first} derivative of each ventilatory equivalent increased above zero (for the last time). In contrast, the majority of other literature we are aware of that uses regression splines or smoothing splines finds the \emph{second} derivative. This function finds the second derivative by default to match the majority of other literature.
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
#'
#' @returns A slice of the original data frame at the threshold index with a new `algorithm` column.
#'
#' @importFrom Deriv Deriv
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
                            degree = NULL,
                            vo2 = "vo2",
                            vco2 = "vco2",
                            ve = "ve",
                            time = "time",
                            alpha_linearity = 0.05 # change to just alpha?
                            ) {
    # check if there is crucial missing data
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))

    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])

    lm_poly <- loop_poly_reg(.data = .data, .x = .x, .y = .y,
                             degree = degree,
                             alpha_linearity = alpha_linearity)

    # TODO what about if the best model is linear? Then this method probably
    # won't work well and you should get a warning or an error about that
    # also, I think you need a minimum of a 4th order equation if you want to
    # find the maxima in the acceleration (2nd derivative) of the relationship
    # b/c you essentially need to take a third derivative and still have an x term

    # lm_poly <- lm_list[[3]] # keep for testing, delete later
    # this definitely breaks if the best-fit polynomial order is too low
    poly_expr <- expr_from_coefs(lm_poly$coefficients)
    deriv1 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 1) # slope
    deriv2 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 2) # acceleration
    deriv3 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 3) # jerk
    deriv4 <- Deriv::Deriv(poly_expr, x = "x", nderiv = 4) # snap. This kinda feels hacky and dumb at this point

    # find x-values of roots
    roots_deriv3 <- deriv3 %>%
        Ryacas::y_fn("Simplify") %>% # simplify derivative for ease of extraction coefficients
        Ryacas::yac_str() %>%
        stringr::str_extract_all("(?<!\\^)-?\\d+\\.?\\d*e?-?\\d*") %>% # extract coefficients
        purrr::map(as.numeric) %>%
        unlist() %>%
        rev() %>% # reverse order for polyroot()
        polyroot() %>%
        find_real()

    # filter by roots within range of x values
    roots_deriv3 <- roots_deriv3[roots_deriv3 >= min(.data[[.x]]) &
                                     roots_deriv3 <= max(.data[[.x]])]

    # use 4th derivative to find which are maxima (y @ roots < 0)
    local_maxima_x <- roots_deriv3[eval(deriv4,
                                        envir = list(x = roots_deriv3)) < 0]
    local_maxima_y <- eval(deriv2, envir = list(x = local_maxima_x))

    # you would generally expect the threshold to be the highest of any possible
    # local maxima (at least for VT2)

    # select higher of the two local maxima b/c we expect just one large change
    # I feel like this goes against my visual detection
    # according to Cross et al. (2012), one should select the lower of any two
    # maxima for VT1 and the higher of any two maxima for VT2
    if(length(local_maxima_x) > 1) {
        if(bp == "vt1") {
            min_local_maxima_x <- local_maxima_x[which.min(local_maxima_y)]
            threshold_idx <-
                which.min(abs(.data[[.x]] - min_local_maxima_x))
        }
        if(bp == "vt2") {
            max_local_maxima_x <- local_maxima_x[which.max(local_maxima_y)]
            threshold_idx <- which.min(abs(.data[[.x]] - max_local_maxima_x))
        }
    } else {
        threshold_idx <- which.min(abs(.data[[.x]] - local_maxima_x))
    }

    # get values at threshold
    bp_dat <- .data[threshold_idx,] %>%
        dplyr::mutate(bp = bp,
               algorithm = "d2_poly_reg_maxima",
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
loop_poly_reg <- function(.data, .x, .y,
                          degree = NULL, alpha_linearity = 0.05) {
    # browser()
    # if the user specifies a degree, find that and be done with it
    if (!is.null(degree)) {
        lm_poly <- stats::lm(.data[[.y]] ~ 1 + poly(.data[[.x]], degree = degree,
                                             raw = TRUE), data = .data)
        # if the user does NOT specify a degree, find the best degree using
        # likelihood ratio test
    } else {
        degree = 5 # start at degree = 5 so you can take 4 derivatives
        lm_list = vector(mode = "list", length = 0) # hold lm data
        # keep raw = TRUE for now b/c it's way easier for me to test
        # if the derivative expressions are working
        # starting with degree = 5 b/c that's the minimum I've seen in
        # other papers that seemed to use polynomial fits
        lm_poly <- stats::lm(.data[[.y]] ~ 1 + poly(.data[[.x]], degree = degree,
                                             raw = TRUE), data = .data)
        lm_list <- append(lm_list, list(lm_poly))

        cont <- TRUE
        i <- 1 # start at 2 b/c we already made linear (degree = 1) model
        while(cont == TRUE) {
            lm_poly <- stats::lm(.data[[.y]] ~ 1 + poly(.data[[.x]],
                                                 degree = degree + i,
                                                 raw = TRUE),
                          data = .data)
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
