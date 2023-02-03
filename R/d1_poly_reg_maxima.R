#' Find a breakpoint using the maxima in the 1st derivative
#'
#' Use polynomial regression of increasing order to find the best-fit line for the data, take the 1st derivative, and find the highest maxima of the 1st derivative. This method is traditionally used with ventilatory equivalents (VE/VO2 or VE/VCO2), PetO2, or PetCO2 vs. time or vs. VO2. This method was specifically implemented by Santos et al. (2004).
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
#' @param make_plot Make plot showing the breakpoint?
#' @param show_plot Show the plot showing the breakpoint? Default = `FALSE`.
#'
#' @returns A slice of the original data frame at the threshold index with a new `algorithm` column.
#'
#' @export
#'
#' @references
#' Santos, E. L., & Giannella-Neto, A. (2004). Comparison of computerized methods for detecting the ventilatory thresholds. European journal of applied physiology, 93, 315-324.
#'
#'
#' @examples
#' # TODO write an example
#'
d1_poly_reg_maxima <- function(.data,
                               .x,
                               .y,
                               bp,
                               degree = 5,
                               vo2 = "vo2",
                               vco2 = "vco2",
                               ve = "ve",
                               time = "time",
                               alpha_linearity = 0.05, # change to just alpha?
                               make_plot = TRUE,
                               show_plot = FALSE)
    {
    # check if there is crucial missing data
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
    # browser()
    .data <- .data %>% # rearrange by x variable. Use time var to break ties.
        dplyr::arrange(.data[[.x]], .data[[time]])

    lm_poly <- loop_d1_poly_reg_maxima(.data = .data, .x = .x, .y = .y,
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

    # find x-values of roots
    roots_deriv2 <- deriv2 %>%
        Ryacas::y_fn("Simplify") %>% # simplify derivative for ease of extraction coefficients
        Ryacas::yac_str() %>%
        stringr::str_extract_all("(?<!\\^)-?\\d+\\.?\\d*e?-?\\d*") %>% # extract coefficients
        purrr::map(as.numeric) %>%
        unlist() %>%
        rev() %>% # reverse order for polyroot()
        polyroot() %>%
        find_real()

    # filter by roots within range of x values
    roots_deriv2 <- roots_deriv2[roots_deriv2 >= min(.data[[.x]]) &
                                     roots_deriv2 <= max(.data[[.x]])]

    # use 4th derivative to find which are maxima (y @ roots < 0)
    local_maxima_x <- roots_deriv2[eval(deriv3,
                                        envir = list(x = roots_deriv2)) < 0]
    local_maxima_y <- eval(deriv1, envir = list(x = local_maxima_x))

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
                      algorithm = "d1_poly_reg_maxima",
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

    if(make_plot) {
        bp_plot <- ggplot2::ggplot(data = .data, aes(x = .data[[.x]], y = .data[[.y]])) +
            geom_point(alpha = 0.5) +
            geom_line(aes(y = lm_poly$fitted.values)) +
            geom_vline(xintercept = bp_dat[[.x]]) +
            theme_minimal()
    }
    if(show_plot & make_plot) {
        plot(bp_plot)
    }
    browser()
    out_list <- list(breakpoint_data = bp_dat,
                     lm_poly = lm_poly)
    if(make_plot) {
        out_list$plot <- bp_plot
    } else {
        out_list$plot <- NULL
    }

    out_list
}

#' @keywords internal
loop_d1_poly_reg_maxima <- function(.data, .x, .y,
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
