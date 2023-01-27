#' Find breakpoints in gas exchange data
#'
#' Breakpoint allows the user to specify particular x-y relationships, algorithms, and criteria for determining if a breakpoint was determinant or indeterminate. It is important to understand the order of steps the function takes.
#' 1. The function preselects x and y variables for the first and sometimes second threshold if the user indicates they wish to use a well-known method. For example, the v-slope method uses VE vs. VCO2 for finding VT2 (RC), and VCO2 vs. VO2 for finding VT1 (GET).
#' 2. Determine the breakpoint at VT2 and if it was determinant.
#' 3. If determinant, truncate the data at VT2.
#' 4. Determine the breakpoint at VT1 and if it was determinant.
#' 5. Return results.
#'
#' @section Warning:
#' It is strongly advised the user selects \emph{absolute VO2 of the same units as VCO2 (e.g. mL/min for both variables)} rather than relative VO2. Using relative VO2 will likely lead to incorrect results.
#'
#' @param .data Gas exchange data
#' @param method Specify \code{excess_co2}, \code{v-slope}, or \code{v-slope_simple} and this fills in the most common variables used for these methods for those not already specified.
#' @param algorithm_vt1 Algorithm to find VT1/GET
#' @param x_vt1 \code{x} variable to use to fine VT1/GET/aerobic threshold
#' @param y_vt1 \code{y} variable to use to fine VT1/GET/aerobic threshold
#' @param algorithm_vt2 Algorithm to find VT2/RC
#' @param x_vt2 \code{x} variable to use to fine VT2/RC/anaerobic threshold
#' @param y_vt2 \code{y} variable to use to fine VT2/RC/anaerobic threshold
#' @param vo2 The name of the vo2 column in \code{.data}
#' @param vco2 The name of the vco2 column in \code{.data}
#' @param ve The name of the ve column in \code{.data}
#' @param time The name of the time column in \code{.data}
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.
#' @param bps Should the function find the breakpoints for VT1, TVT2, or both. Default is \code{both}.
#' @param truncate By default this function truncates the data frame at VT2 prior to finding VT2. Change truncate to \code{FALSE} to use the entire data frame when searching for VT1.
#' @param pos_change Do you expect the change in slope to be positive (default) or negative? If a two-line regression explains significantly reduces the sum square error but the change in slope does not match the expected underlying physiology, the breakpoint will be classified as indeterminate.
#' @param ...
#'
#' @return
#' @export
#' @md
#'
#' @details
#' There are multiple ways to find a breakpoint. The first is by specifying a commonly used \code{method} and the appropriate \code{algorithm} and \code{x}, and \code{y}, variables will be selected for you. For example, entering 'v-slope' for \code{method} will pre-select 'v-slope' for \code{algorithm}, 'vo2' for \code{x_vt1}, and 'vco2' for \code{y_vt1}. If you enter a specific \code{method}, you can override the default \code{algorithm}, \code{x}, and \code{y}.
#'
#' Since there aren't common names for every way to find a breakpoint, you can also leave \code{method} as \code{NULL} and specify the \code{algorithm}, \code{x}, and \code{y}.
#'
#' This function uses the excess CO2 formula as described by Gaskill et al. (2001).
#'
#' @references
#' Gaskill, S. E., Ruby, B. C., Walker, A. V. A. J., Sanchez, O. A., Serfass, R. C., & Leon, A. S. (2001). Validity and reliability of combining three methods to determine ventilatory threshold. Medicine & Science in Sports & Exercise, 33(11), 1841–1848.
#'
#' @examples
#' # TODO write an an example
#'
breakpoint <- function(.data,
                       method = NULL,
                       algorithm_vt1 = NULL,
                       x_vt1 = NULL,
                       y_vt1 = NULL,
                       algorithm_vt2 = NULL,
                       x_vt2 = NULL,
                       y_vt2 = NULL,
                       alpha_linearity = 0.05,
                       vo2 = "vo2",
                       vco2 = "vco2",
                       ve = "ve",
                       time = "time",
                       bps = "both",
                       truncate = TRUE,
                       pos_change = TRUE,
                       # should I add a "front trim" argument? Exercise onset seems to cause
                       # weird issues within fitting lines properly
                       ...) {
    stopifnot(!missing(.data),
              !all(is.null(method), is.null(algorithm_vt1),
                  is.null(x_vt1), is.null(y_vt1)),
              !all(is.null(method), is.null(algorithm_vt2),
                  is.null(x_vt2), is.null(y_vt2)))
    bps <- match.arg(bps, choices = c("vt1", "vt2", "both"), several.ok = FALSE)

    # do I need to restrict the options for x_vt1 etc. to "time", "vo2", "vco2"
    # "speed", "watts", etc? That could serve as the general language to define
    # the relationship graphically. Then I could use the vo2, vco2 arguments to
    # tell the function what column name specifically refers to "vo2", "vco2", etc.
    # or, do I just check to make sure the x_vt1 type arguments can be found
    # in the column names?

    if(!is.null(method)) {
        method <- match.arg(method,
                           choices = c("excess_co2",
                                       "orr",
                                       "v-slope",
                                       "v_slope_simple"))
        if(method == "excess_co2") {
            if(!is.null(x_vt1)) {
                if(x_vt1 != time | x_vt1 != vo2) {
                    warning("x_vt1 is not set to an appropriate value for the excess co2 method. Appropriate x values are time and VO2.")
                }
            }
            if(is.null(x_vt1)) {
                x_vt1 <- time
            }
            if(!is.null(y_vt1)) {
                warning("Overriding calculation of excess co2 for y_vt1")
            }
            if(mean(.data[[vco2]] / .data[[vo2]]) > 2) {
                warning("VO2 and VCO2 columns are unlikely to both be in the same absolute units\nCheck if VO2 column indicated in arguments is absolute and that units match VCO2")
            }
            .data <- .data %>%
                mutate(excess_co2 = .data[[vco2]]^2 / .data[[vo2]] - .data[[vco2]])
            y_vt1 <- "excess_co2"
        }
    }

    # if a specific method was chosen, a function pre-fills algorithms, .x, .y
    # based on the method so long as .x and .y are not already specified

    # orr: VE vs. VO2. This kinda finds RC first if the three-regression model is
    # best, but then only returns vt1. Given that, does this change the truncation
    # decision at RC?
    # V-slope

    algorithm_vt2 <- match.arg(algorithm_vt2,
                  choices = c("dmax", "dmax_mod", "jm",
                              "orr", "v-slope", "simplified_v-slope",
                              "splines"))
    # there's some lactate breakpoints that may be worth adding
    # browser()
    if(bps == "both" | bps == "vt2") {
        params = list(.data = .data,
                      .x = x_vt2,
                      .y = y_vt2,
                      vo2 = vo2,
                      vco2 = vco2,
                      ve = ve,
                      time = time,
                      alpha_linearity = alpha_linearity,
                      bp = "vt2",
                      pos_change = pos_change)
        vt2_out <- switch(algorithm_vt2,
                          "jm" = do.call(what = "jm", args = params),
                          "orr" = do.call(what = "orr", args = params),
                          "v-slope" = do.call(what = "v_slope", args = params),
                          "dmax" = do.call(what = "dmax", args = params))

        if(bps == "vt2") {
            return(vt2_out)
        }

        # truncate if VT2 is found
        if(vt2_out$breakpoint_data$p_val_f < alpha_linearity &
           !is.na(vt2_out$breakpoint_data$p_val_f) & truncate == TRUE) {
            trunc_idx <- which(.data[[time]] == vt2_out$breakpoint_data$time)
            vt1_df <- .data[1:trunc_idx,]
        } else {
            vt1_df <- .data
        }

    } else {
        vt1_df <- .data
    }

    if(bps == "both" | bps == "vt1") {
        params = list(.data = vt1_df,
                      .x = x_vt1,
                      .y = y_vt1,
                      vo2 = vo2,
                      vco2 = vco2,
                      ve = ve,
                      time = time,
                      alpha_linearity = alpha_linearity,
                      bp = "vt1",
                      pos_change = TRUE)
        vt1_out <- switch(algorithm_vt1,
                          "jm" = do.call(what = "jm", args = params),
                          "orr" = do.call(what = "orr", args = params),
                          "v-slope" = do.call(what = "v_slope", args = params),
                          "dmax" = do.call(what = "dmax", args = params))
        if(bps == "vt1") {
            return(vt1_out)
        }
    }

    vt_out <- suppressMessages(full_join(vt1_out$breakpoint_data,
                                         vt2_out$breakpoint_data))
    out <- list(bp_dat = vt_out,
                vt1_dat = vt1_out,
                vt2_dat = vt2_out)
    out
}
