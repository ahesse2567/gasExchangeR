#' Find breakpoints in gas exchange data
#'
#' @param .data Gas exchange data
#' @param method Specify \code{excess_co2}, \code{v-slope}, or \code{v-slope_simple} and this fills in the most common variables used for these methods for those not already specified.
#' @param algorithm_vt1 Algorithm to find VT1/GET
#' @param x_vt1 \code{x} variable to use to fine VT1/GET/aerobic threshold
#' @param y_vt1 \code{y} variable to use to fine VT1/GET/aerobic threshold
#' @param algorithm_rc Algorithm to find VT2/RC
#' @param x_rc \code{x} variable to use to fine VT2/RC/anaerobic threshold
#' @param y_rc \code{y} variable to use to fine VT2/RC/anaerobic threshold
#' @param vo2 The name of the vo2 column in \code{.data}
#' @param vco2 The name of the vco2 column in \code{.data}
#' @param ve The name of the ve column in \code{.data}
#' @param time The name of the time column in \code{.data}
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simplier model.
#'
#' @return
#' @export
#'
#' @details
#' There are multiple ways to find a breakpoint. The first is by specifying a commonly used \code{method} and the appropriate \code{algorithm} and \code{x}, and \code{y}, variables will be selected for you. For example, entering 'v-slope' for \code{method} will pre-select 'v-slope' for \code{algorithm}, 'vo2' for \code{x_vt1}, and 'vco2' for \code{y_vt1}. If you enter a specific \code{method}, you can override the default \code{algorithm}, \code{x}, and \code{y}.
#'
#' Since there aren't common names for every way to find a breakpoint, you can also leave \code{method} as \code{NULL} and specify the \code{algorithm}, \code{x}, and \code{y}.
#'
#' @examples
#' # TODO write an an example
#'
#' @references
#' Beaver, W. L., Wasserman, K., & Whipp, B. J. (1986). A new method for detecting anaerobic threshold by gas exchange. Journal of Applied Physiology, 121(6), 2020–2027.
#'
breakpoint <- function(.data,
                       method = NULL,
                       algorithm_vt1 = NULL,
                       x_vt1 = NULL,
                       y_vt1 = NULL,
                       algorithm_rc = NULL,
                       x_rc = NULL,
                       y_rc = NULL,
                       alpha_linearity = 0.05,
                       vo2 = "vo2",
                       vco2 = "vco2",
                       ve = "ve",
                       time = "time") {
    stopifnot(!missing(.data),
              !all(is.null(method), is.null(algorithm_vt1),
                  is.null(x_vt1), is.null(y_vt1)),
              !all(is.null(method), is.null(algorithm_rc),
                  is.null(x_rc), is.null(y_rc)))

    # browser()
    if(!is.null(method)) {
        method <- match.arg(method,
                           choices = c("excess_co2",
                                       "orr",
                                       "v-slope",
                                       "v_slope_simple"))
        # class(.data) <- append(class(.data), method) # could be useful to
        # use class to pre-select .x and .y variable?
    }

    # if a specific method was chosen, a function pre-fills algorithms, .x, .y
    # based on the method so long as .x and .y are not already specified
    # orr: VE vs. VO2. This kinda finds RC first if the three-regression model is
    # best, but then only returns vt1. Given that, does this change the truncation
    # decision at RC?
    # V-slope

    ##############################
    # Find RC
    ##############################

    algorithm_rc <- match.arg(algorithm_rc,
                  choices = c("dmax", "dmax_mod", "jm",
                              "orr", "v-slope", "v_slope_simple",
                              "splines"))
    # there's some lactate breakpoints that may be worth adding

    rc_dat <- bp_algorithm(.data = .data,
                           algorithm = algorithm_rc,
                           .x = x_rc,
                           .y = y_rc,
                           vo2 = vo2,
                           vco2 = vco2,
                           ve = ve,
                           time = time)

    rc_dat

    ##############################
    # Truncate test if RC is found
    ##############################

    ##############################
    # Find AerTh
    ##############################

    ##############################
    # Return values and stats
    ##############################

}
