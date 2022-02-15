#' Find breakpoints in gas exchange data
#'
#' @param .data Gas exchange data
#' @param method Specify \code{excess_co2}, \code{v-slope}, or \code{v-slope_simple} and this fills in the most common variables used for these methods for those not already specified.
#' @param algorithm_aert Algorithm to find VT1/GET
#' @param x_aert
#' @param y_aert
#' @param algorithm_rc Algorithm to find VT2/RC
#' @param x_rc
#' @param y_rc
#' @param vo2_col
#' @param cco2_col
#' @param ve_col
#' @param time_col
#'
#' @return
#' @export
#'
#' @details
#' There are multiple ways to find a breakpoint. The first is by specifying a commonly used \code{method} and the appropriate \code{algorithm} and \code{.x}, and \code{y}, variables will be selected for you. For example, entering 'v-slope' for \code{method} will pre-select 'v-slope' for \code{algorithm}, 'vo2' for .x, and 'vco2' for .y. If you enter a specific \code{method}, you can override the default \code{algorithm}, \code{.x}, and \code{.y}.
#'
#' Since there aren't common names for every way to find a breakpoint, you can also leave \code{method} as \code{NULL} and specify the \code{algorithm}, \code{.x}, and \code{.y}.
#'
#' @examples
#' # TODO write an an example
#'
#' @references
#' Beaver, W. L., Wasserman, K., & Whipp, B. J. (1986). A new method for detecting anaerobic threshold by gas exchange. Journal of Applied Physiology, 121(6), 2020–2027.
#'
breakpoint <- function(.data,
                       method = NULL,
                       algorithm_aert = NULL,
                       x_aert = NULL,
                       y_aert = NULL,
                       algorithm_rc = NULL,
                       x_rc = NULL,
                       y_rc = NULL,
                       vo2_col = "vo2",
                       cco2_col = "vco2",
                       ve_col = "ve",
                       time_col = "time") { # need to add time_col?
    stopifnot(!missing(.data),
              !all(is.null(method), is.null(algorithm_aert),
                  is.null(x_aert), is.null(y_aert)),
              !all(is.null(method), is.null(algorithm_rc),
                  is.null(x_rc), is.null(y_rc)))

    # browser()
    if(!is.null(method)) {
        method <- match.arg(method,
                           choices = c("excess_co2", "v-slope", "v_slope_simple"))
        # class(.data) <- append(class(.data), method) # could be useful to
        # use class to pre-select .x and .y variable?
    }

    # if a specific method was chosen, a function pre-fills algorithms, .x, .y
    # based on the method so long as .x and .y are not already specified

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
