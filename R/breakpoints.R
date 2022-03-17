#' Find breakpoints in gas exchange data
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
                       algorithm_vt2 = NULL,
                       x_vt2 = NULL,
                       y_vt2 = NULL,
                       alpha_linearity = 0.05,
                       vo2 = "vo2",
                       vco2 = "vco2",
                       ve = "ve",
                       time = "time") {
    stopifnot(!missing(.data),
              !all(is.null(method), is.null(algorithm_vt1),
                  is.null(x_vt1), is.null(y_vt1)),
              !all(is.null(method), is.null(algorithm_vt2),
                  is.null(x_vt2), is.null(y_vt2)))

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

    algorithm_vt2 <- match.arg(algorithm_vt2,
                  choices = c("dmax", "dmax_mod", "jm",
                              "orr", "v-slope", "v_slope_simple",
                              "splines"))
    # there's some lactate breakpoints that may be worth adding

    vt2_dat <- bp_algorithm(.data = .data,
                           algorithm = algorithm_vt2,
                           .x = x_vt2,
                           .y = y_vt2,
                           vo2 = vo2,
                           vco2 = vco2,
                           ve = ve,
                           time = time)

    ss <- loop_jm(.data = .data, .x = .x, .y = .y)
    min_ss_idx <- which.min(ss)

    jm_dat <- jm(.data = .data, .x = .x, .y = .y, vo2 = "vo2_abs")

    jm_dat

    ggplot(data = .data, aes(x = .data[[.x]], y = .data[[.y]])) +
        geom_point(alpha = 0.5) +
        theme_bw() +
        geom_point(data = jm_dat$lm_left$model,
                   aes(x = df_left[[.x]], y = df_left[[.y]]),
                   color = "orange", alpha = 0.5) +
        geom_smooth(data = jm_dat$lm_left$model,
                   aes(x = df_left[[.x]], y = df_left[[.y]]),
                   method = "lm",
                   color = "orange", alpha = 0.5) +
        geom_point(data = jm_dat$lm_right$model,
                   aes(x = df_right[[.x]], y = df_right[[.y]]),
                   color = "red", alpha = 0.5) +
        geom_smooth(data = jm_dat$lm_right$model,
               aes(x = df_right[[.x]], y = df_right[[.y]]),
               method = "lm", se = FALSE,
               color = "red", alpha = 0.5) +
        geom_vline(xintercept = jm_dat$breakpoint_data$vco2)

    vt2_dat


    ##############################
    # Truncate test if RC is found
    ##############################
    if(vt2_dat$p_val_F < alpha_linearity) {
        trunc_idx <- which(.data[[time]] == vt2_dat$breakpoint_data$time)
        vt1_df <-
    }

    ##############################
    # Find AerTh
    ##############################

    ##############################
    # Return values and stats
    ##############################

}
