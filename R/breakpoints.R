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
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.
#' @param bps Should the function find the breakpoints for vt1, vt2, or both. Default is \code{both}.
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
        # class(.data) <- append(class(.data), method) # could be useful to
        # use class to pre-select .x and .y variable?
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

    if(bps == "both" | bps == "vt2") {
        bp <- "vt2"
        vt2_out <- switch(algorithm_vt2,
                          "jm" = jm(.data = .data,
                                    .x = x_vt2,
                                    .y = y_vt2,
                                    vo2 = vo2,
                                    vco2 = vco2,
                                    ve = ve,
                                    time = time,
                                    alpha_linearity = alpha_linearity,
                                    bp = bp),
                          "orr" = orr(.data = .data,
                                      .x = x_vt2,
                                      .y = y_vt2,
                                      vo2 = vo2,
                                      vco2 = vco2,
                                      ve = ve,
                                      time = time,
                                      alpha_linearity = alpha_linearity,
                                      bp = bp),
                          "v-slope" = v_slope(.data = .data,
                                               .x = x_vt2,
                                               .y = y_vt2,
                                               vo2 = vo2,
                                               vco2 = vco2,
                                               ve = ve,
                                               time = time,
                                               alpha_linearity = alpha_linearity,
                                               bp = bp))

        if(bps == "vt2") {
            return(vt2_out)
        }
        # browser()
        # truncate if VT2 is found
        if(vt2_out$breakpoint_data$p_val_f < alpha_linearity &
           !is.na(vt2_out$breakpoint_data$p_val_f)) {
            trunc_idx <- which(.data[[time]] == vt2_out$breakpoint_data$time)
            vt1_df <- .data[1:trunc_idx,]
        } else {
            vt1_df <- .data
        }

    } else {
        vt1_df <- .data
    }

    if(bps == "both" | bps == "vt1") {
        bp <- "vt1"
        vt1_out <- switch(algorithm_vt1,
                          "jm" = jm(.data = vt1_df,
                                    .x = x_vt1,
                                    .y = y_vt1,
                                    vo2 = vo2,
                                    vco2 = vco2,
                                    ve = ve,
                                    time = time,
                                    alpha_linearity = alpha_linearity,
                                    bp = bp),
                          "orr" = orr(.data = vt1_df,
                                      .x = x_vt1,
                                      .y = y_vt1,
                                      vo2 = vo2,
                                      vco2 = vco2,
                                      ve = ve,
                                      time = time,
                                      alpha_linearity = alpha_linearity,
                                      bp = bp),
                          "v-slope" = v_slope(.data = vt1_df,
                                              .x = x_vt1,
                                              .y = y_vt1,
                                              vo2 = vo2,
                                              vco2 = vco2,
                                              ve = ve,
                                              time = time,
                                              alpha_linearity = alpha_linearity,
                                              bp = bp))
        if(bps == "vt1") {
            return(vt1_out)
        }
    }

    vt_out <- rbind(vt1_out$breakpoint_data, vt2_out$breakpoint_data)
    out <- list(bp_dat = vt_out,
                vt1_dat = vt1_out[-1],
                vt2_dat = vt2_out[-1])
    out
}
