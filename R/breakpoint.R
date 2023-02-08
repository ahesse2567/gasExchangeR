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
#' It is strongly advised the user selects \emph{absolute VO2 of the same units as VCO2 (e.g. mL/min or L/min for both variables)} rather than relative VO2 (mL/kg/min). Using relative VO2 will likely lead to incorrect results.
#'
#' @param .data Gas exchange data
#' @param method Specify \code{excess_co2}, \code{v-slope}, or \code{vent_eqs} and this fills in the most common variables used for these methods for those not already specified. The `orr` and `v-slope_simple` methods will be added in the future.
#' @param algorithm_vt1 Algorithm to find VT1/GET
#' @param x_vt1 \code{x} variable to use to fine VT1/GET/aerobic threshold
#' @param y_vt1 \code{y} variable to use to fine VT1/GET/aerobic threshold
#' @param algorithm_vt2 Algorithm to find VT2/RC
#' @param x_vt2 \code{x} variable to use to fine VT2/RC/anaerobic threshold
#' @param y_vt2 \code{y} variable to use to fine VT2/RC/anaerobic threshold
#' @param bp Should the function find the breakpoints for VT1, TVT2, or both. Default is \code{both}.
#' @param front_trim_vt1 The number of seconds to remove from the *beginning* of the data frame prior to finding vt1. See `Details` below.
#' @param front_trim_vt2 The number of seconds to remove from the data prior to finding vt2 Default = `60`.
#' @param vo2 The name of the vo2 column in \code{.data}
#' @param vco2 The name of the vco2 column in \code{.data}
#' @param ve The name of the ve column in \code{.data}
#' @param time The name of the time column in \code{.data}
#' @param truncate By default, this function truncates the data frame at VT2 prior to finding VT1. Change truncate to \code{FALSE} to use the entire data frame when searching for VT1.
#' @param alpha_linearity Cutoff value to determine if a piecewise model significantly reduces the residual sums of squares more than a simple linear regression.
#' @param pos_change Do you expect the change in slope to be positive (default) or negative? If a two-line regression significantly reduces the sum square error but the change in slope does not match the expected underlying physiology, the breakpoint will be classified as indeterminate.
#' @param ... Arguments to pass to functions internal to `breakpoint()`
#'
#' @returns A list that contains a data frame with slices of the original data frame at the threshold index. The data frame new columns describing the methods used and if the breakpoint was truly a breakpoint. Depending on the breakpoint algorithm used, `breakpoint()` also returns fitted values, the left and right sides of the piecewise regression, as well as a simple linear regression.
#' @export
#' @md
#'
#' @details
#' There are multiple ways to find a breakpoint. The first is by specifying a commonly used \code{method} and the appropriate \code{algorithm} and \code{x}, and \code{y}, variables will be selected for you. For example, entering 'v-slope' for \code{method} will pre-select 'v-slope' for \code{algorithm}, \code{vo2} for \code{x_vt1}, and \code{vco2} for \code{y_vt1}. If you enter a specific \code{method}, you can override the default \code{algorithm}, \code{x}, and \code{y}.
#'
#' Since there aren't common names for every way to find a breakpoint, you can also leave \code{method} as \code{NULL} and specify the \code{algorithm}, \code{x}, and \code{y}.
#'
#' This function uses the excess CO2 formula as described by Gaskill et al. (2001).
#'
#' For some relationships, previous research (Beaver et al., 1986) and other gas exchange analysis programs (Epistemic Mindworks, 2003) recommend removing time at the beginning of the test due to differences in the time constants of VO2 and VCO2 due to the "capacitive effects of changing tissue CO2 stores..." that "...distorts the VCO2 vs. VO2 curve" (Beaver et al., 1986). This function will automatically set `front_trim_vt1` to `60` seconds if the VT1 variables are VCO2 vs. VO2 or if the algorithm is a piecewise function. This will similarly set `front_trim_vt2` to 60 seconds if the VT2 variables are VE vs. VCO2, or if the algorithm is a piecewise function. The piecewise functions include `jm`, `orr`, `v-slope`, and `spline_bp`.
#'
#' @references
#' Gaskill, S. E., Ruby, B. C., Walker, A. V. A. J., Sanchez, O. A., Serfass, R. C., & Leon, A. S. (2001). Validity and reliability of combining three methods to determine ventilatory threshold. Medicine & Science in Sports & Exercise, 33(11), 1841–1848.
#' Beaver, W. L., Wasserman, K. A. R. L. M. A. N., & Whipp, B. J. (1986). A new method for detecting anaerobic threshold by gas exchange. Journal of applied physiology, 60(6), 2020-2027.
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
                       bp = c("both", "vt1", "vt2"),
                       vo2 = "vo2",
                       vco2 = "vco2",
                       ve = "ve",
                       time = "time",
                       front_trim_vt1 = NULL,
                       front_trim_vt2 = NULL,
                       alpha_linearity = 0.05,
                       truncate = TRUE, # this may be uncessary for deriv models
                       pos_change = TRUE,
                       ...) {

    bp <- match.arg(bp, several.ok = FALSE)
    # fill in NULL values, primarily based on "method" argument
    resolved_inputs <- resolve_inputs(inputs = as.list(environment(), all = TRUE))
    list2env(resolved_inputs, envir = environment()) # add values from resolved_inputs

    # do I need to restrict the options for x_vt1 etc. to "time", "vo2", "vco2"
    # "speed", "watts", etc? That could serve as the general language to define
    # the relationship graphically. Then I could use the vo2, vco2 arguments to
    # tell the function what column name specifically refers to "vo2", "vco2", etc.
    # or, do I just check to make sure the x_vt1 type arguments can be found
    # in the column names?

    # if a specific method was chosen, a function pre-fills algorithms, .x, .y
    # based on the method so long as .x and .y are not already specified

    # orr: VE vs. VO2. This kinda finds RC first if the three-regression model is
    # best, but then only returns vt1. Given that, does this change the truncation
    # decision at RC?
    # V-slope

    # algorithm_vt2 <- match.arg(algorithm_vt2,
    #               choices = c("dmax", "dmax_mod", "jm",
    #                           "orr", "v-slope", "simplified_v-slope",
    #                           "splines"))
    # there's some lactate breakpoints that may be worth adding

    if(bp == "both" | bp == "vt2") {
        resolved_inputs[["bp"]] <- "vt2"
        params = append(resolved_inputs, list(.x = x_vt2, .y = y_vt2))
        vt2_out <- switch(algorithm_vt2,
                          "jm" = do.call(what = "jm", args = params),
                          "orr" = do.call(what = "orr", args = params),
                          "v-slope" = do.call(what = "v_slope", args = params),
                          "dmax" = do.call(what = "dmax", args = params),
                          stop("Invalid `algorithm_vt2` value"))

        if(bp == "vt2") {
            return(vt2_out)
        }

        # truncate if VT2 is found
        if(vt2_out$breakpoint_data$determinant_bp & truncate == TRUE) {
            trunc_idx <- which(.data[[time]] == vt2_out$breakpoint_data$time)
            vt1_df <- .data[1:trunc_idx,]
        } else {
            vt1_df <- .data
        }

    } else {
        vt1_df <- .data
    }

    if(bp == "both" | bp == "vt1") {
        params[[".data"]] <- vt1_df
        params[[".x"]] <- x_vt1
        params[[".y"]] <- y_vt1
        vt1_out <- switch(algorithm_vt1,
                          "jm" = do.call(what = "jm", args = params),
                          "orr" = do.call(what = "orr", args = params),
                          "v-slope" = do.call(what = "v_slope", args = params),
                          "dmax" = do.call(what = "dmax", args = params),
                          stop("Invalid `algorithm_vt1` value"))
        if(bp == "vt1") {
            return(vt1_out)
        }
    }

    vt_out <- suppressMessages(dplyr::full_join(vt1_out$breakpoint_data,
                                         vt2_out$breakpoint_data))
    out <- list(bp_dat = vt_out,
                vt1_dat = vt1_out,
                vt2_dat = vt2_out)
    out
}
