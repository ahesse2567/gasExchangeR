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
#' @param ... Arguments to pass to functions internal to `breakpoint()`
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param pos_change_vt1 Do you expect the change in slope to be positive (default) or negative? If a two-line regression significantly reduces the sum square error but the change in slope does not match the expected underlying physiology, the breakpoint will be classified as indeterminate.
#' @param pos_change_vt2 Do you expect the change in slope to be positive (default) or negative? If a two-line regression significantly reduces the sum square error but the change in slope does not match the expected underlying physiology, the breakpoint will be classified as indeterminate. *The only time you should set this to `FALSE` is when your y-axis variable is end-tidal carbon dioxide (PetCO2)*.
#' @param pos_slope_after_bp Should the slope after the breakpoint be positive? Default is `TRUE`. This catches cases when the percent change in slope is positive, but the second slope is still negative. Change to `FALSE` when PetCO2 is the y-axis variable.
#' @param plots Should this function generate plots? Set to `FALSE` to save time.
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
                       truncate = TRUE, # may be unnecessary for deriv algs
                       pos_change_vt1 = TRUE,
                       pos_change_vt2 = TRUE,
                       pos_slope_after_bp = TRUE,
                       ordering = c("by_x", "time"),
                       plots = TRUE,
                       ...) {

    bp <- match.arg(bp, several.ok = FALSE)
    ordering <- match.arg(ordering, several.ok = FALSE)
    # fill in NULL values, primarily based on "method" argument
    resolved_inputs <- resolve_inputs(inputs = as.list(environment(),
                                                       all = TRUE))
    list2env(resolved_inputs, envir = environment()) # add values from resolved_inputs

    if(bp == "both" | bp == "vt2") {
        resolved_inputs[["bp"]] <- "vt2"
        params = append(resolved_inputs, c(list(.x = x_vt2, .y = y_vt2,
                                                pos_change = pos_change_vt2),
                                           list(...)))
        vt2_out <- switch(algorithm_vt2,
                          "jm" = do.call("jm", params),
                          "orr" = do.call("orr", params),
                          "v-slope" = do.call("v_slope", params),
                          "dmax" = do.call("dmax", params),
                          "spline_bp" = do.call("spline_bp", params),
                          "d1_crossing" = do.call("d1_crossing", params),
                          "d2_inflection" = do.call("d2_inflection", params),
                          "d2_poly_reg_maxima" = do.call(
                              "d2_poly_reg_maxima", params),
                          "d2_reg_spline_maxima" = do.call(
                              "d2_reg_spline_maxima", params),
                          stop("Invalid `algorithm_vt2` value"))
        if(bp == "vt2") {
            return(vt2_out)
        }

        # check that bp was determinant by using the est_ci value
        det_bp <- vt2_out$breakpoint_data %>%
            dplyr::filter(est_ci == "estimate") %>%
            dplyr::select(determinant_bp) %>%
            dplyr::pull()

        # truncate if VT2 is found
        if(det_bp & truncate == TRUE) {
            trunc_val <- vt2_out$breakpoint_data %>%
                dplyr::filter(est_ci == "estimate") %>%
                dplyr::select(x_vt2) %>%
                dplyr::pull()

            trunc_idx <- which.min(abs(.data[[x_vt2]] - trunc_val))
            vt1_df <- .data[1:trunc_idx,]
        } else {
            vt1_df <- .data
        }

    } else {
        vt1_df <- .data
    }

    if(bp == "both" | bp == "vt1") {
        # should we pass the full df for plotting purposes?
        params[[".data"]] <- vt1_df
        params[[".x"]] <- x_vt1
        params[[".y"]] <- y_vt1
        params[["bp"]] <- "vt1"
        params[["pos_change"]] <- pos_change_vt1
        vt1_out <- switch(algorithm_vt1,
                          "jm" = do.call("jm", params),
                          "orr" = do.call("orr", params),
                          "v-slope" = do.call("v_slope", params),
                          "dmax" = do.call("dmax", params),
                          "spline_bp" = do.call("spline_bp", params),
                          "d1_crossing" = do.call("d1_crossing", params),
                          "d2_inflection" = do.call("d2_inflection", params),
                          "d2_poly_reg_maxima" = do.call("d2_poly_reg_maxima",
                                                         params),
                          "d2_reg_spline_maxima" = do.call("d2_reg_spline_maxima",
                                                           params),
                          stop("Invalid `algorithm_vt1` value"))
        if(bp == "vt1") {
            return(vt1_out)
        }
    }
    # combine breakpoint data
    vt_out <- suppressMessages(dplyr::full_join(vt1_out$breakpoint_data,
                                         vt2_out$breakpoint_data)) %>%
        dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp)

    if(plots) {
        # generate plot showing both thresholds
        vt1_plot <- ggplot2::ggplot(data = .data,
                                    ggplot2::aes(x = .data[[x_vt1]], .data[[y_vt1]])) +
            ggplot2::geom_point() +
            ggplot2::theme_minimal()
        vt1_plot <- add_threshold_lines(vt1_plot, x_var = x_vt1,
                                        vt1_out, vt2_out)

        vt2_plot <- ggplot2::ggplot(data = .data,
                                    ggplot2:: aes(x = .data[[x_vt2]], .data[[y_vt2]])) +
            ggplot2::geom_point() +
            ggplot2::theme_minimal()
        vt2_plot <- add_threshold_lines(plt = vt2_plot, x_var = x_vt2,
                                        vt1_out, vt2_out)

        bp_plots <- list("vt1_plot" = vt1_plot, "vt2_plot" = vt2_plot)
    } else {
        bp_plots <- NULL
    }

    out <- list(bp_dat = vt_out,
                bp_plots = bp_plots,
                vt1_dat = vt1_out,
                vt2_dat = vt2_out)
    out
}
