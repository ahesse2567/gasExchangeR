#' Find breakpoints in gas exchange data (refactored implementation)
#'
#' A behavior-preserving refactor of [breakpoint()]. The control flow, data
#' handling (truncation at VT2), plotting, and return structure are identical to
#' [breakpoint()]; the only structural change is that the duplicated nine-branch
#' algorithm `switch()` is factored into a single internal helper
#' (`dispatch_algorithm()`), and each threshold builds its own argument list.
#'
#' This is currently an internal, non-exported function kept alongside the
#' original so the two can be compared and validated against each other before
#' deciding which to keep. It is intentionally undocumented (`@noRd`) and does
#' not appear in the package NAMESPACE; call it during development with
#' `devtools::load_all()` and then `breakpoint_refactored(...)`, or as
#' `gasExchangeR:::breakpoint_refactored(...)`.
#'
#' Note: unlike [breakpoint()], `bp = "vt1"` works here. In the original,
#' the VT1-only path references a `params` object that is only created inside
#' the VT2 block, so it errors; building each threshold's arguments locally
#' removes that coupling.
#'
#' @keywords internal
#' @noRd
breakpoint_refactored <- function(.data,
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
                                  truncate = TRUE,
                                  pos_change_vt1 = TRUE,
                                  pos_change_vt2 = TRUE,
                                  pos_slope_after_bp = TRUE,
                                  ordering = c("by_x", "time"),
                                  plots = TRUE,
                                  ...) {

    bp <- match.arg(bp, several.ok = FALSE)
    ordering <- match.arg(ordering, several.ok = FALSE)

    # Fill in NULL values, primarily based on the "method" argument
    resolved_inputs <- resolve_inputs(inputs = as.list(environment(),
                                                       all = TRUE))
    list2env(resolved_inputs, envir = environment()) # add resolved values

    dots <- list(...)

    # Dispatch to the requested breakpoint algorithm. Single source of truth for
    # the algorithm -> function mapping (previously duplicated for VT1 and VT2).
    dispatch_algorithm <- function(algorithm, params) {
        switch(algorithm,
               "jm" = do.call("jm", params),
               "orr" = do.call("orr", params),
               "v-slope" = do.call("v_slope", params),
               "dmax" = do.call("dmax", params),
               "spline_bp" = do.call("spline_bp", params),
               "d1_crossing" = do.call("d1_crossing", params),
               "d2_inflection" = do.call("d2_inflection", params),
               "d2_poly_reg_maxima" = do.call("d2_poly_reg_maxima", params),
               "d2_reg_spline_maxima" = do.call("d2_reg_spline_maxima", params),
               stop("Invalid algorithm value: '", algorithm, "'"))
    }

    # Build the argument list for one threshold and run its algorithm. Each call
    # gets a fresh copy of resolved_inputs so VT1 and VT2 never share/mutate the
    # same params object.
    run_threshold <- function(bp_type, .x, .y, algorithm, pos_change, data) {
        params <- append(resolved_inputs,
                         c(list(.x = .x, .y = .y, pos_change = pos_change),
                           dots))
        params[["bp"]] <- bp_type
        params[[".data"]] <- data
        dispatch_algorithm(algorithm, params)
    }

    vt1_out <- NULL
    vt2_out <- NULL
    vt1_df <- .data # default: search VT1 on the full data (used when bp == "vt1")

    # --- VT2 ---------------------------------------------------------------
    if (bp == "both" | bp == "vt2") {
        vt2_out <- run_threshold("vt2", x_vt2, y_vt2, algorithm_vt2,
                                 pos_change_vt2, .data)
        if (bp == "vt2") {
            return(vt2_out)
        }

        # Check that the breakpoint was determinant via the est_ci value
        det_bp <- vt2_out$breakpoint_data %>%
            dplyr::filter(est_ci == "estimate") %>%
            dplyr::select(determinant_bp) %>%
            dplyr::pull()

        # Truncate the data at VT2 before searching for VT1
        if (det_bp & truncate == TRUE) {
            trunc_val <- vt2_out$breakpoint_data %>%
                dplyr::filter(est_ci == "estimate") %>%
                dplyr::select(x_vt2) %>%
                dplyr::pull()

            trunc_idx <- which.min(abs(.data[[x_vt2]] - trunc_val))
            vt1_df <- .data[1:trunc_idx, ]
        } else {
            vt1_df <- .data
        }
    }

    # --- VT1 ---------------------------------------------------------------
    if (bp == "both" | bp == "vt1") {
        vt1_out <- run_threshold("vt1", x_vt1, y_vt1, algorithm_vt1,
                                 pos_change_vt1, vt1_df)
        if (bp == "vt1") {
            return(vt1_out)
        }
    }

    # --- Combine both thresholds (bp == "both") ----------------------------
    vt_out <- suppressMessages(
        dplyr::full_join(vt1_out$breakpoint_data, vt2_out$breakpoint_data)) %>%
        dplyr::relocate(bp, algorithm, x_var, y_var, determinant_bp)

    if (plots) {
        # Plot showing both thresholds
        vt1_plot <- ggplot2::ggplot(
            data = .data,
            ggplot2::aes(x = .data[[x_vt1]], .data[[y_vt1]])) +
            ggplot2::geom_point() +
            ggplot2::theme_minimal()
        vt1_plot <- add_threshold_lines(vt1_plot, x_var = x_vt1,
                                        vt1_out, vt2_out)

        vt2_plot <- ggplot2::ggplot(
            data = .data,
            ggplot2::aes(x = .data[[x_vt2]], .data[[y_vt2]])) +
            ggplot2::geom_point() +
            ggplot2::theme_minimal()
        vt2_plot <- add_threshold_lines(plt = vt2_plot, x_var = x_vt2,
                                        vt1_out, vt2_out)

        bp_plots <- list("vt1_plot" = vt1_plot, "vt2_plot" = vt2_plot)
    } else {
        bp_plots <- NULL
    }

    list(bp_dat = vt_out,
         bp_plots = bp_plots,
         vt1_dat = vt1_out,
         vt2_dat = vt2_out)
}
