#' Finding a breakpoint using Orr's 'bruteforce' algorithm.
#'
#' @param .data Gas exchange data.
#' @param .x The x-axis variable.
#' @param .y the y-axis variable.
#' @param alpha_linearity Significance value to determine if a piecewise model explains significantly reduces the residual sums of squares more than a simpler model.
#' @param bp Is this algorithm being used to find vt1 or vt2?
#' @param vo2 The name of the vo2 column in \code{.data}
#' @param vco2 The name of the vco2 column in \code{.data}
#' @param ve The name of the ve column in \code{.data}
#' @param time The name of the time column in \code{.data}
#' @param pos_change Do you expect the change in slope to be positive (default) or negative? If a two-line regression explains significantly reduces the sum square error but the change in slope does not match the expected underlying physiology, the breakpoint will be classified as indeterminate.
#' @param front_trim_vt1 How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param front_trim_vt2 How much data (in seconds) to remove from the beginning of the test prior to fitting any regressions. The original V-slope paper suggests 1 minute.
#' @param ... Dot dot dot mostly allows this function to work properly if breakpoint() passes arguments that is not strictly needed by this function.
#' @param ordering Prior to fitting any functions, should the data be reordered by the x-axis variable or by time? Default is to use the current x-axis variable and use the time variable to break any ties.
#' @param pos_slope_after_bp Should the slope after the breakpoint be positive? Default is `TRUE`. This catches cases when the percent change in slope is positive, but the second slope is still negative. Change to `FALSE` when PetCO2 is the y-axis variable.
#'
#' @returns A list including slice of the original data frame at the threshold index with new columns `algorithm`, `determinant_bp`, `pct_slope_change`, `f_stat`, and `p_val_f.` The list also includes the fitted values, the left and right sides of the piecewise regression, and a simple linear regression.
#'
#' @importFrom rlang :=
#' @export
#'
#' @references
#' Orr, G. W., Green, H. J., Hughson, R. L., & Bennett, G. W. (1982). A computer linear regression model to determine ventilatory anaerobic threshold. Journal of Applied Physiology Respiratory Environmental and Exercise Physiology, 52(5), 1349–1352. https://doi.org/10.1152/jappl.1982.52.5.1349
#'
#' @examples
#' # TODO write examples
orr <- function(.data,
                .x,
                .y,
                bp,
                ...,
                vo2 = "vo2",
                vco2 = "vco2",
                ve = "ve",
                time = "time",
                ordering = c("by_x", "time"),
                alpha_linearity = 0.05,
                front_trim_vt1 = 60,
                front_trim_vt2 = 60,
                pos_change = TRUE,
                pos_slope_after_bp = TRUE
                ) {

    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
    bp <- match.arg(bp, choices = c("vt1", "vt2"), several.ok = FALSE)
    ordering <- match.arg(ordering, several.ok = FALSE)

    .data <- order_cpet_df(.data, .x = .x , time = time,
                           ordering = ordering)
    plot_df <- .data

    front_trim <- set_front_trim(bp = bp,
                                 front_trim_vt1 = front_trim_vt1,
                                 front_trim_vt2 = front_trim_vt2)
    .data <- .data %>%
        dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))

    ss <- loop_orr(.data = .data, .x = .x, .y = .y)
    min_ss_idx <- which.min(ss)

    df_left <- .data[1:min_ss_idx,]
    df_right <- .data[min_ss_idx:nrow(.data),]

    lm_left <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::as.formula() %>%
        stats::lm(data = df_left)
    lm_right <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::as.formula() %>%
        stats::lm(data = df_right)
    lm_simple <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::as.formula() %>%
        stats::lm(data = .data)

    pw_stats <- piecewise_stats(lm_left, lm_right, lm_simple)
    list2env(pw_stats, envir = environment())

    pct_slope_change <- 100*(lm_right$coefficients[2] - lm_left$coefficients[2]) /
        abs(lm_left$coefficients[2])

    y_hat_left <- tibble::tibble("{.x}" := df_left[[.x]],
                         "{.y}" := lm_left$fitted.values,
                         algorithm = "orr")
    y_hat_right <- tibble::tibble("{.x}" := df_right[[.x]],
                          "{.y}" := lm_right$fitted.values,
                          algorithm = "orr")
    pred <- dplyr::bind_rows(y_hat_left, y_hat_right)

    determinant_bp <- check_if_determinant_bp(p = pf_two,
                                              pct_slope_change = pct_slope_change,
                                              pos_change = pos_change,
                                              pos_slope_after_bp =
                                                  pos_slope_after_bp,
                                              slope_after_bp = coef(lm_right)[2],
                                              alpha = alpha_linearity)

    int_point <- intersection_point(lm_left, lm_right)

    bp_dat <- find_threshold_vals(.data = .data, thr_x = int_point["x"],
                                  thr_y = int_point["y"], .x = .x,
                                  .y = .y, ...)
    bp_dat <- bp_dat %>%
        dplyr::mutate(algorithm = "orr",
                      x_var = .x,
                      y_var = .y,
                      determinant_bp = determinant_bp,
                      bp = bp,
                      pct_slope_change = pct_slope_change,
                      f_stat = f_stat,
                      p_val_f = pf_two) %>%
        dplyr::relocate(bp, algorithm, determinant_bp,
                        pct_slope_change, f_stat, p_val_f)

    bp_plot <- make_piecewise_bp_plot(plot_df, .x, .y, lm_left, lm_right, bp_dat)

    return(list(breakpoint_data = bp_dat,
                fitted_vals = pred,
                lm_left = lm_left,
                lm_right = lm_right,
                lm_simple = lm_simple,
                bp_plot = bp_plot))

    # Should we add the three line regression code later?
}
#'
#' @keywords internal
loop_orr <- function(.data, .x, .y) {
    # browser()
    ss_both <- vector(length = nrow(.data))

    for(i in 1:nrow(.data)) {
        if(i == 1 | i == nrow(.data)) {
            ss_both[i] <- NA
            next
        }
        # split data into left and right halves
        # should these share the same point? Or should, they be different by 1 point?
        # Ekkekakis said WB makes them differ by 1 pt
        # I wonder exactly what the word 'contiguous' means in this context
        df_left <- .data[1:i,]
        df_right <- .data[i:nrow(.data),]
        # NOTE: WB'S implementation of Orr does NOT share a point at i

        # make linear models of the two halves
        lm_left <- paste0(.y, " ~ ", "1 + ", .x) %>%
            stats::as.formula() %>%
            stats::lm(data = df_left)
        lm_right <- paste0(.y, " ~ ", "1 + ", .x) %>%
            stats::as.formula() %>%
            stats::lm(data = df_right)

        ss_left <- sum((lm_left$residuals)^2)
        ss_right <- sum((lm_right$residuals)^2)
        #ss_both <- c(ss_both, ss_left + ss_right)
        ss_both[i] <- (ss_left + ss_right)
    }
    return(ss_both)
}
