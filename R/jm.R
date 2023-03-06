#' Finding a breakpoint using the Jones-Molitoris algorithm.
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
#'
#' @export
#'
#' @references
#' Jones, R. H., & Molitoris, B. A. (1984). A statistical method for determining the breakpoint of two lines. Analytical Biochemistry, 141(1), 287–290. https://doi.org/10.1016/0003-2697(84)90458-5
#'
#' @examples
#' # TODO write examples
jm <- function(.data,
               .x,
               .y,
               bp,
               ...,
               vo2 = "vo2",
               vco2 = "vco2",
               ve = "ve",
               time = "time",
               alpha_linearity = 0.05,
               front_trim_vt1 = 60,
               front_trim_vt2 = 60,
               pos_change = TRUE,
               pos_slope_after_bp = TRUE,
               ordering = c("by_x", "time")
               ) {
    stopifnot(!any(missing(.data), missing(.x), missing(.y), missing(bp)))
    bp <- match.arg(bp, choices = c("vt1", "vt2"), several.ok = FALSE)

    .data <- order_cpet_df(.data, .x = .x , time = time,
                           ordering = ordering)
    plot_df <- .data

    front_trim <- set_front_trim(bp = bp,
                                 front_trim_vt1 = front_trim_vt1,
                                 front_trim_vt2 = front_trim_vt2)

    .data <- .data %>%
        dplyr::filter(.data[[time]] >= min(.data[[time]] + front_trim))

    ss <- loop_jm(.data = .data, .x = .x, .y = .y)
    min_ss_idx <- which.min(ss)

    df_left <- .data[1:min_ss_idx,] # x < x0
    df_right <- .data[(min_ss_idx):nrow(.data),] # x >= x0

    x0 <- .data[[.x]][min_ss_idx]

    lm_left <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::as.formula() %>%
        stats::lm(data = df_left)

    # get y value we will force lm_right to pass through
    b0_plus_b1x0 <- lm_left$coefficients[1] +
        lm_left$coefficients[2] * x0

    # use I() function to force model through (x0, b0_plus_b1x0)
    lm_right <- paste0("I(", .y, "-b0_plus_b1x0)",
                       " ~ 0 + I(", .x, " - ", x0, ")") %>%
        stats::as.formula() %>%
        stats::lm(data = df_right)

    lm_simple <- paste0(.y, " ~ ", "1 + ", .x) %>%
        stats::as.formula() %>%
        stats::lm(data = .data)

    # check for a significant departure from linearity
    pw_stats <- piecewise_stats(lm_left, lm_right, lm_simple)
    list2env(pw_stats, envir = environment())

    y_hat_left <- tibble::tibble("{.x}" := df_left[[.x]],
                         "{.y}" := lm_left$fitted.values,
                         algorithm = "jm")
    y_hat_right <- tibble::tibble("{.x}" := df_right[[.x]],
                          "{.y}" := lm_right$fitted.values + b0_plus_b1x0,
                          algorithm = "jm")
    pred <- dplyr::bind_rows(y_hat_left, y_hat_right)
    # need to use lm_right$coefficients[1] instead of [2] b/c there is only
    # a slope given we forced the line through the origin to satisfy the algorithm
    pct_slope_change <- 100*(lm_right$coefficients[1] - lm_left$coefficients[2]) /
        abs(lm_left$coefficients[2])

    determinant_bp <- check_if_determinant_bp(p = pf_two,
                                              pct_slope_change = pct_slope_change,
                                              pos_change = pos_change,
                                              pos_slope_after_bp =
                                                  pos_slope_after_bp,
                                              slope_after_bp = coef(lm_right)[1],
                                              alpha = alpha_linearity)

    bp_dat <- find_threshold_vals(.data = .data, thr_x = x0,
                                  thr_y = b0_plus_b1x0, .x = .x,
                                  .y = .y, ...)
    bp_dat <- bp_dat %>%
        dplyr::mutate(algorithm = "jm",
                      x_var = .x,
                      y_var = .y,
                      determinant_bp = determinant_bp,
                      bp = bp,
                      pct_slope_change = pct_slope_change,
                      f_stat = f_stat,
                      p_val_f = pf_two) %>%
        dplyr::relocate(bp, algorithm, determinant_bp,
                        pct_slope_change, f_stat, p_val_f)

    plot_x <- plot_df[[.x]]
    plot_df <- tibble::tibble("{.x}" := plot_x) %>%
        dplyr::mutate(
            y_hat_left = stats::predict(lm_left,
                                        tibble::tibble("{.x}" := plot_x)),
            y_hat_right = stats::predict(lm_right,
                                         tibble::tibble("{.x}" := plot_x)) +
                b0_plus_b1x0)

    bp_plot <- ggplot2::ggplot(data = .data,
                               ggplot2::aes(x = .data[[.x]], y = .data[[.y]])) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::geom_line(data = plot_df,
                           ggplot2::aes(x = get(.x), y = y_hat_left)) +
        ggplot2::geom_line(data = plot_df,
                           ggplot2::aes(x = get(.x), y = y_hat_right)) +
        ggplot2::geom_vline(xintercept = bp_dat[[.x]]) +
        ggplot2::theme_minimal()

    return(list(breakpoint_data = bp_dat,
                fitted_vals = pred,
                lm_left = lm_left,
                lm_right = lm_right,
                lm_simple = lm_simple,
                bp_plot = bp_plot))
}

#' @keywords internal
loop_jm <- function(.data, .x, .y) {

    ss_both <- numeric(length = nrow(.data))
    for(i in 1:nrow(.data)) {
        if(i %in% c(1, nrow(.data), nrow(.data) - 1)) {
            ss_both[i] <- NA
            next
        }
        # split data into left and right halves. x0 = .x at index i
        df_left <- .data[1:i,] # x <= x0
        df_right <- .data[(i + 1):nrow(.data),] # x > x0. i+1 correct?
        # does this mean x0 is actually i+1?

        lm_left <- paste0(.y, " ~ ", "1 + ", .x) %>%
            stats::as.formula() %>%
            stats::lm(data = df_left)
        if(is.na(lm_left$coefficients[2])) {
            # avoids a strange corner case when there are only a few
            # data points (e.g. two) and they are both (nearly) identical
            ss_both[i] <- 0
            next
        }

        # according to the JM algorithm, both the left and right regressions
        # are constrained to pass through (x0, y), where x0 is the breakpoint

        # get y value at point x0 with y = b0 + b1*x0,
        b0_plus_b1x0 <- lm_left$coefficients[1] +
            lm_left$coefficients[2] * .data[[.x]][i]
        x0 <- .data[[.x]][i]
        # big thanks to this Stack Overflow post
        # https://stats.stackexchange.com/questions/12484/constrained-linear-regression-through-a-specified-point
        # the JM formula of the right regression is y = b0 + b1*x0 + b3(x-x0)
        # Therefore, we will subtract x0 from each value of x
        # use I() function to force model through (x0, b0_plus_b1x0)
        # This workaround requires forcing the intercept through zero, or at least using that notation
        lm_right <- paste0("I(", .y, "-b0_plus_b1x0)",
                           " ~ 0 + I(", .x, " - ", x0, ")") %>%
            stats::as.formula() %>%
            stats::lm(data = df_right)

        ss_left <- sum((lm_left$residuals)^2)
        ss_right <- sum((lm_right$residuals)^2)
        ss_both[i] <- (ss_left + ss_right)
    }
    ss_both[which(ss_both == 0)] <- NA
    ss_both
}
