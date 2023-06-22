#' Identify and remove outliers from exercise testing data
#'
#' The current implementation only offers rolling-breath averages.
#'
#' @param .data Gas exchange data frame or tibble.
#' @param outlier_cols Which columns should be used to find outliers. Default is \code{vo2} based on prevalent methodology.
#' @param time Name of the time column. Only used if plotting outliers
#' @param sd_lim How wide should the outlier interval be? Lamarra et al. (1987) originally suggested 3 standard deviations, but 4 is slightly more common.
#' @param width How many values should be considered when calculating the rolling window? Odd numbers are recommended.
#' @param mos Measure of center. Default is the \code{mean}.
#' @param align Should the rolling window be center, left, or right aligned? Default is \code{center}.
#' @param use_global_sd Calculate the rolling window using the global or local standard deviation? Local standard deviations can fluctuate in width considerably, so default is \code{TRUE}.
#' @param global_sd_mos Measure of center used to calculate the global standard deviation, if using. Default is \code{median}.
#' @param exclude_test_val Should the potential outlier be used or excluded from the rolling window? Default is to exclude (\code{TRUE}), because if an outlier affects the rolling window, it may not be appropriately flagged as an outlier.
#' @param remove_outliers Should the function return a data frame with outliers removed or merely noted in a new column?
#' @param max_passes By default, this function repeats itself (\code{max_passes = Inf}) until all values fit inside the rollowing window. However, users may wish to allow only a certain number of passes.
#' @param plot_outliers Plot outliers every iteration of the filter? Default is \code{FALSE}.
#' @param print_outliers Print outliers removed after running function? Default is `FALSE`.
#'
#' @returns A data frame or tibble with the rows containing outliers removed or retained according to the `remove_outliers` parameter. If outliers are removed, then the data frame returns a new column named `outliers`.
#'
#' @export
#'
#' @references
#' Lamarra, N., Whipp, B. J., Ward, S. A., & Wasserman, K. (1987). Effect of interbreath fluctuations on characterizing exercise gas exchange kinetics. Journal of Applied Physiology, 62(5), 2003-2012.
#'
#' @examples
#'
#' # TODO add an example
#'
ventilatory_outliers <- function(.data,
                                 outlier_cols = "vo2",
                                 time = "time",
                                 sd_lim = 3, # add option for conf level later
                                 width = 5,
                                 mos = "mean",
                                 align = "center",
                                 use_global_sd = TRUE,
                                 global_sd_mos = "median",
                                 exclude_test_val = TRUE,
                                 remove_outliers = TRUE,
                                 max_passes = Inf,
                                 print_outliers = FALSE,
                                 plot_outliers = FALSE) {
    # confirm good user input
    mos <- match.arg(mos, choices = c("mean", "median"))
    align <- match.arg(align, choices = c("right", "left", "center"))
    global_sd_mos <- match.arg(global_sd_mos, choices = c("mean", "median"))

    # create key column for flagging all outliers later
    .data <- .data %>%
        dplyr::mutate(key = paste0(.data[[outlier_cols]],
                            "-",
                            1:length(.data[[outlier_cols]])))

    # set proper func if excluding the value to consider
    func <- dplyr::if_else(exclude_test_val & align == "center",
                           "roll_func_exclude_middle", mos)

    # create data frame copy
    copy_df <- .data
    any_outliers <- TRUE
    # if(is.null(max_passes))
    n_passes <- 0

    # TODO do I use lead/lag to allow exclude_test_val to work with
    # right and left-aligned rolling averages?

    while(any_outliers & n_passes < max_passes) {
        rolling_mos <- zoo::rollapply(data = copy_df[[outlier_cols]],
                                 width = width,
                                 FUN = get(func), # need get() to avoid error
                                 f = mos,
                                 fill = NA,
                                 align = align)
        rolling_sd <- zoo::rollapply(data = copy_df[[outlier_cols]],
                                width = width,
                                FUN = stats::sd,
                                fill = NA,
                                align = align)
        if(use_global_sd) {
            rolling_sd <- do.call(global_sd_mos,
                                  args = list(x = rolling_sd, na.rm = TRUE))
        }

        lower_lim <- rolling_mos - rolling_sd * sd_lim
        upper_lim <- rolling_mos + rolling_sd * sd_lim

        # lead or lag if exclude_test_val == TRUE with alignment consideration
        if (exclude_test_val & align == "right") {
            # compares current to preceding breaths
            outlier_idx <- copy_df[[outlier_cols]] < dplyr::lag(lower_lim) |
                copy_df[[outlier_cols]] > dplyr::lag(upper_lim)
        } else if (exclude_test_val & align == "left") {
            # compares against subsequent breaths...I can't really think of why but I guess you can
            outlier_idx <- copy_df[[outlier_cols]] < dplyr::lead(lower_lim) | copy_df[[outlier_cols]] > upper_lim
        } else { # no lead or lag necessary.
            outlier_idx <- copy_df[[outlier_cols]] < lower_lim |
                copy_df[[outlier_cols]] > upper_lim
        }

        n_reassign <- dplyr::if_else(align == "center", floor(width / 2), width - 1)

        # TODO Should we try to use a method to determine if the head
        # and tail of the values are outliers? E.g., use the nearby values
        # to see if they are outliers? However, we need a way to keep the test
        # values from being used to calculate the mean and sd

        offset <- dplyr::if_else(exclude_test_val & align != "center",
                                 1, 0)
        if(align == "center") {
            idx <- c(1:n_reassign,
                     (length(outlier_idx) - n_reassign):length(outlier_idx))
            outlier_idx[idx] <- FALSE
        } else if (align == "left") {
            idx <- (length(outlier_idx) - n_reassign - offset):length(outlier_idx)
            outlier_idx[idx] <- FALSE
        } else if (align == "right") {
            idx <- 1:(n_reassign + offset)
            outlier_idx[idx] <- FALSE
        }

        copy_df <- copy_df %>%
            dplyr::mutate(outlier = outlier_idx)

        if(plot_outliers) {
            outlier_plot <- ggplot2::ggplot(data = copy_df, ggplot2::aes(x = .data[[time]],
                                                       y = .data[[outlier_cols]])) +
                ggplot2::geom_point(ggplot2::aes(color = outlier)) +
                ggplot2::geom_line(ggplot2::aes(y = lower_lim), linetype = "dashed") +
                ggplot2::geom_line(ggplot2::aes(y = upper_lim), linetype = "dashed") +
                ggplot2::ggtitle(glue::glue("Pass {n_passes + 1}, {length(which(outlier_idx))} outliers removed.")) +
                ggplot2::theme_minimal()
            plot(outlier_plot)
        }

        copy_df <- copy_df %>%
            dplyr::filter(!outlier)

        n_passes <- n_passes + 1
        if(all(!outlier_idx)) any_outliers <- FALSE

    }

    # check which keys are still in the original .data
    outliers <- dplyr::if_else(.data$key %in% copy_df$key, FALSE, TRUE)
    .data <- .data %>%
        dplyr::mutate(outlier = outliers)

    if(remove_outliers) {
        .data <- .data %>%
            dplyr::filter(!outlier) %>%
            dplyr::select(-outlier)
    }

    if(any(outliers) & print_outliers) {
        event <- dplyr::if_else(remove_outliers, "removed", "detected")
        outs <- paste(which(outliers), collapse = ", ")
        print(glue::glue("{length(which(outliers))} outliers {event} at indicies {outs}"))
    }

    .data %>%
        dplyr::select(-key)
}

#' @keywords internal
roll_func_exclude_middle <- function(x, f = "mean") {
    # browser()
    # change to more general form that takes "align" variable?
    # if window is even, return average b/c there is no middle value
    if(length(x) %% 2 == 0L) {return(do.call(what = f, args = list(x = x)))}
    # -ceiling(0.5*length(x)) identifies and excludes the middle value
    if(length(x) %% 2 == 1L) {
        return(do.call(what = f,
                       args = list(x = x[-ceiling(0.5*length(x))])))
    }
}
