#' Average gas exchange data from an exercise test
#'
#' Averages breath-by-breath gas exchange data using breath-based averages,
#' time-based averages, or digital (Butterworth) filtering. For breath and time
#' methods, choose between rolling or bin averages.
#'
#' To combine bin and rolling averages, call this function twice: first with
#' \code{calc_type = "bin"}, then pass the result into a second call with
#' \code{calc_type = "rolling"} (or vice versa).
#'
#' @param .data Breath-by-breath gas exchange data.
#' @param method Choose between \code{breath} averages, \code{time} averages,
#'   or \code{digital} filtering.
#' @param time_col The name of the column(s) with time. For
#'   \code{method = "digital"}, these columns are excluded from filtering and
#'   passed through unchanged. Accepts a character vector to exclude multiple
#'   columns (e.g. \code{c("time", "clock_time")}).
#' @param calc_type Choose \code{rolling} or \code{bin}. Ignored when
#'   \code{method = "digital"}.
#' @param window How many breaths or seconds to include in each rolling window
#'   or bin. For example, \code{window = 15} with \code{method = "breath"} and
#'   \code{calc_type = "rolling"} gives a 15-breath rolling average;
#'   \code{window = 10} with \code{method = "time"} and \code{calc_type = "bin"}
#'   gives 10-second bin averages.
#' @param align If using a rolling method, how to align the rolling average.
#'   Default is \code{"center"}. Other choices include \code{"left"} and
#'   \code{"right"}.
#' @param mos 'Measure of center'. Choices include \code{"mean"} (default) or
#'   \code{"median"}.
#' @param trim Number of extreme values to remove before averaging. Must be a
#'   non-negative, even integer. For example, \code{trim = 2} with
#'   \code{window = 7} gives a "mid-5-of-7" trimmed average (removes the
#'   highest and lowest values).
#' @param cutoff The cutoff frequency in Hz. Only used by digital filter.
#' @param fs The sampling frequency in Hz. Only used by digital filter.
#' @param order The Butterworth low-pass filter order. Only used by digital
#'   filter.
#'
#' @import magrittr
#'
#' @returns A data frame.
#'
#' @export
#'
#' @references
#' Robergs, R. A., Dwyer, D., & Astorino, T. (2010). Recommendations for
#' improved data processing from expired gas analysis indirect calorimetry.
#' Sports Medicine, 40(2), 95-111.
#' https://doi.org/10.2165/11319670-000000000-00000
#'
#' @examples
#'
#' # Load breath-by-breath graded exercise testing file
#' cpet_bbb <- utils::read.csv(
#'     system.file("extdata", "anton_vo2max_clean.csv", package = "gasExchangeR")
#' )
#'
#' # 10-second (time bin) average
#' cpet_10s_bin <- avg_exercise_test(
#'     cpet_bbb,
#'     method = "time",
#'     calc_type = "bin",
#'     time_col = "time",
#'     window = 10
#' )
#'
#' # 15-breath (breath bin) average
#' cpet_15b_bin <- avg_exercise_test(
#'     cpet_bbb,
#'     method = "breath",
#'     calc_type = "bin",
#'     time_col = "time",
#'     window = 15
#' )
#'
#' # 20-second (rolling time) average
#' cpet_20s_roll <- avg_exercise_test(
#'     cpet_bbb,
#'     method = "time",
#'     calc_type = "rolling",
#'     time_col = "time",
#'     window = 20,
#'     align = "center",
#'     mos = "mean"
#' )
#'
#' # 15-breath (rolling breath) average
#' cpet_15b_roll <- avg_exercise_test(
#'     cpet_bbb,
#'     method = "breath",
#'     calc_type = "rolling",
#'     time_col = "time",
#'     window = 15,
#'     align = "center",
#'     mos = "mean"
#' )
#'
#' # 3rd-order Butterworth low-pass filter using recommendations from
#' # Robergs et al. (2010). time_col excludes time columns from filtering.
#' cpet_dbf <- avg_exercise_test(
#'     cpet_bbb,
#'     method = "digital",
#'     time_col = c("time", "clock_time"),
#'     cutoff = 0.04,
#'     fs = 1,
#'     order = 3
#' )
#'
#' # Middle 5 of 7 breaths per MGC Diagnostics (trimmed, rolling breath
#' # average). trim = 2 removes the highest and lowest values before averaging.
#' cpet_m5o7 <- avg_exercise_test(
#'     cpet_bbb,
#'     method = "breath",
#'     calc_type = "rolling",
#'     time_col = "time",
#'     window = 7,
#'     trim = 2
#' )
#'
avg_exercise_test <- function(.data,
                              method = "breath",
                              calc_type = "rolling",
                              time_col = "time",
                              window = 15,
                              align = "center",
                              mos = "mean",
                              trim = 0,
                              cutoff = 0.04,
                              fs = 1,
                              order = 3) {
    rlang::check_required(.data)
    if (window < 1 || window %% 1 != 0) {
        rlang::abort("`window` must be a positive integer.")
    }
    if (trim < 0 || trim %% 2 != 0) {
        rlang::abort("`trim` must be a non-negative, even integer.")
    }

    method <- rlang::arg_match(method, values = c("breath", "time", "digital"))

    switch(method,
        breath  = avg_breath(.data, calc_type, time_col, window,
                             align, mos, trim),
        time    = avg_time(.data, calc_type, time_col, window,
                           align, mos, trim),
        digital = avg_digital(.data, time_col, cutoff, fs, order)
    )
}

# -- Helpers -------------------------------------------------------------------

#' Separate character columns from numeric data before aggregation.
#'
#' Character columns break rollapply/summarise_all. This helper removes them
#' from .data and returns constant-valued character columns (e.g. subject_id)
#' as a one-row data frame to reattach after aggregation. Non-constant character
#' columns (e.g. time columns read in as characters) are dropped entirely.
#'
#' @returns A list with `data` (character columns removed) and `char_cols`
#'   (one-row data frame of constant character columns, or NULL).
#' @keywords internal
#' @noRd
separate_char_cols <- function(.data) {
    all_chars <- .data[, purrr::map(.data, class) == "character", drop = FALSE]
    if (ncol(all_chars) > 0) {
        constant <- vapply(all_chars,
                           function(x) length(unique(x)) == 1,
                           logical(1))
        char_cols <- all_chars[1, constant, drop = FALSE]
        .data <- .data[, !colnames(.data) %in% colnames(all_chars),
                       drop = FALSE]
    } else {
        char_cols <- NULL
    }
    list(data = .data, char_cols = char_cols)
}

#' Coerce all non-character columns to numeric.
#' @keywords internal
#' @noRd
coerce_numeric <- function(.data) {
    .data %>%
        dplyr::mutate(
            dplyr::across(
                tidyselect::where(purrr::negate(is.character)),
                as.numeric))
}

# -- Method implementations ---------------------------------------------------

#' @keywords internal
#' @noRd
avg_breath <- function(.data, calc_type, time_col, window,
                       align, mos, trim) {
    calc_type <- rlang::arg_match(calc_type, values = c("rolling", "bin"))

    separated <- separate_char_cols(.data)
    char_cols <- separated$char_cols
    data_num <- coerce_numeric(separated$data)

    switch(calc_type,
        rolling = {
            align <- rlang::arg_match(align, values = c("left", "right", "center"))
            mos <- rlang::arg_match(mos, values = c("mean", "median"))

            out <- data_num %>%
                zoo::rollapply(data = .,
                               width = window,
                               align = align,
                               FUN = mos,
                               trim = trim / window / 2) %>%
                dplyr::as_tibble()

            dplyr::bind_cols(char_cols, out)
        },
        bin = {
            out <- data_num %>%
                dplyr::mutate(bin = (1:nrow(.) - 1) %/% window) %>%
                dplyr::group_by(bin) %>%
                dplyr::summarize_all(.funs = list(mos),
                                     na.rm = TRUE,
                                     trim = trim / window / 2) %>%
                dplyr::select(-bin)

            dplyr::bind_cols(char_cols, out)
        }
    )
}

#' @keywords internal
#' @noRd
avg_time <- function(.data, calc_type, time_col, window,
                     align, mos, trim) {
    calc_type <- rlang::arg_match(calc_type, values = c("rolling", "bin"))

    separated <- separate_char_cols(.data)
    char_cols <- separated$char_cols
    data_num <- coerce_numeric(separated$data)

    switch(calc_type,
        rolling = {
            align <- rlang::arg_match(align, values = c("left", "right", "center"))
            mos <- rlang::arg_match(mos, values = c("mean", "median"))

            if (align == "center") {
                a <- window / 2
                b <- window / 2
            } else if (align == "right") {
                a <- window
                b <- 0
            } else {
                a <- 0
                b <- window
            }

            data_num %>%
                dplyr::mutate(
                    dplyr::across(
                        tidyselect::everything(),
                        ~ slider::slide_index_dbl(
                            .,
                            .i = time,
                            .f = mos,
                            na.rm = TRUE,
                            trim = trim / window / 2,
                            .before = b,
                            .after = a,
                            .complete = FALSE
                        ))) %>%
                dplyr::filter(dplyr::if_any(tidyselect::everything(),
                                            ~ !is.na(.)))
        },
        bin = {
            out <- data_num %>%
                dplyr::group_by_at(
                    .vars = time_col,
                    function(x) ceiling(x / window) * window) %>%
                dplyr::summarise_all(.funs = mos,
                                     na.rm = TRUE,
                                     trim = trim / window / 2)

            dplyr::bind_cols(char_cols, out)
        }
    )
}

#' @keywords internal
#' @noRd
avg_digital <- function(.data, time_col, cutoff, fs, order) {
    exclude <- colnames(.data) %in% time_col
    time_df <- .data[, exclude, drop = FALSE]
    .data <- .data[, !exclude, drop = FALSE]

    separated <- separate_char_cols(.data)
    char_cols <- separated$char_cols
    data_num <- coerce_numeric(separated$data)

    bf <- butter_lowpass(cutoff = cutoff, fs = fs, order = order)

    out <- purrr::map(.x = data_num,
                      .f = function(.x, bf) signal::filter(bf, .x),
                      bf = bf)

    dplyr::bind_cols(time_df, char_cols, out)
}

#' @keywords internal
butter_lowpass <- function(cutoff, fs, order = 3){
    nyq <- 0.5 * fs # nyquist frequency is half the sampling rate (fs) b/c you need
    # at a minimum two data points per wave in order to construct the wave
    normal_cutoff <- cutoff / nyq
    signal::butter(n = order, W = normal_cutoff, type = "low", plane = "z")
}
