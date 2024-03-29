#' Average gas exchange data from an exercise test
#'
#' This function averages first by either breath, time, or digital filtering
#' If averaging by breath or time averages, it can also perform rolling or bin averages. Furthermore, you can specify if you want a whole or trimmed mean.
#'
#' @param .data Breath-by-breath gas exchange data.
#' @param method Choose between \code{breath} averages, \code{time} averages, or \code{digital} filtering.
#' @param time_col The name of the column with time.
#' @param calc_type Choose \code{rolling}, \code{bin}, or \code{bin-roll}.
#' @param roll_window How many seconds or breaths to include if rolling.
#' @param bin_w Bin size of breaths or time.
#' @param align If using a rolling method, how to align the rolling average. Default is \code{"center"} Other choices include \code{"left"}, and \code{"right"}.
#' @param mos 'Measure of center'. Choices include \code{"mean"} (default) or \code{"median"}.
#' @param roll_trim Indicate if you want a trimmed mean. Roll_trim removes a number of data points equal to \code{roll_trim} This is used to emulate MCG's "mid-5-of-7" averaging method. \code{roll_trim} must be a positive, even integer.
#' @param bin_trim See \code{roll_trim}.
#' @param cutoff The cutoff frequency in Hz. Only used by digital filter.
#' @param fs The sampling frequency in Hz. Only used by digital filter.
#' @param order The Butterworth low-pass filter order. Only used by digital filter.
#'
#' @import magrittr
#'
#' @returns A data frame.
#'
#' @export
#'
#' @details
#' If you combine rolling and bin averages with \code{calc_type = bin_roll} it is important to note how \code{roll_window} and \code{bin_w} interact. This first creates bin averages that are evenly divisible by the \code{roll_window}. For example, if your bin_w is 5, the function first computes bin averages every 5 breaths or seconds, depending on the \code{method} parameter. Then, if the roll_window was \code{15}, the rolling average would include 3 points in that average because 15 ÷ 5 = 3.
#'
#' \code{roll_window} must be evenly divisible by \code{bin_w}
#'
#' @examples
#'
#' # TODO write an example
#'
#' # TODO force people to use this function twice if doing bin-roll or rolling-bin?
#'
avg_exercise_test <- function(.data,
                              method = "breath",
                              calc_type = "rolling",
                              time_col = "time",
                              roll_window = 15,
                              bin_w = 15,
                              align = "center",
                              mos = "mean",
                              roll_trim = 0,
                              bin_trim = 0,
                              cutoff = 0.04,
                              fs = 1,
                              order = 3) {
    stopifnot(!missing(.data),
              roll_window >= 1,
              bin_w >= 1,
              roll_window %% 1 == 0,
              bin_w %% 1 == 0,
              roll_trim >= 0 & roll_trim %% 2 == 0,
              bin_trim >= 0 & bin_trim %% 2 == 0)

    method <- match.arg(method, choices = c("breath", "time", "digital"))
    class(.data) <- append(class(.data), method)
    UseMethod("avg_exercise_test", .data)
}

#' @export
avg_exercise_test.breath <- function(.data,
                                     method = "breath",
                                     calc_type = "rolling",
                                     time_col = "time",
                                     roll_window = 15,
                                     bin_w = 15,
                                     align = "center",
                                     mos = "mean",
                                     roll_trim = 0,
                                     bin_trim = 0,
                                     cutoff = 0.04,
                                     fs = 1,
                                     order = 3) {
    calc_type <- match.arg(calc_type, choices = c("rolling", "bin", "bin_roll"))

    # save character cols for later
    # I'm beginning to think this is less important and may cause issues
    char_cols <- .data[, purrr::map(.data, class) == "character"]
    if(dim(char_cols)[1] > 0 & dim(char_cols)[2] > 0) {
        char_cols <- unique(char_cols)
        # delete char col b/c they don't play well with rollapply(). Add back later
        .data <- .data[,-which(colnames(.data) %in% names(char_cols))]
        # should I automatically add them back to the original positions later?
    } else {
        char_cols <- NULL
    }

    data_num <- .data %>% # coerce to numeric b/c time may not be of another class
        dplyr::mutate(dplyr::across(tidyselect::where(purrr::negate(is.character)),
                                    as.numeric))
    if(calc_type == "rolling") {
        # rm comments if you want to exactly replicate how breeze does rolling avgs
        # roll_i <- tibble::tibble()
        # for(i in 1:(roll_window-1)) {
        #     # calc rolling average for each row including up to the ith row
        #     # this is how breeze does this
        #     temp <- purrr::map(data_num, function(data_num) {sum(data_num[1:i])/i})
        #     roll_i <- dplyr::bind_rows(roll_i, temp)
        # }

        align <- match.arg(align, choices = c("left", "right", "center"))
        mos <- match.arg(mos, choices = c("mean", "median"))

        out <- data_num %>%
            zoo::rollapply(data = .,
                           width = roll_window,
                           align = align,
                           FUN = mos,
                           trim = roll_trim / roll_window / 2) %>%
            dplyr::as_tibble()

        # out <- rbind(roll_i, out) # if using breeze rolling
        out <- dplyr::bind_cols(char_cols, out)
        return(out)
    } else if (calc_type == "bin") {
        # because of piping, group_by_at(1, ...) means group by the first column.
        # that's probably okay b/c we need to group by breath. For time averaging
        # we should actually find the time_col
        out <- data_num %>%
            dplyr::mutate(bin = (1:nrow(.) - 1) %/% bin_w) %>%
            dplyr::group_by(bin) %>%
            dplyr::summarize_all(.funs = list(mos),
                          na.rm = TRUE,
                          trim = bin_trim / bin_w / 2) %>%
            dplyr::select(-bin)
        out <- dplyr::bind_cols(char_cols, out)
        return(out)
    } else {
        if(roll_window %% bin_w != 0) {
            stop("roll_window is not evenly divisible by bin_w")
        }
        align <- match.arg(align, choices = c("left", "right", "center"))
        mos <- match.arg(mos, choices = c("mean", "median"))

        block <- data_num %>%
            dplyr::mutate(bin = (1:nrow(.) - 1) %/% bin_w) %>%
            dplyr::group_by(bin) %>%
            dplyr::summarize_all(.funs = list(mos),
                          na.rm = TRUE,
                          trim = bin_trim / bin_w / 2) %>%
            dplyr::select(-bin)
        rolled_block <- block %>%
            zoo::rollapply(data = .,
                           width = roll_window / bin_w,
                           align = align,
                           FUN = mos,
                           trim = roll_trim / roll_window / 2) %>%
            dplyr::as_tibble()
        out <- dplyr::bind_cols(char_cols, rolled_block)
        return(out)
    }
}

#' @export
avg_exercise_test.time <- function(.data,
                                   method = "breath",
                                   calc_type = "rolling",
                                   time_col = "time",
                                   roll_window = 15,
                                   bin_w = 15,
                                   align = "center",
                                   mos = "mean",
                                   roll_trim = 0,
                                   bin_trim = 0,
                                   cutoff = 0.04,
                                   fs = 1,
                                   order = 3) {
    # browser()
    calc_type <- match.arg(calc_type, choices = c("rolling", "bin", "bin_roll"))

    # save character cols for later
    # I'm beginning to think this is less important and may cause issues
    char_cols <- .data[, purrr::map(.data, class) == "character"]
    if(dim(char_cols)[1] > 0 & dim(char_cols)[2] > 0) {
        char_cols <- unique(char_cols)
        # delete char col b/c they don't play well with rollapply(). Add back later
        .data <- .data[,-which(colnames(.data) %in% names(char_cols))]
    } else {
        char_cols <- NULL
    }

    data_num <- .data %>% # coerce to numeric b/c time may not be of another class
        dplyr::mutate(dplyr::across(tidyselect::where(purrr::negate(is.character)),
                                    as.numeric))
    if(calc_type == "rolling") {
        align <- match.arg(align, choices = c("left", "right", "center"))
        mos <- match.arg(mos, choices = c("mean", "median"))

        if(align == "center") {
            a <- roll_window / 2
            b <- roll_window / 2
        } else if (align == "right") {
            a <- roll_window
            b <- 0
        } else if (align == "left") {
            a <- 0
            b <- roll_window
        }
        # see https://www.tidyverse.org/blog/2020/02/slider-0-1-0/#index-sliding
        # for details on irregular time series rolling averages

        out <- data_num %>%
            dplyr::mutate(
                dplyr::across(
                    tidyselect::everything(),
                    ~ slider::slide_index_dbl(
                            .,
                            .i = time,
                            .f = mos,
                            na.rm = TRUE,
                            trim = roll_trim / roll_window / 2,
                            .before = b,
                            .after = a,
                            .complete = FALSE
                        ))) %>%
            dplyr::filter(dplyr::if_any(tidyselect::everything(),
                                        ~ !is.na(.)))

        return(out)

    } else if (calc_type == "bin") {
        out <- data_num %>%
            dplyr::group_by_at(.vars = time_col,
                               function(x) ceiling(x / bin_w) * bin_w) %>%
            dplyr::summarise_all(.funs = mos,
                                 na.rm = TRUE,
                                 trim = bin_trim / bin_w / 2)
        # ceiling(x / roll_window) puts the values into groups. * roll_window
        # scales it back to the original time values.
        out <- dplyr::bind_cols(char_cols, out)
        return(out)
    } else {
        if(roll_window %% bin_w != 0) {
            stop("roll_window is not evenly divisible by bin_w")
        }
        align <- match.arg(align, choices = c("left", "right", "center"))
        mos <- match.arg(mos, choices = c("mean", "median"))

        block <- data_num %>%
            dplyr::group_by_at(.vars = time_col,
                               function(x) ceiling(x / bin_w) * bin_w) %>%
            dplyr::summarise_all(.funs = mos,
                                 na.rm = TRUE,
                                 trim = bin_trim / bin_w / 2)
        rolled_block <- block %>%
            zoo::rollapply(data = .,
                           width = roll_window / bin_w,
                           align = align,
                           FUN = mos,
                           trim = roll_trim / roll_window / 2) %>%
            dplyr::as_tibble()
        out <- dplyr::bind_cols(char_cols, rolled_block)
        return(out)
    }
}

#' @export
avg_exercise_test.digital <- function(.data,
                                      method = "breath",
                                      calc_type = "rolling",
                                      time_col = "time",
                                      roll_window = 15,
                                      bin_w = 15,
                                      align = "center",
                                      mos = "mean",
                                      roll_trim = 0,
                                      bin_trim = 0,
                                      cutoff = 0.04,
                                      fs = 1,
                                      order = 3) {
    # browser()

    # save character cols for later
    # I'm beginning to think this is less important and may cause issues
    char_cols <- .data[, purrr::map(.data, class) == "character"]
    if(dim(char_cols)[1] > 0 & dim(char_cols)[2] > 0) {
        char_cols <- unique(char_cols)
        # delete char col b/c they don't play well with rollapply(). Add back later
        .data <- .data[,-which(colnames(.data) %in% names(char_cols))]
    } else {
        char_cols <- NULL
    }

    data_num <- .data %>% # coerce to numeric b/c time may not be of another class
        dplyr::mutate(dplyr::across(tidyselect::where(purrr::negate(is.character)),
                                    as.numeric))

    bf <- butter_lowpass(cutoff = cutoff, fs = fs, order = order)

    # should this exclude the time column?
    out <- purrr::map(.x = data_num,
                      .f = function(.x, bf) signal::filter(bf, .x),
                      bf = bf)

    out <- dplyr::bind_cols(char_cols, out)
    return(out)
}

#' @keywords internal
butter_lowpass <- function(cutoff, fs, order = 3){
    nyq <- 0.5 * fs # nyquist frequency is half the sampling rate (fs) b/c you need
    # at a minimum two data points per wave in order to construct the wave
    normal_cutoff <- cutoff / nyq
    bf <- signal::butter(n = order, W = normal_cutoff, type = "low", plane = "z")
    bf
}
