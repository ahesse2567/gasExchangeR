#' Average gas exchange data from an exercise test
#'
#' This function averages first by either breath, time, or digital filtering
#' If averaging by breath or time averages, it can also perform rolling or bin averages. Furthermore, you can specify if you want a whole or trimmed mean.
#'
#' @param .data Gas exchange data.
#' @param type Choose between breath averages, time averages, or digital filtering.
#' @param subtype Choose rolling, bin, or combo.
#' @param rolling_window How many seconds or breaths to include if rolling.
#' @param bins Bin size of breaths or time.
#' @param align If using a rolling method, how to align the rolling average. Other choices include "left", and "right".
#'
#' @import magrittr
#'
#' @return
#' @export
#'
#' @examples
#'
#' # TODO write an example later b/c I was getting errors earlier despite getting the functions to work when they were in the global environment

avg_exercise_test <- function(.data,
                              type = c("breath", "time", "digital_filter"),
                              subtype = c("rolling", "bin", "combo"),
                              rolling_window = 15,
                              bins = 15,
                              align = "center") {
    stopifnot(!missing(.data),
              !missing(type),
              !missing(subtype),
              rolling_window >= 1,
              bins >= 1,
              rolling_window %% 1 == 0,
              bins %% 1 == 0)

    type <- match.arg(type)
    class(.data) <- append(class(.data), type)
    UseMethod("avg_exercise_test", .data)
}

#' @export
avg_exercise_test.breath <- function(.data,
                                     type = c("breath", "time", "digital_filter"),
                                     subtype = c("rolling", "bin", "combo"),
                                     rolling_window = 15,
                                     bins = 15,
                                     align = "center") {
    subtype <- match.arg(subtype)
    class(.data) <- append(class(.data), subtype)
    UseMethod("avg_exercise_test.breath", .data)
}

#' @export
avg_exercise_test.breath.rolling <- function(.data,
                                     type = c("breath", "time", "digital_filter"),
                                     subtype = c("rolling", "bin", "combo"),
                                     rolling_window = 15,
                                     bins = 15,
                                     align = "center") {
    # browser()
    # get character columns and their values
    char_cols <- .data[, purrr::map(.data, class) == "character"]
    if(dim(char_cols)[1] > 0 & dim(char_cols)[2] > 0) {
        char_cols <- unique(char_cols)
        # delete char col b/c they don't play well with rollapply(). Add back later
        .data <- .data[,-which(colnames(x) %in% names(char_cols))]
    } else {
        char_cols <- NULL
    }

    data_num <- .data %>% # coerce to numeric b/c time may not be of another class
        dplyr::mutate(dplyr::across(where(purrr::negate(is.character)),
                                    as.numeric))

    roll_i <- tibble::tibble()
    for(i in 1:(rolling_window-1)) {
        # calc rolling average for each row including up to the ith row
        # this is how breeze does this
        temp <- purrr::map(data_num, function(data_num) {sum(data_num[1:i])/i})
        roll_i <- dplyr::bind_rows(roll_i, temp)
    }

    align <- match.arg(align, choices = c("left", "right", "center"))

    out <- data_num %>%
        zoo::rollmean(x = ., k = rolling_window, align = align) %>%
        dplyr::as_tibble()

    out <- rbind(roll_i, out)
    out <- dplyr::bind_cols(char_cols, out)
    out
}

#' @export
avg_exercise_test.breath.bin <- function(
    .data,
    type = c("breath", "time", "digital_filter"),
    subtype = c("rolling", "bin", "combo"),
    rolling_window = 15,
    bins = 15,
    align = "center") {

    # browser()
    # get character columns and their values
    char_cols <- .data[, purrr::map(.data, class) == "character"]
    if(dim(char_cols)[1] > 0 & dim(char_cols)[2] > 0) {
        char_cols <- unique(char_cols)
        # delete char col b/c they don't play well with rollapply(). Add back later
        .data <- .data[,-which(colnames(x) %in% names(char_cols))]
    } else {
        char_cols <- NULL
    }

    time_col <- stringr::str_which(colnames(.data), "time")

    data_num <- .data %>% # coerce to numeric b/c time may not be of another class
        dplyr::mutate(dplyr::across(where(purrr::negate(is.character)),
                                    as.numeric))

    # because of piping, group_by_at(1, ...) means group by the first column,
    # which is time. This is probablya vulnerability
    out <- data_num %>%
        dplyr::group_by_at(time_col,
                           function(x) round(x / rolling_window) * rolling_window) %>%
        dplyr::summarise_all(mean, na.rm = TRUE)
    # round(x / rolling_window) puts the values into groups. * rolling_window
    #scales it back to the original time values.
    #Currently I don't think this adjusts to 0, 15, 30 seconds etc.
    out <- dplyr::bind_cols(char_cols, out)
    out
}
