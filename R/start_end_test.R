start_end_test <- function(.data,
                           intensity_col,
                           start_speed,
                           end_speed = 0) {
    start <- which(lag(diff(.data[[intensity_col]])) == start_speed)
    end <- which(diff(.data[[intensity_col]]) < end_speed)
    if(length(end) == 0) { # test was terminated before clicking recovery on the computer
        end <- nrow(.data) # instead, set end to last data point in file
    }
    .data <- .data[start:end,]
}
