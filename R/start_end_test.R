start_end_test <- function(.data, start_speed) {
    start <- which(lag(diff(.data$speed)) == 3)
    end <- which(diff(.data$speed) < 0)
    if(length(end) == 0) { # test was terminated before clicking recovery on the computer
        end <- nrow(.data) # instead, set end to last data point in file
    }
    .data <- .data[start:end,]
}
