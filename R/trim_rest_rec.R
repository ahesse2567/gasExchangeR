#' Delete any 'resting/pre-exercise' and 'recovery' portions from a graded exercise test
#'
#' @param .data Unaveraged, breath-by-breath exercise test
#' @param intensity_col Speed or wattage values that help determine start and end of a test
#' @param start_intensity Initial intensity just after 'pre-exercise' data
#' @param end_intensity Intensity that denotes the test has ended or entered the recovery phase
#'
#' @return
#' @export
#'
#' @examples
#' df <- data.frame(time = c(14, 18, 61, 78, 88, 100, 120, 150, 220, 231),
#' speed = c(0, 0, 0, 0, 3, 3, 3, 7.4, 7.4, 0))
#' trim_rest_rec(.data = df, intensity_col = "speed", start_intensity = 3)
trim_rest_rec <- function(.data,
                           intensity_col,
                           start_intensity,
                           end_intensity = 0) {
    start <- which(dplyr::lag(diff(.data[[intensity_col]])) == start_intensity)
    # I think there's an issue here if the speed/grade doesn't repeat
    end <- which(diff(.data[[intensity_col]]) < end_intensity)
    if(length(end) == 0) { # test was terminated before clicking recovery on the computer
        end <- nrow(.data) # instead, set end to last data point in file
    }
    .data <- .data[start:end,]
    .data
}



