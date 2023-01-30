#' Interpolate graded exercise test data with linear or cubic spline interpolation
#'
#' @param .data Data frame to interpolate
#' @param time_col Time column from data frame to reference for interpolation
#' @param every_s Should the data be interpolated every second? Every 2 seconds?
#' @param method Use linear or cubic spline interpolation
#'
#' @returns A data frame or tibble with numeric values interpolated to the specified time frame and by linear or cubic interpolation.
#' @export
#'
#' @details
#' The time column must be in seconds.
#'
#' @examples
#' time <- sample.int(100, 20)
#' time <- sort(time)
#' y <- rnorm(length(time), mean = mean(time), sd = sd(time))
#' df <- data.frame(time = time, y = y)
#' interpolate(.data = df, time_col = "time", every_s = 2, method = "linear")

interpolate <- function(.data,
                        time_col,
                        method = "linear",
                        every_s = 1) {
    # TODO add stopifnot()
    method <- match.arg(method,
                        choices = c("linear", "cubic"),
                        several.ok = FALSE)
    # turn this into a method later
    # class(.data) <- method
    # UseMethod("interpolate", .data)

    data_num <- dplyr::mutate(.data,
                              dplyr::across(tidyselect::where(purrr::negate(is.character)),
                                            as.numeric))
    # add removal of NA vals to the numeric coercion above?

    per_every <- seq(from = min(as.integer(.data[[time_col]])),
                     to = as.integer(max(.data[[time_col]])), by = every_s)

    # do we need to coerce everything to numeric and save the character cols?

    if(method == "cubic") {
        out <- purrr::map(.x = data_num, .f = function(i) stats::spline(
            x = data_num[[time_col]],
            y = i,
            xout = per_every)$y)
        out <- dplyr::as_tibble(out)
    } else {
        out <- purrr::map(.x = data_num, .f = function(i) stats::approx(
            x = data_num[[time_col]],
            y = i,
            xout = per_every)$y)
        out <- dplyr::as_tibble(out)
    }

    out
}
