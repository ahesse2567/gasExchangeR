interpolate <- function(.data, time_col, every_s = 1, method = c("linear", "cubic")) {
    # TODO add stopifnot()
    method <- match.arg(method, choices = c("linear", "cubic"))
    # turn this into a method later
    # class(.data) <- method
    # UseMethod("interpolate", .data)

    data_num <- .data %>% # coerce to numeric b/c time may not be of another class
        dplyr::mutate(dplyr::across(where(purrr::negate(is.character)),
                                    as.numeric))
    # add removal of NA vals to the numeric coercion above?

    per_every <- seq(from = min(as.integer(.data[[time_col]])),
                     to = as.integer(max(.data[[time_col]])), by = every_s)

    # do we need to coerce everything to numeric and save the character cols?

    if(method == "cubic") {
        out <- purrr::map(.x = data_num, .f = function(i) stats::spline(
            x = data_num[[time_col]],
            y = i,
            xout = per_every)$y) %>%
            dplyr::as_tibble()
    } else {
        out <- purrr::map(.x = data_num, .f = function(i) stats::approx(
            x = data_num[[time_col]],
            y = i,
            xout = per_every)$y) %>%
            dplyr::as_tibble()
    }

    out
}
