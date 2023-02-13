#' Find Other Variables at Ventilatory Thresholds
#'
#' Choose between inverse distance weighted k-nearest neighbors interpolation (`inv-dist`, (default), `dist_x-y`, or `dist_x` to obtain the other values at ventilatory thresholds besides those found directly by the algorithm.
#'
#' @param .data Gas exchange data frame or tibble.
#' @param thr_x X coordinate at the threshold.
#' @param thr_y Y coordinate at the threshold.
#' @param .x Name of the x variable used to find the threshold
#' @param .y Name of the y variable used to find the threshold
#' @param thr_calc_method Name of the threshold calculation method. Default is `inv_dist`. The next recommended method is `dist_x_y`, but `dist_x` is also an option.
#' @param k How many of the nearest data points should be used to calculate the inverse distance weights?
#' @param ... Dot dot dot so this function plays nice with other functions that call it.
#'
#' @details
#' #' Most algorithms that identify ventilatory thresholds are NOT constrained to identify the threshold at an existing data point. Instead, most algorithms locate a point in space that *very* rarely coincides exactly with an existing data point. Given that, what are the values at the threshold *besides* the x and y coordinates you just found? The simplest approach is to find the closest observed data point in one dimension, usually x, and extract all the other variables at that observed data point.
#'
#' The disadvantage of using the x-axis alone to determine the values at the threshold is that the closest observed data point in the x-dimension is often *not* the closest data point to the x-y coordinates of the threshold value. Using values from the closest-x data point may therefore produce misleading values. A better procedure finds the distance (Euclidean) from the threshold point to all observed data points. The nearest data point in the x and y dimensions thus provides the threshold data. The rationale is that an observed data point that is similar to the threshold in two dimensions is likely also more representative in the other dimensions that were not used to find the threshold.
#'
#' However, what if the x-y coordinate of the threshold is close to several observed data points? Wouldn't an average be even better? The default method borrows from geospatial research and uses a version of inverse distanced weighted k-nearest neighbors interpolation. In geospatial research, given a data set of latitude, longitude, and elevation data, one can use this approach to predict the elevation at a new latitude and longitude value. The nearest elevations at known latitude and longitude values contribute more to the predicted elevation because "everything is related to everything else, but near things are more related than distant things" (Tobler, 1970).
#'
#' In a similar fashion, this inverse-distance weight function calculates the Euclidean distance between all observed data points and the x-y threshold coordinate. Then, it selects from the nearest `k` points, defined by the user. Then, a "weight" is assigned proportional to the inverse of the distance of each point divided by the farthest distance of the farthest point ((d / max(d)^-1). The sum of the weights multiplied by their respective observed values become the values for each variable at the threshold.
#'
#' @references
#' Tobler, W. R. (1970). A computer movie simulating urban growth in the Detroit region. Economic geography, 46(sup1), 234-240.
#'
#' @returns A tibble
#' @export
#'
#' @examples
#'
#' # TODO write an example
#'
find_threshold_vals <- function(.data, thr_x, thr_y, .x = .x, .y = .y,
                                thr_calc_method = c("inv_dist",
                                           "dist_x_y",
                                           "dist_x"),
                                k = 5,
                                ...) {
    thr_calc_method <- match.arg(thr_calc_method, several.ok = FALSE)
    class(.data) <- append(class(.data), thr_calc_method)
    UseMethod("find_threshold_vals", .data)
}

#' @keywords internal
find_threshold_vals.inv_dist <- function(.data,
                                         thr_x,
                                         thr_y,
                                         .x = .x,
                                         .y = .y,
                                         thr_calc_method = c("inv_dist",
                                                             "dist_x-y",
                                                             "dist_x"),
                                         k = 5,
                                         ...) {
    thr_calc_method <- match.arg(thr_calc_method, several.ok = FALSE)
    # calculate distance between observed x-y threshold coordinate
    d <- dplyr::bind_rows(tibble::tibble("{.x}" := thr_x, "{.y}" := thr_y),
                   .data %>%
                       dplyr::select(as.name(.x), as.name(.y))) %>%
        dplyr::mutate(dplyr::across(tidyselect::everything(), normalize01)) %>%
        stats::dist()

    cols_to_interpolate <- colnames(.data)[!(colnames(.data) %in%
                                                 c(.x, .y))]
    threshold_data <- numeric(length = length(cols_to_interpolate)) %>%
        rlang::set_names(cols_to_interpolate)
    for(i in seq_along(threshold_data)) {
        if(!is.numeric(.data[[i]])) threshold_data[i] <- NA
        else {
            var_name <- names(threshold_data[i])
            val <- .data %>%
                dplyr::select(as.name(var_name)) %>%
                dplyr::mutate(normalized_dist =
                        d[dist_idx_conv(2:attr(d, "Size"), 1, d)]) %>%
                dplyr::slice_min(normalized_dist, n = k) %>%
                dplyr::mutate(inv_weight = (normalized_dist /
                                      max(normalized_dist))^-1,
                    pct_tot_weight = inv_weight / sum(inv_weight),
                    weighted_contrib = pct_tot_weight *
                        get(var_name)) %>%
                dplyr::summarize("{var_name}" := sum(weighted_contrib)) %>%
                dplyr::pull()
            threshold_data[i] <- val
        }
    }
    threshold_data <- tibble::enframe(threshold_data) %>%
        tidyr::pivot_wider(names_from = name) %>%
        dplyr::mutate("{.x}" := thr_x,
               "{.y}" := thr_y) %>%
        dplyr::relocate(colnames(.data)[colnames(.data) %in%
                                     colnames(.)])
    threshold_data
}

#' @keywords internal
find_threshold_vals.dist_x_y <- function(.data,
                                         thr_x,
                                         thr_y,
                                         .x = .x,
                                         .y = .y,
                                         thr_calc_method = c("inv_dist",
                                                             "dist_x-y",
                                                             "dist_x"),
                                         k = 5,
                                         ...) {
    # calculate distance between observed x-y threshold coordinate
    d <- dplyr::bind_rows(tibble::tibble("{.x}" := thr_x, "{.y}" := thr_y),
                   .data %>%
                       dplyr::select(as.name(.x), as.name(.y))) %>%
        dplyr::mutate(dplyr::across(tidyselect::everything(), normalize01)) %>%
        stats::dist()
    # select row of data closest to the x-y threshold coordiante
    threshold_data <- .data %>%
        dplyr::mutate(normalized_dist =
                d[dist_idx_conv(2:attr(d, "Size"), 1, d)]) %>%
        dplyr::slice_min(normalized_dist, n = 1)

    threshold_data
}

#' @keywords internal
find_threshold_vals.dist_x <- function(.data,
                                       thr_x,
                                       thr_y,
                                       .x = .x,
                                       .y = .y,
                                       thr_calc_method = c("inv_dist",
                                                           "dist_x-y",
                                                           "dist_x"),
                                       k = 5,
                                       ...) {

    threshold_idx <- which.min(abs(.data[[.x]] - thr_x))
    .data[threshold_idx,]
}

#' Convert a 2d index to a 1d index to subset dist objects
#'
#' @keywords internal
#' @noRd
dist_idx_conv <- function (i, j, dist_obj) {
    # convert a 2D to a 1D index
    # i = row (remember this starts at 2), j = column
    if (!inherits(dist_obj, "dist")) stop("please provide a 'dist' object")
    n <- attr(dist_obj, "Size")
    valid <- (i >= 1) & (j >= 1) & (i > j) & (i <= n) & (j <= n)
    k <- (2 * n - j) * (j - 1) / 2 + (i - j)
    k[!valid] <- NA_real_
    k
}
