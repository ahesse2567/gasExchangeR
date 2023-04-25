##########
# This file contains only internal functions
##########

#' Find real values in a vector that contains complex numbers
#'
#' @return a numeric vector
#' @keywords internal
#' @noRd
find_real <- function(v, threshold = 1e-6) {
    # find real roots by fixing rounding errors
    Re(v)[abs(Im(v)) < threshold]
}

#' Create an expression from a vector of polynomial coefficients
#'
#' @return An expression
#' @keywords internal
#' @noRd
expr_from_coefs <- function(poly_coefs, expr = TRUE) {
    string_expr <- paste("x", seq_along(poly_coefs) - 1, sep = "^")
    string_expr <- paste(string_expr, poly_coefs, sep = " * ")
    string_expr <- paste(string_expr, collapse = " + ")
    if (expr) {
        return(parse(text = string_expr))
    } else {
        return(string_expr)
    }
}

#' Find when the slope of a vector changes sign
#'
#' @returns A vector of sign change indices
#' @keywords internal
#' @noRd
slope_sign_changes <- function(y, change = "both") {
    # return the index of sign changes in a vector
    change <- match.arg(change,
                        choices = c("both", "pos_to_neg", "neg_to_pos"),
                        several.ok = FALSE)
    if(change == "both") {
        return(which(diff(sign(diff(y))) != 0) + 1)
    } else if (change == "pos_to_neg") {
        return(which(diff(sign(diff(y))) < 0) + 1)
    } else {
        return(which(diff(sign(diff(y))) > 0) + 1)
    }
}

#' Turn a function into an expression if that makes your life easier
#'
#' @keywords internal
#' @noRd
expr_to_func <- function(expr) {
    if(is.character(expr)) {
        expr <- parse(text = expr) # parse() turns a character string into an expression
    }
    function(x) {eval(expr, envir = list(x = x))} # return a function
}

#' Normalize numeric values between 0 and 1
#'
#' @keywords internal
#' @noRd
normalize01 <- function(x, ...) {
    return((x- min(x)) /(max(x)-min(x)))
}

#' Set the front_trim according to the ventilatory threshold number and front_trim_vtX arguments
#'
#' @keywords internal
#' @noRd
set_front_trim <- function(bp, front_trim_vt1, front_trim_vt2) {
    if (rlang::is_missing("front_trim_vt1")) front_trim_vt1 <- NULL
    if (rlang::is_missing("front_trim_vt2")) front_trim_vt2 <- NULL
    if(bp == "vt1" & !is.null(front_trim_vt1)) {
        front_trim <- front_trim_vt1
    } else if (bp == "vt2" & !is.null(front_trim_vt2)) {
        front_trim <- front_trim_vt2
    } else {
        front_trim <- 0
    }
    front_trim
}

#' Order data by x-variable and time, or by exclusively the time variable
#'
#' @keywords internal
#' @noRd
order_cpet_df <- function(.data, .x, time, ordering = c("by_x", "time")) {
    ordering <- match.arg(ordering, several.ok = FALSE)
    if(ordering == "by_x") {
        .data <- .data %>%
            dplyr::arrange(get(.x), time)
    } else {
        .data <- .data %>%
            dplyr::arrange(time)
    }
    .data
}

#' Order data by x-variable and time, or by exclusively the time variable
#'
#' @keywords internal
#' @noRd
return_indeterminant_findings <- function(.data,
                                 bp,
                                 algorithm,
                                 .x,
                                 .y,
                                 est_ci = "estimate") {
    # extract char/factor columns with unique values to retain ID
    # and related info. Use plot_df since this is a copy
    non_numeric_df <- .data %>%
        dplyr::select(tidyselect::where(
            function(x) is.character(x) |
                is.factor(x) &
                all(x == x[1]))) %>%
        dplyr::slice(1)

    bp_dat <- tibble::tibble(bp = bp,
                   algorithm = algorithm,
                   x_var = .x,
                   y_var = .y,
                   est_ci = "estimate",
                   determinant_bp = FALSE)
    bp_dat <- dplyr::bind_cols(bp_dat, non_numeric_df)
    bp_dat
}
