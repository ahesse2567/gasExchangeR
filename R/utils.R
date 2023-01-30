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
#' @return an expression
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
