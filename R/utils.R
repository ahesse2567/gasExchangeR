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

# ---------------------------------------------------------------------------
# Exact, closed-form polynomial arithmetic (no symbolic engine).
#
# For f(x) = c[1] + c[2]*x + c[3]*x^2 + ... + c[n+1]*x^n (coefficients in the
# order returned by stats::coef(lm(y ~ poly(x, n, raw = TRUE)))), calculus
# gives the derivative's coefficients directly: f'(x)'s coefficient vector is
# c[-1] * 1:n. This is EXACT (a closed-form rule, not a numerical
# approximation) and needs no symbolic differentiation engine.
#
# Using this instead of Deriv::Deriv() + Ryacas::y_fn("Simplify") for the
# per-bootstrap-replicate hot path in d1_crossing()/d2_inflection() is ~800x
# faster for a single differentiation (Deriv::Deriv alone; Ryacas, an
# external computer algebra system, is slower still) and gives identical
# results to machine precision (~1e-12 to 1e-15 in testing) since both
# methods compute the same exact derivative.
# ---------------------------------------------------------------------------

#' Coefficients of the derivative of a polynomial (exact, no symbolic engine).
#' @keywords internal
#' @noRd
poly_deriv_coefs <- function(coefs) {
  if (length(coefs) <= 1) {
    return(0) # derivative of a constant is 0
  }
  # coefs[-1] removes the intercept term
  # seq_along(coefs[-1]) gets the exponent terms to multiply by when taking
  # the derivative
  coefs[-1] * seq_along(coefs[-1])
}

#' Evaluate a polynomial, given its coefficients, at one or more x values.
#' Vectorized over `x`; `coefs[1]` is the constant term.
#' @keywords internal
#' @noRd
poly_eval_coefs <- function(coefs, x) {
  powers <- outer(x, seq_along(coefs) - 1, "^")
  as.vector(powers %*% coefs)
}

#' Subtract two polynomials' coefficient vectors, padding the shorter with 0.
#' @keywords internal
#' @noRd
poly_subtract_coefs <- function(a, b) {
  len_a <- length(a)
  len_b <- length(b)
  # find the highest number of coefficients, in case this subtracts polynomial
  # coefficients of different orders. E.g., a 5th minus a 6th order.
  n <- max(len_a, len_b)
  # pad 0's at the end of each vector of coefficients, if necessary.
  a <- c(a, rep(0, n - len_a))
  b <- c(b, rep(0, n - len_b))
  a - b
}

#' Find when the slope of a vector changes sign
#'
#' @returns A vector of sign change indices
#' @keywords internal
#' @noRd
slope_sign_changes <- function(y, change = "both") {
  # return the index of sign changes in a vector
  change <- rlang::arg_match(
    change,
    values = c("both", "pos_to_neg", "neg_to_pos"),
    multiple = FALSE
  )
  if (change == "both") {
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
  if (is.character(expr)) {
    expr <- parse(text = expr) # parse() turns a character string into an expression
  }
  function(x) {
    eval(expr, envir = list(x = x))
  } # return a function
}

#' Normalize numeric values between 0 and 1
#'
#' @keywords internal
#' @noRd
normalize01 <- function(x, ...) {
  return((x - min(x)) / (max(x) - min(x)))
}

#' Set the front_trim according to the ventilatory threshold number and front_trim_vtX arguments
#'
#' @keywords internal
#' @noRd
set_front_trim <- function(bp, front_trim_vt1, front_trim_vt2) {
  if (rlang::is_missing("front_trim_vt1")) front_trim_vt1 <- NULL
  if (rlang::is_missing("front_trim_vt2")) front_trim_vt2 <- NULL
  if (bp == "vt1" & !is.null(front_trim_vt1)) {
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
  if (ordering == "by_x") {
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
      function(x) {
        is.character(x) |
          is.factor(x) &
            all(x == x[1])
      }
    )) %>%
    dplyr::slice(1)

  bp_dat <- tibble::tibble(
    bp = bp,
    algorithm = algorithm,
    x_var = .x,
    y_var = .y,
    est_ci = "estimate",
    determinant_bp = FALSE
  )
  bp_dat <- dplyr::bind_cols(bp_dat, non_numeric_df)
  bp_dat
}
